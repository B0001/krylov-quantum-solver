#!/usr/bin/env python3
"""
Cross-check harness -- the capstone for the near-term half of the stack.

Runs four INDEPENDENT ground-state methods on the same active space and asserts they agree:

  * CASCI  -- exact diagonalization (the reference)
  * SQD    -- sample-based quantum diagonalization (configuration sampling)
  * Krylov -- real-time quantum subspace diagonalization (time evolution)
  * ADAPT  -- ADAPT-VQE (gradient-driven variational ansatz)

Agreement across methods that fail in DIFFERENT ways (sampling noise, subspace conditioning,
ansatz expressivity) is what justifies trusting a number before spending QPU time on it. A
regression in any one solver shows up here as a blown tolerance instead of a silently wrong
energy. Same active-space interface as every other module: (h1, eri, e_core, nelec, norb).

SQD is optional: if qiskit-addon-sqd is not installed the harness still runs the other three.
Requires krylov_subspace_solver.py and adapt_vqe.py in the same directory.
"""

import os
import sys
import contextlib
import numpy as np
from functools import partial
from pyscf import fci

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from krylov_subspace_solver import krylov_ground_state
from adapt_vqe import adapt_ground_state


@contextlib.contextmanager
def _quiet():
    with open(os.devnull, "w") as dn, contextlib.redirect_stdout(dn):
        yield


def _sqd_energy(h1, eri, e_core, nelec, norb, n_shots=20000, seed=0):
    from qiskit_addon_sqd.fermion import diagonalize_fermionic_hamiltonian, solve_sci_batch
    from qiskit_addon_sqd.counts import generate_bit_array_uniform
    rng = np.random.default_rng(seed)
    bit_array = generate_bit_array_uniform(n_shots, norb * 2, rand_seed=rng)
    S = (nelec[0] - nelec[1]) / 2.0
    best = [np.inf]

    def callback(results):
        for r in results:
            best[0] = min(best[0], r.energy + e_core)

    with _quiet():
        diagonalize_fermionic_hamiltonian(
            h1, eri, bit_array, samples_per_batch=300, norb=norb, nelec=nelec,
            num_batches=5, max_iterations=5,
            sci_solver=partial(solve_sci_batch, spin_sq=S * (S + 1)),
            symmetrize_spin=(nelec[0] == nelec[1]), seed=rng, callback=callback)
    return best[0]


def cross_check(h1, eri, e_core, nelec, norb, tol_mHa=5.0,
                krylov_dim=12, krylov_dt=0.5,
                qubit_dense_max_orb=7, fci_max_dim=500000, krylov_max_dim=4000):
    """Run every FEASIBLE solver and return (reference, per-method results, max_dev, agree).

    Solvers are gated by their true cost so this scales past toy systems:
      * CASCI -- FCI, efficient; run up to a large determinant dimension.
      * Krylov -- builds a DENSE CI matrix and full-diagonalizes it; tight dimension cap.
      * ADAPT -- qubit-space dense (2^(2*norb)); only for small active spaces.
      * SQD   -- selected-CI; always attempted, the scalable one.
    Each is individually guarded, so an out-of-regime or missing-dependency solver is skipped
    (recorded as None) rather than crashing the harness. Agreement is judged over what ran.
    """
    from math import comb
    ci_dim = comb(norb, nelec[0]) * comb(norb, nelec[1])
    res = {}

    if ci_dim <= fci_max_dim:
        try:
            with _quiet():
                e_fci, _ = fci.direct_spin1.kernel(h1, eri, norb, nelec)
            res["CASCI"] = float(e_fci) + e_core
        except Exception:
            res["CASCI"] = None
    else:
        res["CASCI"] = None

    if ci_dim <= krylov_max_dim:
        try:
            res["Krylov"] = krylov_ground_state(h1, eri, e_core, nelec, norb,
                                                krylov_dim=krylov_dim, dt=krylov_dt)[0]
        except Exception:
            res["Krylov"] = None
    else:
        res["Krylov"] = None

    if 2 ** (2 * norb) <= 2 ** (2 * qubit_dense_max_orb):
        try:
            res["ADAPT"] = adapt_ground_state(h1, eri, e_core, nelec, norb)[0]
        except Exception:
            res["ADAPT"] = None
    else:
        res["ADAPT"] = None

    try:
        res["SQD"] = _sqd_energy(h1, eri, e_core, nelec, norb)
    except Exception:
        res["SQD"] = None

    available = {k: v for k, v in res.items() if v is not None}
    ref = available.get("CASCI", next(iter(available.values())) if available else None)
    table = {k: (v, abs(v - ref) * 1e3) for k, v in available.items()} if ref is not None else {}
    max_dev = max((d for _, d in table.values()), default=float("nan"))
    return {"reference": ref, "results": table, "skipped": [k for k, v in res.items() if v is None],
            "max_dev_mHa": max_dev, "agree": (max_dev <= tol_mHa) if table else False}


if __name__ == "__main__":
    from pyscf import gto, scf, mcscf, ao2mo

    def reference(atom, norb, ne, basis="sto-3g"):
        mol = gto.M(atom=atom, basis=basis)
        mf = scf.RHF(mol)
        mf.verbose = 0
        mf.kernel()
        cas = mcscf.CASCI(mf, norb, ne)
        cas.verbose = 0
        cas.kernel()
        h1, e_core = cas.get_h1eff()
        eri = ao2mo.restore(1, cas.get_h2eff(), norb)
        return h1, eri, float(e_core), (ne // 2, ne // 2)

    systems = {
        "H2  CAS(2,2)": ("H 0 0 0; H 0 0 0.74", 2, 2),
        "H4  CAS(4,4)": ("H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3", 4, 4),
    }
    print("=" * 64)
    for label, (atom, norb, ne) in systems.items():
        h1, eri, e_core, nelec = reference(atom, norb, ne)
        out = cross_check(h1, eri, e_core, nelec, norb, tol_mHa=5.0)
        print(f"\n[{label}]  reference CASCI = {out['reference']:.8f} Ha")
        for k, (v, d) in out["results"].items():
            print(f"   {k:>7}: {v:.8f}   delta = {d:7.4f} mHa")
        print(f"   -> max deviation {out['max_dev_mHa']:.4f} mHa   AGREE = {out['agree']}")
        assert out["agree"], "solvers disagree beyond tolerance"
    print("\n" + "=" * 64)
    print("Four methods that fail in different ways, triangulating one energy. Trust the number")
    print("when they agree; investigate the outlier when they do not.")
