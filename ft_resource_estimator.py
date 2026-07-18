#!/usr/bin/env python3
"""
Future-proof fault-tolerant resource estimator.

Brute-force Pauli lambda (qubitization_blueprint / lambda_ladder) is exact but dies past ~4
orbitals. This module computes the factorization-native 1-norms and full FT cost estimates
that are valid at ANY size, by wrapping the authoritative reference implementation
(openfermion.resource_estimates -- the Lee/Babbush et al. code), rather than re-deriving the
published formulas and risking wrong prefactors.

For a given active space it reports:
  * lambda_SF  -- single-factorization 1-norm
  * lambda_DF  -- double-factorization 1-norm (the lever; ~halves SF at real sizes)
  * DF rank / eigenvector count after truncation at a threshold
  * Toffoli count and logical-qubit count for chemical-accuracy QPE
  * T-count ~ 4 x Toffoli

It accepts either a pyscf mean-field object or raw active-space integrals (h1, eri, ecore,
na, nb) -- the same tensors your run_nbn_sqd_sweep.py already produces from a CIF -- so you can
estimate the FT cost of the exact Nb active space you are sweeping on the near-term side.

Requires: openfermion (with resource_estimates), pyscf, pytest installed.
    pip install openfermion pytest
"""

import os
import sys
import contextlib
import numpy as np
from openfermion.resource_estimates import df, sf
from openfermion.resource_estimates.molecule import (
    pyscf_to_cas, cas_to_pyscf, factorized_ccsd_t,
)


@contextlib.contextmanager
def _suppress():
    """File-descriptor-level silence. pyscf's logger caches sys.stdout, so a Python-level
    redirect doesn't catch its (and the C extension's) output -- dup2 over fds 1/2 does.
    Flush FIRST so already-buffered prints reach the real fd before the swap; otherwise a
    mid-region auto-flush dumps them into /dev/null."""
    sys.stdout.flush()
    sys.stderr.flush()
    devnull = os.open(os.devnull, os.O_WRONLY)
    old1, old2 = os.dup(1), os.dup(2)
    os.dup2(devnull, 1)
    os.dup2(devnull, 2)
    try:
        yield
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(old1, 1)
        os.dup2(old2, 2)
        os.close(devnull)
        os.close(old1)
        os.close(old2)


def estimate_from_mf(mf, cas_orbitals=None, cas_electrons=None,
                     thresh=1e-3, dE=0.001, chi=10, beta=20):
    """Full FT resource estimate for a pyscf mean-field object (optionally an active space)."""
    h1, eri, ecore, na, nb = pyscf_to_cas(mf, cas_orbitals, cas_electrons)
    return _estimate(mf, h1, eri, ecore, na, nb, thresh, dE, chi, beta)


def estimate_from_integrals(h1, eri, ecore, na, nb,
                            thresh=1e-3, dE=0.001, chi=10, beta=20):
    """Full FT resource estimate from active-space integrals (e.g. cas.get_h1eff/get_h2eff).

    eri must be the full 4-index tensor in chemist notation, shape (norb,)*4.
    """
    _, mf = cas_to_pyscf(h1, eri, ecore, na, nb)
    return _estimate(mf, h1, eri, ecore, na, nb, thresh, dE, chi, beta)


def _estimate(mf, h1, eri, ecore, na, nb, thresh, dE, chi, beta):
    n = h1.shape[0]
    _, df_facs, rank, nev = df.factorize(eri, thresh)
    lam_df = df.compute_lambda(mf, df_facs)
    _, sf_facs = sf.factorize(eri, n)
    lam_sf = sf.compute_lambda(mf, sf_facs)
    toff_step, toff_total, qubits = df.compute_cost(
        n, lam_df, dE=dE, L=rank, Lxi=nev, chi=chi, beta=beta, stps=20000)
    return {
        "n_orbitals": n, "n_qubits": 2 * n, "nelec": (na, nb),
        "lambda_SF": float(lam_sf), "lambda_DF": float(lam_df),
        "df_rank": int(rank), "df_eigenvectors": int(nev),
        "toffoli_total": int(toff_total), "t_count_approx": int(4 * toff_total),
        "logical_qubits": int(qubits),
    }


def threshold_sweep(mf, thresholds=(1e-1, 1e-2, 1e-3, 1e-4, 1e-5)):
    """Show how the truncation threshold trades DF rank against lambda_DF."""
    h1, eri, ecore, na, nb = pyscf_to_cas(mf)
    rows = []
    for t in thresholds:
        _, df_facs, rank, nev = df.factorize(eri, t)
        rows.append({"thresh": t, "df_rank": int(rank),
                     "lambda_DF": float(df.compute_lambda(mf, df_facs))})
    return rows


def cross_validate_df(eri, norb):
    """Confirm the hand-rolled double_factorize (df_factorization.py) agrees with the
    reference: same exact rank and an exact reconstruction at full rank."""
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from df_factorization import double_factorize, reconstruct_eri
    leaves, _, my_rank = double_factorize(eri, norb)
    recon_err = np.linalg.norm(reconstruct_eri(leaves, norb) - eri)
    _, _, ref_rank, _ = df.factorize(eri, 1e-12)
    return {"my_full_rank": my_rank, "reference_rank": int(ref_rank),
            "my_reconstruction_err": float(recon_err)}


def accuracy_gate(mf, thresholds=(1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 1e-4),
                  target_mHa=1.0, dE=0.001, chi=10, beta=20, no_triples=False):
    """Pick the DF truncation threshold by measured accuracy, not by guessing.

    For each threshold it rank-reduces the ERIs and measures the CCSD(T) energy error that
    truncation actually incurs (vs the untruncated CCSD(T) reference). It then recommends the
    LOOSEST threshold (smallest DF rank -> lowest FT cost) whose error is within target_mHa,
    and returns the full FT cost at that recommended point.

    Note: the truncation error is generally non-monotonic in the threshold, so measuring is
    the point -- a tighter threshold is not guaranteed to be more accurate.
    """
    h1, eri_full, ecore, na, nb = pyscf_to_cas(mf)
    n = h1.shape[0]
    with _suppress():
        _, _, etot_exact = factorized_ccsd_t(mf, None, no_triples=no_triples)

    rows = []
    for t in thresholds:
        with _suppress():
            eri_rr, df_facs, rank, nev = df.factorize(eri_full, t)
            lam = df.compute_lambda(mf, df_facs)
            _, _, etot_rr = factorized_ccsd_t(mf, eri_rr, no_triples=no_triples)
        err = abs(etot_rr - etot_exact) * 1e3
        rows.append({"thresh": float(t), "df_rank": int(rank), "df_eigenvectors": int(nev),
                     "lambda_DF": float(lam), "ccsd_t_err_mHa": float(err),
                     "passes": bool(err <= target_mHa)})

    passing = [r for r in rows if r["passes"]]
    recommended = min(passing, key=lambda r: r["df_rank"]) if passing else None
    if recommended is not None:
        _, toff_total, qubits = df.compute_cost(
            n, recommended["lambda_DF"], dE=dE, L=recommended["df_rank"],
            Lxi=recommended["df_eigenvectors"], chi=chi, beta=beta, stps=20000)
        recommended = {**recommended, "toffoli_total": int(toff_total),
                       "t_count_approx": int(4 * toff_total), "logical_qubits": int(qubits)}

    return {"exact_ccsd_t": float(etot_exact), "target_mHa": target_mHa,
            "rows": rows, "recommended": recommended}




if __name__ == "__main__":
    from pyscf import gto, scf

    print("=" * 78)
    print("Fault-tolerant resource scaling for N2 across basis sets (DF, chemical accuracy)")
    print("-" * 78)
    print(f"{'basis':>10} {'norb':>5} {'qubits':>7} {'lam_SF':>9} {'lam_DF':>9} "
          f"{'DF_rank':>8} {'Toffoli':>11} {'log_qubits':>11}")
    for basis in ["sto-3g", "6-31g", "cc-pVDZ"]:
        mol = gto.M(atom="N 0 0 0; N 0 0 1.0977", basis=basis)
        mf = scf.RHF(mol)
        mf.verbose = 0
        mf.kernel()
        r = estimate_from_mf(mf)
        print(f"{basis:>10} {r['n_orbitals']:>5} {r['n_qubits']:>7} "
              f"{r['lambda_SF']:>9.1f} {r['lambda_DF']:>9.1f} {r['df_rank']:>8} "
              f"{r['toffoli_total']:>11.3e} {r['logical_qubits']:>11}")

    print("-" * 78)
    # cross-validate the hand-rolled DF against the reference on a small ERI
    mol = gto.M(atom="N 0 0 0; N 0 0 1.10", basis="sto-3g")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    h1, eri, ecore, na, nb = pyscf_to_cas(mf)
    cv = cross_validate_df(eri, h1.shape[0])
    print(f"DF cross-check (my module vs reference): {cv}")

    # the integrals path -- exactly what run_nbn_sqd_sweep.py feeds on the near-term side
    r = estimate_from_integrals(h1, eri, ecore, na, nb)
    print(f"from-integrals path OK: lambda_DF={r['lambda_DF']:.2f}, "
          f"Toffoli={r['toffoli_total']:.3e}, logical_qubits={r['logical_qubits']}")
    print("=" * 78)

    # CCSD(T) accuracy gate: choose the threshold by measured error, not by guessing
    print("\nCCSD(T) accuracy gate for N2 / 6-31g (target 1 mHa)")
    print("-" * 78)
    mol = gto.M(atom="N 0 0 0; N 0 0 1.0977", basis="6-31g")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    gate = accuracy_gate(mf, target_mHa=1.0)
    print(f"{'thresh':>8} {'rank':>5} {'lambda_DF':>10} {'ccsd_t_err_mHa':>15} {'pass':>6}")
    for row in gate["rows"]:
        print(f"{row['thresh']:>8.0e} {row['df_rank']:>5} {row['lambda_DF']:>10.2f} "
              f"{row['ccsd_t_err_mHa']:>15.4f} {str(row['passes']):>6}")
    rec = gate["recommended"]
    print("-" * 78)
    if rec:
        print(f"RECOMMENDED threshold {rec['thresh']:.0e}: rank {rec['df_rank']}, "
              f"lambda_DF {rec['lambda_DF']:.2f}, err {rec['ccsd_t_err_mHa']:.3f} mHa "
              f"-> Toffoli {rec['toffoli_total']:.2e}, logical_qubits {rec['logical_qubits']}")
        print("Loosest threshold that still clears chemical accuracy = lowest FT cost. The error")
        print("is non-monotonic in the threshold, which is exactly why it is measured, not assumed.")
    else:
        print("No threshold met the target; tighten the sweep or relax target_mHa.")
    print("=" * 78)

