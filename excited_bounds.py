#!/usr/bin/env python3
"""
Certified lower bounds on the FIRST EXCITED state of the Krylov solve -- Lehmann's optimal
subspace bounds, and the separator problem that makes self-certification impossible.

`temple_bounds.py` brackets E_0 from below; `certified_gaps.py` brackets the gap using the
*unverified* floor eps_1 = theta_1 - sigma_1 for E_1. This module replaces that floor with a real
lower bound on E_1 -- and then shows, by construction, that the floor can never be made
self-certifying from Krylov data alone.

THE BOUND (Lehmann 1949/1950; Maehly; modern treatment Beattie & Goerisch, Numer. Math. 72, 143,
1995). Choose a separator beta with theta_1 < beta <= E_2, where E_2 is the third eigenvalue of the
HF-REACHABLE sector. On the span of the Ritz vectors with theta_i < beta form

    A1 = W^dag (H - beta) W   (negative definite by that selection),
    A2 = W^dag (H - beta)^2 W,   and solve   A2 y = tau A1 y .

Sorted descending, the d roots beta + tau are lower bounds on the d eigenvalues immediately below
beta:  beta + tau_(1) <= E_1  and  beta + tau_(2) <= E_0  when exactly two eigenvalues lie below
beta. The upper bounds are free: theta_i >= E_i by Poincare separation (Cauchy interlacing), the
same fact `solve_excited` already rests on. At d = 1 the pencil reduces algebraically to Temple's
inequality; at d = 2 it is strictly tighter (measured +0.07 mHa on LiH at M = 4).

NOT Lehmann's *enclosure*, and not what SPEC v1.0 wrote down. v1.0's formula
E_i >= theta_i - sigma_i^2/(beta - theta_i) is the single-vector Temple/Kato inequality applied
per Ritz vector, not Lehmann's multi-dimensional theorem. Both are implemented here
(`temple_lower_1` carries the per-vector number for comparison); the Lehmann pencil is the one the
brackets use, because it is the optimal bound the subspace supports.

THE FINDING (specs/SPEC_excited_state_certification.md v2.0):

  * With an oracle separator the bound is rigorous and it works: containment of E_0 AND E_1 on
    LiH / H4 / N2 CAS(6,6) at every M in 4..24 (zero escapes), plus zero escapes over 300 random
    Hermitian matrices with random subspaces. The E_1 bracket closes 88 -> 0.026 mHa (LiH,
    M = 4 -> 20) and 216 -> 0.002 mHa (H4, M = 6 -> 20) -- but NOT monotonically.
  * Certification is cheap in width, not only in cost: the certified gap interval is 1.2-1.3x the
    (unknowable) raw Ritz gap error on LiH, 2.3-5.4x on H4, 2.6-11.7x on N2. The spec's "< 5x"
    holds on LiH, fails on H4 at M >= 20 and fails by 2.3x on N2.
  * SELF-CERTIFICATION IS FALSE. beta = theta_2 - sigma_2 is an estimate, not a bound: Kato's
    interval says only that SOME eigenvalue lies within sigma_2 of theta_2, never that it is E_2.
    When the subspace has not resolved level 2 the estimate overshoots E_2, three eigenvalues sit
    below beta instead of two, and the pencil's top root -- a true bound on E_2 -- is reported as a
    bound on E_1. It escapes on real molecules, not just in principle: LiH at M = 4, 6, 8 (E_1
    under-bounded by 2.0 mHa) and N2 CAS(6,6) at M = 16 (by 17.4 mHa).
  * The impossibility is structural, not a depth problem. A reachable level with tiny HF overlap is
    invisible to the Krylov data: the constructed witness in `tests/test_excited_bounds_spec.py`
    hides a level of amplitude 1e-4 (population 1e-8, ABOVE this repo's 1e-10 reachability cut) and
    the self-certified bound overshoots E_1 by 0.5 Ha while every residual says the subspace is
    converged. No function of the subspace data can see it. Same class of defect as
    SPEC_subspace_floor_resolvability's ~1e-4-amplitude level and SPEC_reachability_tolerance's
    threshold divergence.
  * Be2 -- the spec's showcase -- is out of scope for a structural reason. At CAS(4,8)/cc-pVDZ its
    lowest excited singlets are a DEGENERATE pi pair (E_1 = E_2 = -29.085033 Ha). A separator needs
    theta_1 < beta <= E_2 while theta_1 >= E_1 = E_2, so no beta exists at any M, in any basis.
    On top of that those states are not HF-reachable by real-time evolution (Krylov rank saturates
    at 3, theta_1 sits 147 mHa above E_1^CASCI at every M <= 14), so the optically bright
    transition the spec asks to certify is not in the subspace at all. The benchmark keeps Be2 as
    a refusal row.

HONEST SCOPE. (1) The precondition is exactly theta_1 < beta <= E_2 over the REACHABLE spectrum
(the Krylov space lies in that sector, so unreachable levels between E_1 and beta are harmless);
given a valid beta the inequality holds for any states, so all the risk lives in beta. (2) A
degenerate E_1 admits no separator -- the bound is vacuous by construction, not by failure.
(3) An oracle beta means an independent lower bound on E_2 (CASCI/DMRG/experiment); self mode is
an estimate and is labelled as one -- do not quote it as a certificate. (4) -inf is a valid,
deliberate refusal (the precondition failed); check `isfinite` before quoting a width.
(5) Exact statevector: the hardware shot cost of <H^2> (a ~lambda^2 Pauli expansion) is not
modelled, and shot noise on H/S is not propagated into the bound. (6) ARITHMETIC FLOOR ~1e-8 Ha:
the Ritz vectors of a thresholded canonical orthogonalisation are orthonormal only to ~1e-9, the
pencil divides that by beta - theta_1, and the default dt comes from ARPACK with a random start
vector -- so L_1 moves by up to ~1e-8 Ha between identical runs. The bound is rigorous in exact
arithmetic; below ~1e-8 Ha this implementation reports noise, and a vanishing beta - theta_1 makes
it worse (the constructed witness in the gate file has beta - theta_1 = 3e-9 and leaks 1e-6 Ha).
(7) These certificates do NOT
plug into `certkit_bridge.py`: certkit v0.2.0's only claim kind is `lambda_min_enclosure` and its
Temple rule takes one vector plus a beta -- a two-vector pencil claim about the SECOND eigenvalue
has no rule there, so emitting one would be unchecked by construction.

REPRODUCTION vs NOVELTY: the bound is textbook (Lehmann 1949; Beattie & Goerisch 1995; Zimmermann
& Mertins, Z. Angew. Math. Mech. 75, 1995) and its numerical use is a well-established literature.
What is new here is only the composition -- Lehmann bounds on a real-time QKSD basis -- and the
negative result that the separator cannot be self-certified from that basis.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np
from scipy.linalg import eigh

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from reachability import reachable_eigenpairs
from temple_bounds import mean_and_variance


@dataclass(frozen=True)
class ExcitedBracket:
    """Certified brackets on E_0 and E_1 at one Krylov dimension. TOTAL energies in Ha.

    ``lower_*`` is -inf whenever the separator precondition fails -- a refusal, and a valid but
    vacuous bound. ``upper_*`` are the Ritz values, which are always finite.
    """
    m: int                    # Krylov dimension
    lower_0: float            # Lehmann lower bound on E_0
    upper_0: float            # theta_0 (Poincare/variational upper bound)
    lower_1: float            # Lehmann lower bound on E_1
    upper_1: float            # theta_1 (Poincare upper bound, Cauchy interlacing)
    temple_lower_1: float     # per-vector Temple bound on E_1 (comparison; <= lower_1)
    beta: float               # separator actually used
    beta_source: str          # "oracle" (caller-supplied E_2 floor) | "self" (theta_2 - sigma_2)
    sigma1: float             # sqrt(variance) of the first-excited Ritz state
    rank: int                 # effective subspace rank kept by the overlap threshold

    @property
    def best_estimate_0(self) -> float:
        return self.upper_0

    @property
    def best_estimate_1(self) -> float:
        return self.upper_1

    @property
    def width_1(self) -> float:
        """Excited-state bracket width U_1 - L_1 (inf when the separator was refused)."""
        return self.upper_1 - self.lower_1

    @property
    def gap_lower(self) -> float:
        """Certified minimum possible gap E_1 - E_0."""
        return self.lower_1 - self.upper_0

    @property
    def gap_upper(self) -> float:
        """Certified maximum possible gap E_1 - E_0."""
        return self.upper_1 - self.lower_0

    @property
    def gap_width(self) -> float:
        return self.gap_upper - self.gap_lower


def lehmann_lower_bounds(H, states: Sequence[np.ndarray], beta: float) -> np.ndarray:
    """Lehmann lower bounds on the eigenvalues immediately below ``beta``, descending.

    ``states``: rows are orthonormal Ritz vectors (``QuantumKrylovSolver.eigenstates``). Only those
    with theta_i < beta enter the pencil -- that selection is what makes ``A1`` negative definite,
    which is Lehmann's hypothesis. Returns ``[bound on E_{nu-1}, bound on E_{nu-2}, ...]`` where nu
    is the number of eigenvalues strictly below beta. **nu is an input the caller must know**; this
    function cannot check it and neither can the subspace (see the module docstring).

    Empty array when no Ritz value lies below beta -- a refusal, never a guess.
    """
    W = np.asarray(states).T                                  # (N, k), columns are |u_i>
    if not np.isfinite(beta) or W.size == 0:
        return np.empty(0)
    HW = H @ W
    theta = np.real(np.einsum("ij,ij->j", W.conj(), HW))
    keep = theta < beta
    if not keep.any():
        return np.empty(0)
    W, HW = W[:, keep], HW[:, keep]
    A1 = W.conj().T @ HW - beta * np.eye(W.shape[1])           # W^dag (H - beta) W
    R = HW - beta * W
    A2 = R.conj().T @ R                                        # W^dag (H - beta)^2 W
    neg = -0.5 * (A1 + A1.conj().T)                            # positive definite by construction
    if float(np.linalg.eigvalsh(neg).min()) <= 0.0:            # non-Ritz input; refuse
        return np.empty(0)
    tau = eigh(0.5 * (A2 + A2.conj().T), neg, eigvals_only=True)   # ascending, positive
    return beta - tau                                          # descending, all < beta


def bracket_from_states(H, states: Sequence[np.ndarray], beta: Optional[float] = None,
                        m: int = 0, offset: float = 0.0) -> ExcitedBracket:
    """Brackets on E_0/E_1 from Ritz vectors of any subspace -- the matrix-level entry point.

    ``beta``: an ELECTRONIC-frame separator satisfying theta_1 < beta <= E_2 (oracle mode), or None
    to estimate it from the same data as theta_2 - sigma_2 (self mode -- an ESTIMATE, gated false
    as a certificate by G4/G4b). ``offset`` is added to every energy on the way out.
    """
    states = list(states)
    th, var = zip(*(mean_and_variance(H, s) for s in states))
    if beta is None:
        beta_e = th[2] - math.sqrt(var[2]) if len(states) > 2 else -np.inf
        source = "self"
    else:
        beta_e = float(beta)
        source = "oracle"
    bounds = lehmann_lower_bounds(H, states, beta_e)
    lo1 = float(bounds[0]) if len(bounds) >= 1 else -np.inf
    lo0 = float(bounds[1]) if len(bounds) >= 2 else -np.inf
    th1 = th[1] if len(states) > 1 else np.inf
    var1 = var[1] if len(states) > 1 else np.inf
    temple1 = th1 - var1 / (beta_e - th1) if beta_e > th1 else -np.inf
    return ExcitedBracket(
        m=m,
        lower_0=lo0 + offset if np.isfinite(lo0) else -np.inf,
        upper_0=th[0] + offset,
        lower_1=lo1 + offset if np.isfinite(lo1) else -np.inf,
        upper_1=th1 + offset if np.isfinite(th1) else np.inf,
        temple_lower_1=temple1 + offset if np.isfinite(temple1) else -np.inf,
        beta=beta_e + offset if np.isfinite(beta_e) else -np.inf,
        beta_source=source,
        sigma1=float(np.sqrt(var1)) if np.isfinite(var1) else np.inf,
        rank=len(states),
    )


def lehmann_excited_brackets(mh: MolecularHamiltonian, krylov_dim: int,
                             beta_oracle: Optional[float] = None,
                             solver: Optional[QuantumKrylovSolver] = None) -> ExcitedBracket:
    """Certified brackets on E_0 and E_1 from a ``krylov_dim``-dimensional QKSD solve.

    ``beta_oracle``: a TOTAL-energy separator known to satisfy beta <= E_2 of the reachable sector
    (``reachable_separator`` supplies it on systems small enough to diagonalise). None selects the
    self-certified estimate, which is NOT a certificate -- see the module docstring. Pass a shared
    ``solver`` to reuse the cached Krylov basis.
    """
    solver = solver if solver is not None else QuantumKrylovSolver(mh)
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    _, states = solver.eigenstates(krylov_dim, n_states=3)
    beta = None if beta_oracle is None else float(beta_oracle) - mh.energy_offset
    return bracket_from_states(H, states, beta=beta, m=krylov_dim, offset=mh.energy_offset)


def bracket_ladder(mh: MolecularHamiltonian, dims: Sequence[int],
                   beta_oracle: Optional[float] = None,
                   solver: Optional[QuantumKrylovSolver] = None) -> List[ExcitedBracket]:
    """Brackets at each Krylov dimension in ``dims`` (Krylov basis built once and reused)."""
    solver = solver if solver is not None else QuantumKrylovSolver(mh)
    return [lehmann_excited_brackets(mh, m, beta_oracle=beta_oracle, solver=solver) for m in dims]


def reachable_separator(mh: MolecularHamiltonian) -> float:
    """REFERENCE ONLY (dense, O(2^n)): E_2 of the HF-reachable sector, as a TOTAL energy.

    The largest valid separator. +inf when the sector holds fewer than three levels -- then E_1 is
    the top of the sector and every beta > theta_1 is admissible; -inf is never returned.
    Validation oracle, never a live path: it diagonalises H exactly.
    """
    w, _ = reachable_eigenpairs(mh)
    return float(w[2]) + mh.energy_offset if len(w) > 2 else np.inf


if __name__ == "__main__":
    import csv
    import os

    from pyscf import ao2mo, fci, gto, mcscf, scf

    from hybrid_quantum_solver.molecular_hamiltonian import (build_hamiltonian_from_integrals,
                                                             build_molecular_hamiltonian)

    OUTPUT = "data/excited_certification_bench.csv"
    DIMS = (4, 6, 8, 10, 12, 16, 20, 24)
    CASES = {
        "LiH CAS(2,5)": dict(atom="Li 0 0 0; H 0 0 1.6", active_electrons=2, active_orbitals=5),
        "H4 chain": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
        "N2 CAS(6,6)": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
        "H2 STO-3G": dict(atom="H 0 0 0; H 0 0 0.74"),
    }
    rows = []

    def record(name, br, e0, e1, beta_is_oracle):
        gap = e1 - e0
        ritz_err = abs((br.upper_1 - br.upper_0) - gap)
        rows.append(dict(
            system=name, m=br.m, rank=br.rank, mode=br.beta_source,
            theta0=br.upper_0, theta1=br.upper_1, lower_0=br.lower_0, lower_1=br.lower_1,
            temple_lower_1=br.temple_lower_1, width_1=br.width_1,
            gap_lower=br.gap_lower, gap_upper=br.gap_upper, gap_exact=gap,
            ritz_gap_error=ritz_err,
            overhead=(br.gap_width / ritz_err) if ritz_err > 0 and np.isfinite(br.gap_width)
            else float("inf"),
            beta=br.beta, e0_ref=e0, e1_ref=e1,
            contains_e0=bool(br.lower_0 <= e0 <= br.upper_0),
            contains_e1=bool(br.lower_1 <= e1 <= br.upper_1),
        ))

    for name, spec in CASES.items():
        mh = build_molecular_hamiltonian(**spec)
        w, _ = reachable_eigenpairs(mh)
        e0, e1 = float(w[0]) + mh.energy_offset, float(w[1]) + mh.energy_offset
        beta = reachable_separator(mh)
        solver = QuantumKrylovSolver(mh)
        print("=" * 100)
        print(f"{name}: reachable sector holds {len(w)} levels; E_0={e0:.6f} E_1={e1:.6f} "
              f"gap={(e1 - e0) * 1e3:.3f} mHa; separator E_2={beta:.6f}")
        print("   M | U_1 - L_1 (mHa) | gap bracket (mHa)          | overhead | E_0,E_1 inside? "
              "| self-mode E_1 inside?")
        for m in DIMS:
            br = lehmann_excited_brackets(mh, m, beta_oracle=beta, solver=solver)
            sb = lehmann_excited_brackets(mh, m, beta_oracle=None, solver=solver)
            record(name, br, e0, e1, True)
            record(name, sb, e0, e1, False)
            ok = (br.lower_0 <= e0 <= br.upper_0) and (br.lower_1 <= e1 <= br.upper_1)
            gap = e1 - e0
            ritz_err = abs((br.upper_1 - br.upper_0) - gap)
            over = br.gap_width / ritz_err if ritz_err > 0 else float("nan")
            print(f"  {m:2d} | {br.width_1 * 1e3:15.4f} | "
                  f"[{br.gap_lower * 1e3:9.3f}, {br.gap_upper * 1e3:9.3f}] | {over:8.2f} | "
                  f"{str(ok):>5s} | {str(sb.lower_1 <= e1):>5s}")

    # Be2 -- the spec's showcase, kept as the refusal it is. Too large to diagonalise densely, so
    # the separator comes from a singlet-only CASCI (the Krylov space stays in the singlet sector
    # H conserves S^2 and |HF> is a singlet, so any singlet root is a conservative separator).
    print("=" * 100)
    R = 2.45
    mol = gto.M(atom=f"Be 0 0 0; Be 0 0 {R}", basis="ccpvdz", spin=0, verbose=0)
    mf = scf.RHF(mol).run()
    cas = mcscf.CASCI(mf, 8, 4)
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), 8)
    ne = (int(cas.nelecas[0]), int(cas.nelecas[1]))
    roots = fci.direct_spin0.kernel(h1, eri, 8, ne, ecore=float(e_core), nroots=3)[0]
    mh = build_hamiltonian_from_integrals(h1, eri, num_particles=ne, energy_offset=float(e_core))
    solver = QuantumKrylovSolver(mh)
    print(f"Be2 CAS(4,8)/cc-pVDZ at R={R} A: singlet CASCI roots {np.round(roots, 6)} -- "
          f"E_1 = E_2 (degenerate pi pair), so NO separator theta_1 < beta <= E_2 exists")
    print("   M | theta_1 - E_1^CASCI (mHa) | L_1        | rank")
    for m in (6, 8, 10, 12, 14):
        br = lehmann_excited_brackets(mh, m, beta_oracle=float(roots[2]), solver=solver)
        record("Be2 CAS(4,8)", br, float(roots[0]), float(roots[1]), True)
        print(f"  {m:2d} | {(br.upper_1 - roots[1]) * 1e3:25.3f} | {br.lower_1:10.4f} | {br.rank}")
    print("Be2 refuses at every depth: the bright pi states are not HF-reachable by real-time "
          "evolution, and a degenerate E_1 admits no separator in any case.")

    os.makedirs("data", exist_ok=True)
    with open(OUTPUT, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {OUTPUT} ({len(rows)} rows)")
