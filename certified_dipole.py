#!/usr/bin/env python3
"""
Certified error bars on a molecular PROPERTY (dipole moment) from Krylov data -- no FCI oracle.

The certified arc so far bounds energies (`temple_bounds`), gaps (`certified_gaps`), and
self-verifies them (`gap_selfcheck`). But the observables experiments report -- dipole moments,
and any <psi_0|A|psi_0> -- were still bare point estimates (`qksd_properties`). This closes that:
a rigorous interval around the ground-state dipole, from the SAME Krylov data.

The mechanism composes the arc. If the ground Ritz state psi_0 (Rayleigh quotient rho_0, residual
sigma_0 = ||(H-rho_0)psi_0||) has angle theta to the exact eigenstate, the eigenvector-perturbation
(Davis-Kahan / gap) theorem gives sin theta <= sigma_0 / (E_1 - rho_0), and E_1 - rho_0 >= Delta_lo,
the certified GAP lower bound from `certified_gaps`. So

    sin theta <= s := sigma_0 / Delta_lo                                  (inherits the gap premise)

and, for a bounded Hermitian A, the SHARP first-order bound (fluctuation, not operator norm):

    | <psi_0|A|psi_0> - <exact|A|exact> |  <=  2 sigma_A(psi_0) sin theta + W_A sin^2 theta
                                            <=  2 sigma_A s + W_A s^2  =: half_width

with sigma_A(psi_0) = sqrt(<A^2> - <A>^2) the dipole fluctuation in the trial state and W_A =
lambda_max - lambda_min the operator's spectral width. Using sigma_A instead of ||A|| is the whole
game: for LiH sigma_A ~ 1.1 vs ||mu_z|| ~ 6.9, a ~6x tighter (and non-vacuous) interval.

THE FINDING (specs/SPEC_certified_dipole.md): the exact FCI dipole lies inside
[mu +/- half_width] at every depth (zero escapes), and the interval closes -- LiH reaches
-1.818 +/- 0.065 a.u. at M=24 (exact -1.817). **The property certificate INHERITS the gap
certificate:** half_width is finite iff s < 1, i.e. iff sigma_0 < Delta_lo, so it is vacuous exactly
at the depths where the certified gap lower bound is weak (LiH M=8-16) and sharp where Delta_lo is
healthy (M >= 20). Certified properties are only as good as the certified gap beneath them; pair
with `gap_selfcheck` to know when Delta_lo is trustworthy.

HONEST SCOPE: sector-restricted ground state (QKSD's scope); the sin theta bound rests on the same
premise as `certified_gaps` (valid M >= 6, oracle-free-checkable via `gap_selfcheck`); exact
statevector (the shot cost of sigma_0, sigma_A is not modeled); full orbital space (the dipole
convention of `build_dipole_operators`).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np
from scipy.sparse.linalg import eigsh

from certified_gaps import gap_bracket
from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from reachability import reachable_eigenpairs
from temple_bounds import mean_and_variance


@dataclass
class CertifiedDipole:
    """A certified interval [mu - half_width, mu + half_width] on a ground-state property."""
    m: int                  # Krylov dimension
    mu: float               # <psi_0|A|psi_0> point estimate
    half_width: float       # certified half-width (inf when s >= 1: valid but vacuous)
    sin_theta_bound: float  # s = sigma_0 / Delta_lo (eigenvector rotation bound)
    sigma_A: float          # dipole fluctuation in psi_0
    gap_lower: float        # certified gap lower bound Delta_lo (Ha) feeding s
    finite: bool            # half_width < inf (s < 1)


def spectral_width(a_sparse) -> float:
    """W_A = lambda_max - lambda_min of a Hermitian operator (two Lanczos extremal solves)."""
    hi = float(eigsh(a_sparse, k=1, which="LA", return_eigenvectors=False)[0])
    lo = float(eigsh(a_sparse, k=1, which="SA", return_eigenvectors=False)[0])
    return hi - lo


def certified_dipole(mh: MolecularHamiltonian, a_sparse, m: int, width: Optional[float] = None,
                     solver: Optional[QuantumKrylovSolver] = None,
                     e1: Optional[float] = None) -> CertifiedDipole:
    """Certified dipole interval from an M-dim Krylov solve. ``a_sparse`` is the dipole operator as
    a sparse matrix (same qubit basis as ``mh``); ``width`` its spectral width (computed if None).
    ``e1``: exact E_1 (oracle) for Delta_lo, else the self estimate (see certified_gaps)."""
    solver = solver if solver is not None else QuantumKrylovSolver(mh)
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    _, states = solver.eigenstates(m, n_states=1)
    psi0 = states[0]
    _, var0 = mean_and_variance(H, psi0)
    sigma0 = float(np.sqrt(var0))
    ax = a_sparse @ psi0
    mu = float((psi0.conj() @ ax).real)
    sigma_A = float(np.sqrt(max((ax.conj() @ ax).real - mu * mu, 0.0)))
    W = spectral_width(a_sparse) if width is None else width
    dlo = gap_bracket(mh, m, e1=e1, solver=solver).gap_lower
    s = sigma0 / dlo if dlo > 0 else np.inf
    half = 2.0 * sigma_A * s + W * s * s if s < 1.0 else np.inf
    return CertifiedDipole(m=m, mu=mu, half_width=half, sin_theta_bound=s, sigma_A=sigma_A,
                           gap_lower=dlo, finite=np.isfinite(half))


def certified_dipole_ladder(mh: MolecularHamiltonian, a_sparse, dims: Sequence[int],
                            solver: Optional[QuantumKrylovSolver] = None) -> List[CertifiedDipole]:
    """Certified dipole at each Krylov dimension (spectral width computed once)."""
    solver = solver if solver is not None else QuantumKrylovSolver(mh)
    width = spectral_width(a_sparse)
    return [certified_dipole(mh, a_sparse, m, width=width, solver=solver) for m in dims]


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import (
        build_dipole_operators,
        build_molecular_hamiltonian,
    )

    cases = {
        "HeH+": dict(atom="He 0 0 0; H 0 0 0.772", charge=1),
        "LiH": dict(atom="Li 0 0 0; H 0 0 1.6"),
    }
    for name, spec in cases.items():
        mh = build_molecular_hamiltonian(**spec)
        Az = build_dipole_operators(**spec)[2].to_matrix(sparse=True)
        psi_ex = reachable_eigenpairs(mh)[1][:, 0]   # HF-reachable ground (right charge sector)
        mu_exact = float((psi_ex.conj() @ (Az @ psi_ex)).real)
        solver = QuantumKrylovSolver(mh)
        print("=" * 72)
        print(f"{name}: exact mu_z = {mu_exact:+.4f} a.u. (certified WITHOUT it)")
        print("   M |    mu_z   | certified half-width (a.u.) | inside?")
        for cd in certified_dipole_ladder(mh, Az, (6, 8, 12, 16, 20, 24), solver=solver):
            hw = f"{cd.half_width:.4f}" if cd.finite else "inf (Delta_lo weak)"
            inside = cd.mu - cd.half_width <= mu_exact <= cd.mu + cd.half_width
            print(f"  {cd.m:2d} | {cd.mu:+.4f} | {hw:>26} | {inside}")
    print("=" * 72)
    print("Certified dipole = point estimate +/- Davis-Kahan(residual / certified gap lower bound).")
