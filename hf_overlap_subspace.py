#!/usr/bin/env python3
"""
Certified Hartree-Fock SUBSPACE overlap on molecules -- SPEC-21b meets the pipeline.

The d=1 sibling (hf_overlap_certificate.py) certifies |<HF|psi_0>| for the reachable ground
STATE. But on a strongly-multireference molecule HF is a poor proxy for the ground state and a
strong proxy for the low-lying ground EIGENSPACE: its weight spreads across the two lowest
reachable levels, its overlap with either one alone drops below the residual-to-gap ratio, and
the single-vector certificate goes vacuous -- while its overlap with the two-level eigenspace
stays large.

This feeds the repo's own premise-gated self-mode Krylov E_d floor (theta_d - sigma_d, the
Weinstein floor generalized from the d=1 gap_bracket) into the SPEC-21b block Davis-Kahan
machinery and certifies gamma_min <= ||P_S u|| for the lowest-d reachable eigenspace, with NO
oracle. Worked system: square H4 (see specs/SPEC_hf_overlap_subspace.md).

SECTOR HONESTY: reachable-sector restricted exactly like QKSD/temple_bounds/certified_gaps/the
d=1 path -- P_S spans the lowest d REACHABLE eigenstates, beta floors the (d+1)-th reachable
level; unreachable levels have zero HF amplitude and are irrelevant to both r^2 and sin^2(theta).
Exact-statevector caveat inherited from temple_bounds (the <H^2> shot cost is not modeled).
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from hybrid_quantum_solver.certified_overlap import (
    ClusterGapCertificate,
    OverlapCertificate,
    certify_subspace_overlap,
    rayleigh_quotient,
    residual_norm,
)
from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

# The self-mode Weinstein floor cannot verify its own premise below M = 6 (the gated
# temple_bracket boundary). Inherited here as a hard raise, not re-derived.
_SELF_MODE_MIN_M = 6
_REACHABLE_TOL = 1e-8


def _weinstein_intervals_disjoint(centers, sigmas) -> bool:
    """Are the Weinstein intervals [theta_k - sigma_k, theta_k + sigma_k] pairwise disjoint?

    The in-band resolvability signal (SPEC_subspace_floor_resolvability.md): when the Krylov space
    has NOT resolved the size-d cluster, the residuals sigma_k are large and adjacent intervals
    overlap, so the self-mode floor theta_d - sigma_d is untrustworthy (it can exceed the true
    E_d). ``centers`` must be ascending. A necessary, oracle-free check -- NOT proven sufficient.
    """
    for k in range(len(centers) - 1):
        if centers[k] + sigmas[k] >= centers[k + 1] - sigmas[k + 1]:
            return False
    return True


def certify_hf_subspace_overlap(
    mh: MolecularHamiltonian, cluster_size: int, m: int = 8,
    e_d: Optional[float] = None, solver: Optional[QuantumKrylovSolver] = None,
) -> OverlapCertificate:
    """Certified lower bound on ||P_S u|| for the lowest-d reachable eigenspace, u = HF.

    ``cluster_size`` d: dimension of the reachable ground eigenspace to certify.
    ``e_d``: exact (d+1)-th reachable level as a TOTAL energy (oracle mode -- validation), or None
    for the self-mode Krylov floor theta_d - sigma_d (production path; requires m >= 6). Pass a
    shared ``solver`` to reuse the Krylov basis.

    Raises on: self mode below M = 6, a Krylov subspace too small to expose level d (rank
    < d+1), or any SPEC-21b invariant. Returns a possibly-VACUOUS certificate -- check
    ``.vacuous`` before quoting gamma_min. In self mode the certificate is also VACUOUS when the
    Krylov space has not resolved the cluster (overlapping Weinstein intervals), so the self-mode
    floor is never emitted unsound -- see SPEC_subspace_floor_resolvability.md.
    """
    if e_d is None and m < _SELF_MODE_MIN_M:
        raise ValueError(
            f"self-mode E_d floor requires m >= {_SELF_MODE_MIN_M} (the gated temple_bracket "
            f"premise boundary), got m = {m}. Supply an oracle e_d or increase m."
        )
    solver = solver if solver is not None else QuantumKrylovSolver(mh)
    offset = mh.energy_offset
    Hs = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    u = np.asarray(mh.hf_state().data, dtype=complex)

    energies, states = solver.eigenstates(m, n_states=cluster_size + 1)
    if len(energies) < cluster_size + 1:
        raise ValueError(
            f"Krylov subspace at M = {m} exposes only {len(energies)} reachable levels; "
            f"need d + 1 = {cluster_size + 1} to floor E_d. Increase m."
        )

    if e_d is None:
        # Resolvability guard: the self-mode floor theta_d - sigma_d is only trustworthy when the
        # Krylov space has resolved the cluster. Unresolved => overlapping Weinstein intervals =>
        # return VACUOUS rather than a possibly-unsound positive floor (fail-safe, not fail-silent).
        centers = [energies[k] - offset for k in range(cluster_size + 1)]
        sigmas = [residual_norm(Hs, states[k], centers[k]) for k in range(cluster_size + 1)]
        if not _weinstein_intervals_disjoint(centers, sigmas):
            lam = rayleigh_quotient(Hs, u)
            return OverlapCertificate(
                gamma_min=0.0, lambda_u=lam, residual_norm=residual_norm(Hs, u, lam),
                gap_certificate_id=f"hf_subspace:self:M={m}:d={cluster_size}:UNRESOLVED",
                vacuous=True, cluster_size=cluster_size,
                vacuous_reason=(
                    f"self-mode floor unresolved at M={m}: the {cluster_size + 1} lowest Weinstein "
                    "intervals overlap, so theta_d - sigma_d is not a trustworthy E_d floor. "
                    "Increase M or supply an oracle e_d."
                ),
            )
        e_above_floor = centers[cluster_size] - sigmas[cluster_size]   # Weinstein floor, electronic
        source = "krylov_self_eps"
        cert_id = f"hf_subspace:self:M={m}:d={cluster_size}"
    else:
        e_above_floor = float(e_d) - offset              # oracle total -> electronic
        source = "oracle"
        cert_id = f"hf_subspace:oracle:d={cluster_size}"

    gap_cert = ClusterGapCertificate(
        e_above_floor=e_above_floor, cluster_size=cluster_size,
        certificate_id=cert_id, source=source,
    )
    return certify_subspace_overlap(Hs, u, gap_cert, n_qubits=mh.qubit_hamiltonian.num_qubits)


def exact_hf_subspace_overlap(mh: MolecularHamiltonian, cluster_size: int,
                              tol: float = _REACHABLE_TOL) -> float:
    """REFERENCE ONLY (dense, O(2^n)): exact ||P_S u|| for the lowest-d REACHABLE eigenspace --
    the killable check. Reachable = nonzero HF amplitude (the QKSD sector). Never the live path."""
    H = mh.qubit_hamiltonian.to_matrix()
    w, V = np.linalg.eigh(H)
    u = np.asarray(mh.hf_state().data, dtype=complex)
    reach = np.where(np.abs(V.conj().T @ u) ** 2 > tol)[0]   # reachable indices, ascending energy
    P_S = V[:, reach[:cluster_size]]
    return float(np.linalg.norm(P_S.conj().T @ u))


def _reachable_e_d_total(mh: MolecularHamiltonian, cluster_size: int,
                         tol: float = _REACHABLE_TOL) -> float:
    """Oracle (d+1)-th reachable level as a TOTAL energy (dense; validation only)."""
    H = mh.qubit_hamiltonian.to_matrix()
    w, V = np.linalg.eigh(H)
    u = np.asarray(mh.hf_state().data, dtype=complex)
    reach = np.where(np.abs(V.conj().T @ u) ** 2 > tol)[0]
    return float(w[reach[cluster_size]]) + mh.energy_offset


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
    from hf_overlap_certificate import certify_hf_overlap

    print("=" * 92)
    print("Certified HF subspace overlap on square H4 (d=2 reachable cluster, self-mode, no oracle)")
    print("  side a |  d=1 cert (single-vector) |  d=2 gamma_min (self) | d=2 (oracle) | exact ||P_S u||")
    for a in (1.4, 1.3, 1.2, 1.1, 1.0):
        atom = f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0"
        mh = build_molecular_hamiltonian(atom=atom)
        solver = QuantumKrylovSolver(mh)
        c1 = certify_hf_overlap(mh, m=8, solver=solver)
        c2s = certify_hf_subspace_overlap(mh, 2, m=8, solver=solver)
        c2o = certify_hf_subspace_overlap(mh, 2, m=8, e_d=_reachable_e_d_total(mh, 2), solver=solver)
        exact = exact_hf_subspace_overlap(mh, 2)
        d1 = "VACUOUS" if c1.vacuous else f"{c1.gamma_min:.4f}"
        print(f"   {a:.1f}   | {d1:>24s} | {c2s.gamma_min:>21.4f} | {c2o.gamma_min:>12.4f} | {exact:.4f}")
    print("=" * 92)
