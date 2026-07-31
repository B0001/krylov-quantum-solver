#!/usr/bin/env python3
"""
Certified Hartree-Fock guiding-state overlap on molecules -- SPEC-21 meets the pipeline.

The guided-LH results (arXiv:2111.09079 and successors) make ground-state energy estimation
BQP-complete GIVEN a guiding state with overlap gamma >= 1/poly(n) -- an assumption that is
asserted, never certified, in that pipeline. This module produces the certificate: it feeds
the repo's own premise-gated Krylov E_1 floor (certified_gaps.gap_bracket, self mode eps_1 =
theta_1 - sigma_1, valid for M >= 6) into the SPEC-21 Davis-Kahan machinery and returns a
rigorous floor gamma_min <= |<HF|psi_0>| from Krylov data alone -- no oracle.

SECTOR HONESTY (specs/SPEC_hf_overlap_certificate.md): the certified object is the overlap
with the lowest HF-REACHABLE eigenstate, the same scope as QKSD/temple_bounds/certified_gaps.
The Davis-Kahan proof carries over with "spectrum" read as "reachable spectrum" because the
HF state's unreachable components vanish by construction (reachability is defined by the HF
state itself). Exact-statevector caveat inherited from temple_bounds: the shot cost of <H^2>
is not modeled.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from certified_gaps import gap_bracket
from hybrid_quantum_solver.certified_overlap import (
    GapCertificate,
    OverlapCertificate,
    certify_overlap,
)
from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

# The self-mode floor eps_1 = theta_1 - sigma_1 cannot verify its own premise below M = 6
# (the gated temple_bracket boundary). Inherited here as a hard raise, not re-derived.
_SELF_MODE_MIN_M = 6

# The HF-reachable sector is defined by |<HF|psi_k>|^2 > tol. This is the CERTIFIED ARC's value,
# shared with certified_gaps.py, certified_dipole.py and certified_noise.py -- but NOT with
# hf_overlap_subspace.py, which uses 1e-8. That divergence is not cosmetic: at square H4 a = 1.1 A
# the two thresholds select DIFFERENT ground states (HF overlap 2.25e-5 vs 0.667), so the d=1 and
# d=2 certificates whose head-to-head is SPEC_hf_overlap_subspace's headline are certifying
# different targets there. Named here so the divergence is visible and gated; unifying it decides
# which of two recorded findings is right and is deliberately left to a follow-up.
# See specs/SPEC_reachability_tolerance.md.
REACHABLE_TOL_CERTIFIED = 1e-10


def certify_hf_overlap(mh: MolecularHamiltonian, m: int = 8, e1: Optional[float] = None,
                       solver: Optional[QuantumKrylovSolver] = None) -> OverlapCertificate:
    """Certified lower bound on |<HF|psi_0>| for the reachable-sector ground state.

    ``e1``: exact E_1 as a TOTAL energy (oracle mode -- validation only), or None for the
    self-mode Krylov floor (production path; requires m >= 6). Pass a shared ``solver`` to
    reuse the cached Krylov basis. Everything runs in the electronic frame (overlaps are
    frame-independent; the E_1 floor and Rayleigh quotient must just agree, and do).

    Raises on: self mode below M = 6, a non-finite E_1 floor (rank-1 Krylov subspace), or
    any SPEC-21 invariant violation. Returns a possibly-VACUOUS certificate -- check
    ``.vacuous`` before quoting gamma_min.
    """
    if e1 is None and m < _SELF_MODE_MIN_M:
        raise ValueError(
            f"self-mode E1 floor requires m >= {_SELF_MODE_MIN_M} (the gated temple_bracket "
            f"premise boundary), got m = {m}. Supply an oracle e1 or increase m."
        )
    solver = solver if solver is not None else QuantumKrylovSolver(mh)
    gb = gap_bracket(mh, m, e1=e1, solver=solver)
    gap_cert = GapCertificate(
        e1_floor=gb.eps1,   # electronic frame; GapCertificate raises if -inf (rank-1)
        certificate_id=f"gap_bracket:{gb.eps1_source}:M={m}",
        source="oracle" if gb.eps1_source == "oracle" else "krylov_self_eps",
    )
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    u = np.asarray(mh.hf_state().data, dtype=complex)
    return certify_overlap(H, u, gap_cert, n_qubits=mh.qubit_hamiltonian.num_qubits)


def exact_reachable_overlap(mh: MolecularHamiltonian) -> float:
    """REFERENCE ONLY (dense, O(2^n)): exact |<HF|psi_0>| against the lowest HF-reachable
    eigenstate -- the killable check for the certificate. Never the live path."""
    H = mh.qubit_hamiltonian.to_matrix()
    w, V = np.linalg.eigh(H)
    hf = np.asarray(mh.hf_state().data, dtype=complex)
    amps = np.abs(V.conj().T @ hf)
    reachable = np.where(amps ** 2 > REACHABLE_TOL_CERTIFIED)[0]
    return float(amps[reachable[0]])


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    cases = {
        "H2 (0.74 A, equilibrium)": dict(atom="H 0 0 0; H 0 0 0.74"),
        "H2 (2.0 A, stretched)": dict(atom="H 0 0 0; H 0 0 2.0"),
        "H4 chain (1.0 A)": dict(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"),
    }
    print("=" * 88)
    print("Certified HF guiding-state overlap floors (self-mode Krylov E1 floor, no oracle)")
    print("   system                      |  M | gamma_min (self) | gamma_min (oracle) | exact")
    for name, spec in cases.items():
        mh = build_molecular_hamiltonian(**spec)
        solver = QuantumKrylovSolver(mh)
        exact = exact_reachable_overlap(mh)
        # oracle E1: exact reachable E1 as a total energy
        Hd = mh.qubit_hamiltonian.to_matrix()
        w, V = np.linalg.eigh(Hd)
        hf = np.asarray(mh.hf_state().data, dtype=complex)
        reach = np.where(np.abs(V.conj().T @ hf) ** 2 > REACHABLE_TOL_CERTIFIED)[0]
        e1_total = float(w[reach[1]]) + mh.energy_offset
        for m in (6, 8, 12):
            c_self = certify_hf_overlap(mh, m, solver=solver)
            c_orac = certify_hf_overlap(mh, m, e1=e1_total, solver=solver)
            gs = "VACUOUS" if c_self.vacuous else f"{c_self.gamma_min:.6f}"
            go = "VACUOUS" if c_orac.vacuous else f"{c_orac.gamma_min:.6f}"
            print(f"  {name:27s} | {m:2d} | {gs:>16s} | {go:>18s} | {exact:.6f}")
    print("=" * 88)
