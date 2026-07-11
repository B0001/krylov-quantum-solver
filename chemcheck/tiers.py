"""Frozen ChemCheck tier registry (task 1).

T0–T3 ship as immutable data — geometry, active space, exact FCI reference, accuracy
thresholds, reference two-qubit-gate count, and a canonical Hamiltonian SHA-256 — so scores are
comparable within a benchmark version. Loading :data:`TIERS` pulls in **no** solver dependency;
:func:`build_tier_hamiltonian` / :func:`canonical_hamiltonian_sha256` import PySCF/qiskit lazily
and exist so a test can re-verify the frozen numbers on any machine.

Frozen values were produced with ``basis="sto3g"`` and the same build path used here (full
active space for T0/T1, explicit CAS for T2/T3), matching ``benchmark_resources.py``.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from typing import Any, Optional

#: Version string; scores are only comparable within a version (schema pattern chemcheck-YYYY.N).
BENCHMARK_VERSION = "chemcheck-2026.1"

#: Chemical accuracy and the marginal band, in milli-Hartree.
ACCURACY_PASS_MHA = 1.6
ACCURACY_MARGINAL_MHA = 16.0


@dataclass(frozen=True)
class Tier:
    """One benchmark tier. ``aspirational`` tiers (T4) carry no frozen reference yet."""

    name: str
    system: str
    atom: str
    basis: str
    active_electrons: Optional[int]
    active_orbitals: Optional[int]
    spin_orbitals: Optional[int]
    fci_reference_hartree: Optional[float]
    two_qubit_gates_per_trotter_step: Optional[int]
    hamiltonian_pauli_terms: Optional[int]
    hamiltonian_sha256: Optional[str]
    classically_simulable: bool
    aspirational: bool = False


TIERS: dict[str, Tier] = {
    "T0": Tier(
        name="T0", system="H2 / STO-3G (sanity floor)",
        atom="H 0 0 0; H 0 0 0.74", basis="sto3g",
        active_electrons=None, active_orbitals=None, spin_orbitals=4,
        fci_reference_hartree=-1.137283834,
        two_qubit_gates_per_trotter_step=70, hamiltonian_pauli_terms=15,
        hamiltonian_sha256="26e50a294c09a23e602a52683e0006cf472d9aee6baf1a40e45ca2abf660558e",
        classically_simulable=True,
    ),
    "T1": Tier(
        name="T1", system="H4 chain / STO-3G (first multireference stress)",
        atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7", basis="sto3g",
        active_electrons=None, active_orbitals=None, spin_orbitals=8,
        fci_reference_hartree=-2.180316614,
        two_qubit_gates_per_trotter_step=2174, hamiltonian_pauli_terms=185,
        hamiltonian_sha256="2fe671eae4ae3b1ba43e3f669ac7bd87237c51d36e961d82cbe870ad85f2a1e8",
        classically_simulable=True,
    ),
    "T2": Tier(
        name="T2", system="LiH / STO-3G CAS(2,5) (ionic + covalent)",
        atom="Li 0 0 0; H 0 0 1.6", basis="sto3g",
        active_electrons=2, active_orbitals=5, spin_orbitals=10,
        fci_reference_hartree=-7.882096600,
        two_qubit_gates_per_trotter_step=3858, hamiltonian_pauli_terms=276,
        hamiltonian_sha256="52794491d99533ae5848da27e56ff307198375d191d271c5c82939409d9ff663",
        classically_simulable=True,
    ),
    "T3": Tier(
        name="T3", system="N2 stretched / STO-3G CAS(6,6) (strongly multireference)",
        atom="N 0 0 0; N 0 0 1.1", basis="sto3g",
        active_electrons=6, active_orbitals=6, spin_orbitals=12,
        fci_reference_hartree=-107.623101772,
        two_qubit_gates_per_trotter_step=6534, hamiltonian_pauli_terms=383,
        hamiltonian_sha256="382de579ca32d5e877b5292a169642bed38117085b342bb05772f4166c928460",
        classically_simulable=False,  # ~6,500 CX/step — the current hardware wall
    ),
    "T4": Tier(
        name="T4", system="Nb3X8-class cluster active space (materials target)",
        atom="", basis="sto3g",
        active_electrons=None, active_orbitals=None, spin_orbitals=None,
        fci_reference_hartree=None, two_qubit_gates_per_trotter_step=None,
        hamiltonian_pauli_terms=None, hamiltonian_sha256=None,
        classically_simulable=False, aspirational=True,
    ),
}


def build_tier_hamiltonian(tier: Tier):
    """Build the tier's MolecularHamiltonian (lazy solver import). Raises on aspirational tiers."""
    if tier.aspirational:
        raise ValueError(f"tier {tier.name} is aspirational — no Hamiltonian defined yet")
    from hybrid_quantum_solver import build_molecular_hamiltonian

    kwargs: dict[str, Any] = dict(atom=tier.atom, basis=tier.basis)
    if tier.active_electrons is not None and tier.active_orbitals is not None:
        kwargs["active_electrons"] = tier.active_electrons
        kwargs["active_orbitals"] = tier.active_orbitals
    return build_molecular_hamiltonian(**kwargs)


def canonical_hamiltonian_sha256(mh) -> str:
    """Deterministic content hash of a qubit Hamiltonian.

    Sorted ``(pauli_label, Re(coeff), Im(coeff))`` stream with coefficients quantized to 1e-12,
    JSON-serialized, SHA-256'd — stable across machines and Pauli orderings.
    """
    op = mh.qubit_hamiltonian
    terms = sorted(
        (p.to_label(), round(float(c.real), 12), round(float(c.imag), 12))
        for p, c in zip(op.paulis, op.coeffs)
    )
    return hashlib.sha256(json.dumps(terms, separators=(",", ":")).encode()).hexdigest()


def recompute_tier_reference(tier: Tier) -> dict[str, Any]:
    """Recompute the live FCI / hash / CX-count for cross-machine verification of the freeze."""
    from hybrid_quantum_solver.hardware_krylov import HardwareKrylovSolver

    mh = build_tier_hamiltonian(tier)
    report = HardwareKrylovSolver(mh).resource_report(4, shots=8192)
    return {
        "fci_reference_hartree": mh.ground_state_energy(),
        "hamiltonian_sha256": canonical_hamiltonian_sha256(mh),
        "two_qubit_gates_per_trotter_step": report["trotter_step_cx"],
        "hamiltonian_pauli_terms": report["hamiltonian_pauli_terms"],
    }
