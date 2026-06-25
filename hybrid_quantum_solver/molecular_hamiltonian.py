#!/usr/bin/env python3
"""
molecular_hamiltonian.py -- corrected electronic-structure -> qubit Hamiltonian.

This module REPLACES the hand-rolled ``AdvancedStochasticCompactor`` Jordan-Wigner
mapping in ``orchestrate_hybrid_pipeline.py``. That implementation (a) dropped electron
spin by mapping spatial orbitals directly to qubits, (b) used an incomplete two-body
expansion, and (c) did not reproduce Full Configuration Interaction (FCI) for even the
smallest molecule. See ``tests/test_reference_energies.py`` and ``REFACTOR_PLAN.md``.

Here we use Qiskit Nature's vetted ``PySCFDriver`` + ``JordanWignerMapper`` and expose
everything the quantum subspace solver needs:

  * ``qubit_hamiltonian`` -- a qiskit ``SparsePauliOp`` (electronic part; the constant
    nuclear/core energy is kept separately, NOT folded into the operator).
  * ``energy_offset``     -- add this to ANY eigenvalue of the qubit Hamiltonian to get
    a physical total energy (nuclear repulsion + frozen-core energy for active spaces).
  * ``hf_circuit``        -- the Hartree-Fock reference state |psi_0>, the correct
    starting point for a real-time Krylov subspace (the old code started from the
    particle-number-zero vacuum |00..0>, where <0|H|0> = 0).

Validation (full configuration space, cross-checked against an independent PySCF FCI):

    build_molecular_hamiltonian("H 0 0 0; H 0 0 0.74").ground_state_energy()  == -1.137284 Ha
    H4 chain                                                                  == -2.166387 Ha
    LiH (12 qubits)                                                           == -7.882324 Ha

and for an active space, it reproduces PySCF CASCI to < 1e-6 Ha.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit_nature.second_q.circuit.library import HartreeFock
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.hamiltonians import ElectronicEnergy
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer


@dataclass
class MolecularHamiltonian:
    """Container for a qubit Hamiltonian plus the bookkeeping a solver needs.

    The physical total energy of any eigenstate is ``eigenvalue + energy_offset``.
    """

    qubit_hamiltonian: SparsePauliOp     # electronic operator (offset excluded)
    energy_offset: float                 # nuclear repulsion (+ frozen-core) constant
    hf_circuit: QuantumCircuit           # Hartree-Fock reference |psi_0>
    num_qubits: int
    num_particles: Tuple[int, int]       # (n_alpha, n_beta)
    num_spatial_orbitals: int
    hf_energy: float                     # <HF|H|HF> + offset  (== RHF total energy)

    def total_energy(self, electronic_eigenvalue: float) -> float:
        """Lift an eigenvalue of ``qubit_hamiltonian`` into the physical energy frame."""
        return float(electronic_eigenvalue) + self.energy_offset

    def hf_state(self) -> Statevector:
        """The Hartree-Fock reference as a statevector."""
        return Statevector(self.hf_circuit)

    def ground_state_energy(self) -> float:
        """Exact lowest total energy via dense diagonalisation.

        Reference/validation helper only -- O(4^n) memory. The quantum subspace
        solver (Phase 2) approximates this without forming the full matrix.
        """
        eigenvalue = float(np.linalg.eigvalsh(self.qubit_hamiltonian.to_matrix())[0])
        return self.total_energy(eigenvalue)


def build_molecular_hamiltonian(
    atom: str,
    basis: str = "sto3g",
    charge: int = 0,
    spin: int = 0,
    active_electrons: Optional[int] = None,
    active_orbitals: Optional[int] = None,
    mapper: Optional[JordanWignerMapper] = None,
) -> MolecularHamiltonian:
    """Build the correct qubit Hamiltonian for a molecule.

    Args:
        atom: PySCF-style geometry string, e.g. ``"H 0 0 0; H 0 0 0.74"`` (Angstrom).
        basis: Atomic-orbital basis set.
        charge, spin: Molecular charge and ``2S`` spin (number of unpaired electrons).
        active_electrons, active_orbitals: If both given, restrict to a CAS active space
            (the frozen-core energy is folded into ``energy_offset``). If omitted, the
            full orbital space is used.
        mapper: Fermion-to-qubit mapper; defaults to Jordan-Wigner.

    Returns:
        A populated :class:`MolecularHamiltonian`.
    """
    mapper = mapper or JordanWignerMapper()

    problem = PySCFDriver(atom=atom, basis=basis, charge=charge, spin=spin).run()

    if active_electrons is not None and active_orbitals is not None:
        problem = ActiveSpaceTransformer(active_electrons, active_orbitals).transform(problem)
    elif (active_electrons is None) ^ (active_orbitals is None):
        raise ValueError("Provide BOTH active_electrons and active_orbitals, or neither.")

    qubit_hamiltonian = mapper.map(problem.hamiltonian.second_q_op())

    # Every constant the driver/transformer set aside (nuclear repulsion, and for an
    # active space the inactive/frozen-core energy) must be re-added to recover totals.
    energy_offset = float(sum(problem.hamiltonian.constants.values()))

    hf_circuit = HartreeFock(problem.num_spatial_orbitals, problem.num_particles, mapper)
    hf_energy = float(np.real(Statevector(hf_circuit).expectation_value(qubit_hamiltonian)))
    hf_energy += energy_offset

    return MolecularHamiltonian(
        qubit_hamiltonian=qubit_hamiltonian,
        energy_offset=energy_offset,
        hf_circuit=hf_circuit,
        num_qubits=qubit_hamiltonian.num_qubits,
        num_particles=problem.num_particles,
        num_spatial_orbitals=problem.num_spatial_orbitals,
        hf_energy=hf_energy,
    )


def build_hamiltonian_from_integrals(
    h1: np.ndarray,
    eri: np.ndarray,
    num_particles: Tuple[int, int],
    energy_offset: float = 0.0,
    mapper: Optional[JordanWignerMapper] = None,
) -> MolecularHamiltonian:
    """Build a qubit Hamiltonian from precomputed active-space integrals.

    This is the bridge for the CIF / CASCI materials path (see ``chemistry_gateway``):
    feed the PySCF active-space one-body integrals and two-body ERIs plus the core energy
    from ``cas.get_h1eff()``. Reproduces PySCF CASCI to <1e-6 Ha (verified for H2 and LiH).

    Args:
        h1: active-space one-body integrals in the MO basis, shape ``(n_orb, n_orb)``.
        eri: active-space two-body integrals in chemist notation ``(pq|rs)`` restored to
            four indices, shape ``(n_orb,) * 4`` (e.g. ``ao2mo.restore(1, cas.get_h2eff(), n)``).
        num_particles: ``(n_alpha, n_beta)`` electrons in the active space.
        energy_offset: constant added to recover total energies -- typically ``e_core`` from
            ``cas.get_h1eff()``, which already bundles nuclear repulsion + frozen-core energy.
        mapper: fermion-to-qubit mapper; defaults to Jordan-Wigner.

    Note:
        Assumes spin-restricted active-space integrals (the standard closed-shell CASCI case).
    """
    mapper = mapper or JordanWignerMapper()
    h1 = np.asarray(h1)
    eri = np.asarray(eri)
    n_orb = h1.shape[0]
    if eri.shape != (n_orb, n_orb, n_orb, n_orb):
        raise ValueError(f"eri must have shape {(n_orb,) * 4}, got {eri.shape}")

    qubit_hamiltonian = mapper.map(ElectronicEnergy.from_raw_integrals(h1, eri).second_q_op())
    offset = float(energy_offset)

    hf_circuit = HartreeFock(n_orb, tuple(num_particles), mapper)
    hf_energy = float(np.real(Statevector(hf_circuit).expectation_value(qubit_hamiltonian)))
    hf_energy += offset

    return MolecularHamiltonian(
        qubit_hamiltonian=qubit_hamiltonian,
        energy_offset=offset,
        hf_circuit=hf_circuit,
        num_qubits=qubit_hamiltonian.num_qubits,
        num_particles=tuple(num_particles),
        num_spatial_orbitals=n_orb,
        hf_energy=hf_energy,
    )


if __name__ == "__main__":
    for name, spec in {
        "H2": dict(atom="H 0 0 0; H 0 0 0.74"),
        "H4 chain": dict(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"),
        "LiH": dict(atom="Li 0 0 0; H 0 0 1.6"),
    }.items():
        mh = build_molecular_hamiltonian(**spec)
        print(f"{name:9s}  qubits={mh.num_qubits:2d}  terms={len(mh.qubit_hamiltonian):4d}  "
              f"E_HF={mh.hf_energy:.6f}  E_ground(FCI)={mh.ground_state_energy():.6f} Ha")
