"""
Sprint 10: Quantum Krylov Oracle
Executes qDRIFT compiled sequences into Qiskit circuits and measures 
subspace expectation values for the generalized eigenvalue problem.
"""

import numpy as np
from typing import List, Dict, Any
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector, SparsePauliOp
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.synthesis import LieTrotter


class QiskitKrylovSampler:
    def __init__(self, n_qubits: int):
        self.n_qubits = n_qubits
        print(f"[QUANTUM NODE] Initializing Qiskit Oracle for {n_qubits} qubits.")

    def _build_hamiltonian_operator(self, pauli_strings: List[str], coefficients: List[float]) -> SparsePauliOp:
        """Reconstructs the full target Hamiltonian for expectation measurements."""
        # Qiskit reads Pauli strings right-to-left (endianness).
        # We reverse the strings to match the standard left-to-right JW notation.
        reversed_paulis = [p[::-1] for p in pauli_strings]
        return SparsePauliOp(reversed_paulis, coefficients)

    def execute_subspace_sampling(
        self, 
        compiled_circuits: List[Dict[str, Any]], 
        subspace_dim: int,
        full_pauli_strings: List[str],
        full_coefficients: List[float]
    ) -> List[Dict[str, Any]]:
        """
        Generates Krylov basis states from the qDRIFT circuit sequence 
        and calculates the H and S matrix elements.
        """
        print(f"[QUANTUM NODE] Synthesizing {subspace_dim}-dimensional Krylov basis...")
        
        # 1. Reconstruct the observable Hamiltonian
        H_obs = self._build_hamiltonian_operator(full_pauli_strings, full_coefficients)
        
        # 2. Build the reference state |00...0>
        qc_ref = QuantumCircuit(self.n_qubits)
        basis_states = [Statevector(qc_ref)]
        
        # 3. Apply the qDRIFT evolution steps to generate the subspace states
        current_qc = qc_ref.copy()
        
        # Initialize the Trotter synthesizer (reps=1 because we only have 1 Pauli string per step)
        synthesizer = LieTrotter(reps=1)

        for step in range(subspace_dim - 1):
            if step < len(compiled_circuits):
                circuit_data = compiled_circuits[step]
                pauli_op = SparsePauliOp(circuit_data["target_pauli"][::-1])
                tau = circuit_data["evolution_angle_tau"]
                
                # Create the abstract gate
                evo_gate = PauliEvolutionGate(pauli_op, time=tau)
                
                # SYNTHESIS: Decompose the abstract exponential into native quantum gates (CX, Rz)
                synthesized_step = synthesizer.synthesize(evo_gate)
                
                # Append the decomposed circuit instructions directly
                current_qc.compose(synthesized_step, inplace=True)
                
                basis_states.append(Statevector(current_qc))
            else:
                basis_states.append(basis_states[-1])
                
        # 4. Measure the Subspace Matrices (H_ij and S_ij)
        print(f"[QUANTUM NODE] Measuring expectation values across O(M^2) tensor grid...")
        mapped_samples = []
        for i in range(subspace_dim):
            for j in range(subspace_dim):
                # Overlap: S_ij = <psi_i | psi_j>
                s_ij = basis_states[i].inner(basis_states[j])
                
                # Energy: H_ij = <psi_i | H | psi_j>
                h_ket = basis_states[j].evolve(H_obs)
                h_ij = basis_states[i].inner(h_ket)
                
                # SVD Shifter requires real-valued Hermitian components
                mapped_samples.append({
                    "row": i,
                    "col": j,
                    "h_val": float(np.real(h_ij)),
                    "s_val": float(np.real(s_ij))
                })
                
        return mapped_samples