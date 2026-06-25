#!/usr/bin/env python3
"""
hardware_krylov.py -- on-hardware quantum Krylov: subspace matrices from circuit measurements.

``quantum_krylov_solver.py`` (exact statevector) and ``trotter_krylov.py`` both form the subspace
matrices Sᵢⱼ, Hᵢⱼ from statevector inner products. THIS module measures them the way a real device
would: an ancilla **Hadamard test** with controlled Trotter evolution yields the real and imaginary
parts of ⟨φᵢ|φⱼ⟩ and ⟨φᵢ|H|φⱼ⟩ directly, evaluated through qiskit-aer's Estimator (exact,
finite-shot, or under a device ``NoiseModel``).

For each pair i ≤ j (φₖ = Uᵏ|φ_HF⟩, U = Trotter step):

    prepare |φ_HF⟩ ─ apply Uⁱ ─ H(ancilla) ─ controlled-U^(j−i) ─ measure
        Sᵢⱼ = ⟨Xₐ⊗I⟩ + i⟨Yₐ⊗I⟩ ,   Hᵢⱼ = ⟨Xₐ⊗H⟩ + i⟨Yₐ⊗H⟩

Measuring each pair directly (rather than assuming the Hermitian-Toeplitz structure of exact
evolution) keeps the result exact for the *Trotterised* basis, so the generalized eigenproblem is
a true Rayleigh quotient over genuine states and respects the variational floor (up to shot/device
noise). This is the hardware-faithful counterpart of the statevector solver, and the natural place
to add error mitigation (ZNE / readout) next.
"""
from __future__ import annotations

from typing import List, Optional, Sequence

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import SparsePauliOp
from scipy.sparse.linalg import eigsh

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.noise import extrapolate_zero_noise, fold_global_circuit
from hybrid_quantum_solver.quantum_krylov_solver import KrylovStep, solve_generalized_eig
from hybrid_quantum_solver.trotter_krylov import build_trotter_step

_BASIS_GATES = ["cx", "u", "rz", "h", "x", "sx", "p", "s", "sdg"]


class HardwareKrylovSolver:
    """Quantum Krylov whose subspace matrices are measured by Hadamard tests via qiskit-aer."""

    def __init__(
        self,
        molecular_hamiltonian: MolecularHamiltonian,
        dt: Optional[float] = None,
        trotter_order: int = 2,
        trotter_reps: int = 1,
        threshold: float = 1e-10,
        shots: Optional[int] = None,
        noise_model=None,
        seed: Optional[int] = None,
        zne_scale_factors: Optional[Sequence[int]] = None,
        zne_order: int = 1,
    ):
        self.mh = molecular_hamiltonian
        self.offset = molecular_hamiltonian.energy_offset
        self.threshold = threshold
        self.shots = shots
        self.zne_scale_factors = list(zne_scale_factors) if zne_scale_factors else None
        self.zne_order = zne_order
        self._N = molecular_hamiltonian.num_qubits
        self._Hop = molecular_hamiltonian.qubit_hamiltonian
        self._Hsparse = self._Hop.to_matrix(sparse=True).tocsc()
        self.dt = float(dt) if dt is not None else self._default_dt()
        self.trotter_order = trotter_order
        self.trotter_reps = trotter_reps

        step = build_trotter_step(self._Hop, self.dt, trotter_order, trotter_reps)
        self._ustep = step.to_gate(label="U")           # uncontrolled Trotter step
        self._cstep = self._ustep.control(1)            # ancilla-controlled Trotter step
        self._observables = self._ancilla_observables()
        self._estimator = self._build_estimator(noise_model, seed)
        self._noise_floor = 5.0 / np.sqrt(float(shots)) if shots else 0.0

    # -- setup ---------------------------------------------------------------
    def _default_dt(self) -> float:
        try:
            e_max = float(eigsh(self._Hsparse, k=1, which="LA", return_eigenvectors=False)[0])
            e_min = float(eigsh(self._Hsparse, k=1, which="SA", return_eigenvectors=False)[0])
            width = e_max - e_min
        except Exception:
            width = 2.0 * float(np.sum(np.abs(self._Hop.coeffs)))
        if not np.isfinite(width) or width <= 0:
            width = 2.0 * float(np.sum(np.abs(self._Hop.coeffs)))
        return np.pi / width

    def _ancilla_observables(self) -> List[SparsePauliOp]:
        """[Xₐ⊗I, Yₐ⊗I, Xₐ⊗H, Yₐ⊗H] on N+1 qubits (ancilla = highest index)."""
        n = self._N
        labels = self._Hop.paulis.to_labels()
        coeffs = self._Hop.coeffs
        return [
            SparsePauliOp(["X" + "I" * n], [1.0]),
            SparsePauliOp(["Y" + "I" * n], [1.0]),
            SparsePauliOp(["X" + lbl for lbl in labels], coeffs),
            SparsePauliOp(["Y" + lbl for lbl in labels], coeffs),
        ]

    def _build_estimator(self, noise_model, seed):
        from qiskit_aer.primitives import EstimatorV2

        options = {}
        if noise_model is not None:
            options["backend_options"] = {"noise_model": noise_model}
        if self.shots is not None:
            options["default_precision"] = 1.0 / np.sqrt(float(self.shots))
        if seed is not None:
            options["run_options"] = {"seed_simulator": int(seed)}
        return EstimatorV2(options=options) if options else EstimatorV2()

    # -- measurement ---------------------------------------------------------
    def _pair_circuit(self, i: int, j: int) -> QuantumCircuit:
        n = self._N
        qc = QuantumCircuit(n + 1)
        qc.compose(self.mh.hf_circuit, qubits=range(n), inplace=True)
        for _ in range(i):
            qc.append(self._ustep, range(n))            # U^i (uncontrolled)
        qc.h(n)                                          # ancilla -> |+>
        for _ in range(j - i):
            qc.append(self._cstep, [n] + list(range(n)))  # controlled-U^(j-i)
        return transpile(qc, basis_gates=_BASIS_GATES, optimization_level=0)

    def _measure_observables(self, circuit) -> np.ndarray:
        """Expectation values of [Xₐ⊗I, Yₐ⊗I, Xₐ⊗H, Yₐ⊗H] on one circuit."""
        return np.asarray(
            self._estimator.run([(circuit, self._observables)]).result()[0].data.evs, dtype=float
        )

    def _measure_pair(self, i: int, j: int):
        """Measure (Sᵢⱼ, Hᵢⱼ), applying zero-noise extrapolation if configured."""
        base = self._pair_circuit(i, j)
        if not self.zne_scale_factors:
            evs = self._measure_observables(base)
        else:
            # Global folding C -> C (C^dagger C)^n amplifies noise; extrapolate each observable to 0.
            samples = np.array([
                self._measure_observables(fold_global_circuit(base, s))
                for s in self.zne_scale_factors
            ])  # shape (n_scales, 4)
            evs = np.array([
                extrapolate_zero_noise(self.zne_scale_factors, samples[:, k], self.zne_order)
                for k in range(4)
            ])
        return evs[0] + 1j * evs[1], evs[2] + 1j * evs[3]

    def _measure_matrices(self, dim: int):
        S = np.zeros((dim, dim), dtype=complex)
        H = np.zeros((dim, dim), dtype=complex)
        for i in range(dim):
            for j in range(i, dim):
                S[i, j], H[i, j] = self._measure_pair(i, j)
                S[j, i] = np.conj(S[i, j])
                H[j, i] = np.conj(H[i, j])
        return 0.5 * (H + H.conj().T), 0.5 * (S + S.conj().T)

    # -- public API ----------------------------------------------------------
    def solve(self, krylov_dim: int) -> KrylovStep:
        if krylov_dim < 1:
            raise ValueError("krylov_dim must be >= 1")
        H, S = self._measure_matrices(krylov_dim)
        energy, rank = solve_generalized_eig(H, S, self.threshold, self._noise_floor)
        return KrylovStep(dim=krylov_dim, rank=rank, energy=energy + self.offset)

    def resource_report(self, krylov_dim: int, shots: Optional[int] = None) -> dict:
        """Estimate the quantum resources to run the solver at dimension ``krylov_dim``.

        Reports the qubit count (system + 1 ancilla), the Trotter step depth / 2-qubit-gate cost,
        the depth / 2-qubit cost of the deepest Hadamard-test circuit (the i=0, j=M−1 pair, with
        the most controlled-U powers), the number of distinct pair circuits (M(M+1)/2), the total
        observable evaluations (4 per circuit × ZNE factor), and the total shot budget.
        """
        shots = shots if shots is not None else self.shots
        n_zne = len(self.zne_scale_factors) if self.zne_scale_factors else 1

        step = transpile(
            build_trotter_step(self._Hop, self.dt, self.trotter_order, self.trotter_reps),
            basis_gates=_BASIS_GATES, optimization_level=0,
        )
        deepest = self._pair_circuit(0, krylov_dim - 1)
        if self.zne_scale_factors and max(self.zne_scale_factors) > 1:
            deepest = fold_global_circuit(deepest, max(self.zne_scale_factors))
        deep_ops = deepest.count_ops()

        n_pairs = krylov_dim * (krylov_dim + 1) // 2
        observable_evals = n_pairs * 4 * n_zne
        return {
            "qubits": self._N + 1,
            "krylov_dim": krylov_dim,
            "hamiltonian_pauli_terms": len(self._Hop),
            "trotter_step_depth": step.depth(),
            "trotter_step_cx": step.count_ops().get("cx", 0),
            "deepest_circuit_depth": deepest.depth(),
            "deepest_circuit_cx": deep_ops.get("cx", 0),
            "distinct_pair_circuits": n_pairs,
            "zne_factor": n_zne,
            "observable_evaluations": observable_evals,
            "shots_per_evaluation": shots,
            "total_shots": (shots or 0) * observable_evals,
        }

    def convergence(self, max_dim: int) -> List[KrylovStep]:
        """Measure the full max_dim x max_dim problem once, then solve every leading block."""
        H, S = self._measure_matrices(max_dim)
        out = []
        for d in range(1, max_dim + 1):
            energy, rank = solve_generalized_eig(H[:d, :d], S[:d, :d], self.threshold, self._noise_floor)
            out.append(KrylovStep(dim=d, rank=rank, energy=energy + self.offset))
        return out


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")
    fci = mh.ground_state_energy()
    print(f"H2 FCI = {fci:.6f} Ha  (subspace matrices measured by Hadamard tests on Aer)")
    for s in HardwareKrylovSolver(mh).convergence(4):
        print(f"  M={s.dim} rank={s.rank} E={s.energy:.6f}  err={(s.energy - fci) * 1e3:+.4f} mHa")
