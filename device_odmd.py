#!/usr/bin/env python3
"""
Device-noise ODMD -- eigenphases are depolarizing-immune, until the noise edge.

A global depolarizing channel damps the survival signal geometrically, s_k -> f^k s_k. That
multiplies every DMD eigenvalue by f and leaves its PHASE untouched, so the ODMD energy is
exactly invariant under uniform damping -- while KQD's generalized eigenproblem on the same
damped Toeplitz data has no such protection (the damped matrices are not a similarity transform
of the originals). Local gate-level noise on real circuits is NOT a global channel, so the
residual phase bias must be measured, not assumed zero: the Aer arm here drives the full
hardware stack (transpiled ancilla-controlled Trotter circuits, `hardware_krylov` Hadamard
tests, a qiskit-aer NoiseModel) and references the exact eigenphase of the same step circuit.

Two compositions, no new DMD code (see specs/SPEC_device_odmd.md):
  * ``device_odmd_energy`` -- the damping-robust estimator: noise-edge SVD cutoff + a WIDE
    modulus window (odmd_energy's unit-modulus filter ||lam|-1|<0.2 misidentifies damped signal
    modes at f < 0.8; under noise, spurious near-unimodular modes pass instead) + an amplitude
    floor from the Vandermonde refit.
  * ``measure_survival_signal`` -- the centered-frame Hadamard-test measurement of s_k = S_0k.

HONEST SCOPE: Aer depolarizing + readout is a simulated device (no coherent errors, crosstalk,
drift); the measured energies are eigenphases of the TROTTERIZED unitary (quote against the
circuit eigenphase, or remove the bias per SPEC_trotter_odmd); amplitude floors hide states with
p_n below them (the visibility law); ODMD remains non-variational.
"""
from __future__ import annotations

import dataclasses

import numpy as np
from qiskit.quantum_info import SparsePauliOp

from hybrid_quantum_solver.hardware_krylov import HardwareKrylovSolver
from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from odmd import noise_edge, odmd_spectrum


def centered_frame(mh: MolecularHamiltonian):
    """(mh_centered, tau, mu): shift the qubit Hamiltonian by the HF-reachable spectral center.

    Total energies are invariant (the shift moves into ``energy_offset``); tau = pi / W keeps
    every reachable eigenphase inside (-pi/2, pi/2) -- no wrapping. Same frame conventions as
    ``odmd.build_odmd_problem`` / ``trotter_odmd.build_trotter_odmd_problem`` (dense
    diagonalization: validation scale only).
    """
    w_eig, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    psi0 = np.asarray(mh.hf_state().data, dtype=complex)
    pops = np.abs(V.conj().T @ psi0) ** 2
    reach = w_eig[pops > 1e-8].real
    mu = float(0.5 * (reach.max() + reach.min()))
    tau = float(np.pi / (reach.max() - reach.min()))
    shifted = (mh.qubit_hamiltonian
               - SparsePauliOp("I" * mh.num_qubits, coeffs=[mu])).simplify()
    mh_c = dataclasses.replace(mh, qubit_hamiltonian=shifted,
                               energy_offset=mh.energy_offset + mu)
    return mh_c, tau, mu


def device_odmd_energy(s, tau: float, sigma: float, amp_floor: float = 0.05,
                       c: float = 1.2) -> float:
    """Damping-robust ODMD ground-energy estimate (centered frame) from a device signal.

    ``sigma`` is the per-element measurement noise scale (~1/sqrt(shots); 0 = exact): it sets
    the absolute noise-edge SVD cutoff. The modulus window is wide open (2.0) because damped
    signal modes sit at |lam| = f < 1; mode selection instead uses the amplitude floor (modes
    are kept iff their Vandermonde-refit amplitude exceeds ``amp_floor``). Returns +inf if no
    mode survives -- signal below the noise floor, the visibility boundary.
    """
    n = len(s)
    d = n // 2
    cutoff = noise_edge(max(float(sigma), 1e-10), d, n - d, c)
    energies, _, _ = odmd_spectrum(s, tau, cutoff=cutoff, mod_window=2.0, amp_floor=amp_floor)
    return float(energies[0]) if len(energies) else np.inf


def measure_survival_signal(mh: MolecularHamiltonian, n: int, shots: int | None = None,
                            noise_model=None, seed: int | None = None, trotter_order: int = 2,
                            trotter_reps: int = 1) -> np.ndarray:
    """Measure s_k = <phi0|U_trot^k|phi0>, k = 0..n-1, by ancilla Hadamard tests on qiskit-aer.

    Centered frame; genuinely Trotterized, transpiled, ancilla-controlled circuits (the
    ``hardware_krylov`` stack). ``shots=None, noise_model=None`` -> exact statevector
    expectation values (plumbing reference).
    """
    mh_c, tau, _ = centered_frame(mh)
    solver = HardwareKrylovSolver(mh_c, dt=tau, shots=shots, noise_model=noise_model, seed=seed,
                                  trotter_order=trotter_order, trotter_reps=trotter_reps)
    return solver.measure_signal(n)


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
    from hybrid_quantum_solver.noise import build_depolarizing_noise_model
    from trotter_odmd import build_trotter_odmd_problem

    mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")
    ref = build_trotter_odmd_problem(mh, n=8, reps=1)
    shots = 32768
    print("H2, K=8 Hadamard-test signal (32768 shots/observable, 3 seeds/row);")
    print("errors vs the exact circuit eigenphase (Trotter bias excluded by construction):")
    for cx in (0.0, 1e-4, 3e-4, 1e-3):
        nm = None if cx == 0.0 else build_depolarizing_noise_model(cx / 10, cx, cx)
        errs, damp = [], []
        for seed in range(3):
            s = measure_survival_signal(mh, 8, shots=shots, noise_model=nm, seed=seed)
            errs.append(abs(device_odmd_energy(s, ref.tau, 1 / np.sqrt(shots)) - ref.e_circuit))
            damp.append(abs(s[7]) / abs(ref.s[7]))
        print(f"  cx={cx:7.1e}: |E - e_circ| = {np.median(errs) * 1e3:8.3f} mHa   "
              f"amplitude retained |s7|/|s7_ideal| = {np.median(damp):5.3f}")
