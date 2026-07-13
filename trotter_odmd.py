#!/usr/bin/env python3
"""
Circuit-real ODMD -- Trotter eigenphases, the dt^2 law, and Richardson bias removal.

A fixed-step Suzuki-Trotter circuit is EXACT evolution of an effective Hamiltonian
H_eff = H + O(dt_eff^2), so ODMD on the circuit signal s_k = <phi0| U_trot^k |phi0> returns the
ground eigenphase of the circuit unitary to machine precision -- DMD stacks no approximation of
its own on top of Trotterization. The eigenphase bias vs FCI then follows the clean second-order
law (ratio ~4 per reps doubling), and two-point Richardson extrapolation across reps,

    E = (r^p E_fine - E_coarse) / (r^p - 1),    r = step ratio, p = Trotter order,

removes it below 0.1 mHa. Prior art: Trotter effective-Hamiltonian theory (arXiv:1912.08854);
extrapolation of Trotterized phase estimation (e.g. arXiv:2212.14144). Reproduction of known
theory, composed with the pinned ODMD machinery.

Found while probing this spec (see specs/SPEC_trotter_odmd.md section 2): the repo's Trotter
statevector path was silently EXACT -- Operator()/Statevector.evolve() ignore an opaque
PauliEvolutionGate's synthesis. ``build_trotter_step`` now materializes the synthesized circuit;
this module's problem builder measures the resulting unitary deviation explicitly
(``unitary_deviation``) so the Trotterization can never again be silently absent.

HONEST SCOPE: statevector-simulated circuits (no device noise / native-gate transpilation);
Richardson's price is circuit DEPTH (reps=4 is 4x deeper), not shots; extrapolate from the fine
pair and only when bias > noise (large-dt pairs leave higher-order residue -- H2 reps 1+2 leaves
0.9 mHa).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from qiskit.quantum_info import Operator, SparsePauliOp, Statevector
from scipy.linalg import expm

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.trotter_krylov import build_trotter_step
from odmd import ODMDProblem


@dataclass
class TrotterODMDProblem(ODMDProblem):
    """ODMD problem whose signal comes from a genuinely Trotterized circuit.

    Subclasses :class:`odmd.ODMDProblem`, so ``odmd.odmd_energy`` / ``odmd.sample_odmd_energy``
    apply unchanged. ``ref`` stays the exact (FCI) centered ground energy; ``e_circuit`` is the
    exact ground eigenphase of the circuit unitary -- the value ODMD actually estimates, whose
    difference from ``ref`` is the Trotter bias.
    """
    reps: int = 1
    order: int = 2
    e_circuit: float = 0.0        # ground eigenphase of U_trot (centered frame, HF-reachable)
    unitary_deviation: float = 0.0  # ||U_trot - exp(-i tau H)||_2  (0 would mean: not Trotterized)
    depth: int = 0                # circuit depth of one step


def select_ground_eigenphase(angles: np.ndarray, pops_u: np.ndarray, width: float,
                             pop_cut: float = 1e-8) -> float:
    """Pick the ground (minimum) eigenphase among HF-reachable circuit eigenvalues.

    ``angles`` are ``-angle(lam)/tau`` for every eigenvalue of the (centered) circuit unitary;
    the physical eigenphases live in the closed band ``[-width/2, width/2]`` by construction
    (``tau = pi/width``) -- a periodic image outside it is an artifact of the ``arg`` branch cut,
    not a reachable energy, and must never be chosen even if its population clears ``pop_cut``.
    A small relative tolerance guards the closed boundary against roundoff in the eig/angle
    chain: without it, a candidate that IS the band edge can round a few ULPs outside and get
    dropped, silently flipping which branch is selected.
    """
    band = width / 2 + 1e-9 * max(width, 1.0)
    physical = (pops_u > pop_cut) & (np.abs(angles) <= band)
    return float(np.min(angles[physical]))


def build_trotter_odmd_problem(mh: MolecularHamiltonian, n: int = 24, reps: int = 1,
                               order: int = 2) -> TrotterODMDProblem:
    """Trotterized survival amplitudes + exact circuit-eigenphase references.

    Centered frame identical to ``odmd.build_odmd_problem`` -- the shift by mu is applied to the
    Hamiltonian BEFORE synthesis (an identity term only adds a global phase; it does not change
    the splitting error). References computed exactly: dense diagonalization for FCI, dense
    ``eig`` of the circuit ``Operator`` for the eigenphase, an explicit 2-norm against
    ``expm(-i tau H)`` for the deviation.
    """
    nq = mh.num_qubits
    H_dense = np.asarray(mh.qubit_hamiltonian.to_matrix())
    psi0 = np.asarray(mh.hf_state().data, dtype=complex)
    w_eig, V = np.linalg.eigh(H_dense)
    pops = np.abs(V.conj().T @ psi0) ** 2
    reach = w_eig[pops > 1e-8].real
    width = float(reach.max() - reach.min())
    mu = float(0.5 * (reach.max() + reach.min()))
    tau = float(np.pi / width)

    H_s = (mh.qubit_hamiltonian - SparsePauliOp("I" * nq, coeffs=[mu])).simplify()
    step = build_trotter_step(H_s, tau, order=order, reps=reps)
    U = Operator(step).data
    U_exact = expm(-1j * tau * (H_dense - mu * np.eye(H_dense.shape[0])))
    deviation = float(np.linalg.norm(U - U_exact, 2))
    lam, W = np.linalg.eig(U)
    pops_u = np.abs(W.conj().T @ psi0) ** 2
    angles = -np.angle(lam) / tau
    e_circuit = select_ground_eigenphase(angles, pops_u, width)

    s = np.empty(n, dtype=complex)
    s[0], psi = 1.0, Statevector(psi0)
    for k in range(1, n):
        psi = psi.evolve(step)
        s[k] = psi0.conj() @ psi.data
    return TrotterODMDProblem(n=n, tau=tau, mu=mu, s=s, dim=int(H_dense.shape[0]),
                              offset=mh.energy_offset + mu, ref=float(reach[0] - mu),
                              reps=reps, order=order, e_circuit=e_circuit,
                              unitary_deviation=deviation, depth=step.depth())


def richardson_energy(e_coarse: float, e_fine: float, step_ratio: float = 2.0,
                      order: int = 2) -> float:
    """Two-point Richardson extrapolation of Trotter eigenphases to dt_eff -> 0.

    ``e_coarse``/``e_fine`` from step sizes differing by ``step_ratio`` (fine = coarse /
    step_ratio), with leading bias ~ dt_eff^order.
    """
    w = float(step_ratio) ** order
    return float((w * e_fine - e_coarse) / (w - 1.0))


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
    from odmd import odmd_energy

    mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7")
    print("H4 chain, order-2 Suzuki, tau = pi/W  (bias = circuit eigenphase - FCI):")
    probs = {r: build_trotter_odmd_problem(mh, n=24, reps=r) for r in (1, 2, 4)}
    for r, p in probs.items():
        e, _ = odmd_energy(p.s, p.tau)
        print(f"  reps={r}: ||U-Uexact||={p.unitary_deviation:.2e}  depth={p.depth:5d}  "
              f"|E_odmd-E_U|={abs(e - p.e_circuit):.1e}  bias={(p.e_circuit - p.ref) * 1e3:+8.4f} mHa")
    for c, f in ((1, 2), (2, 4)):
        e_rich = richardson_energy(probs[c].e_circuit, probs[f].e_circuit)
        print(f"  Richardson reps {c}+{f}: residual = {(e_rich - probs[f].ref) * 1e3:+8.5f} mHa")
