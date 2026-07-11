#!/usr/bin/env python3
"""
Trotter resolution floor -- when a circuit eigenphase cannot be trusted.

Grew out of the nb3x8_device_gap G2 coin-flip flake (`make gates` red ~60% of runs). The
assertion tested `circuit_gap(Nb3F8, reps=1)` -- a quantity that lives BELOW the Trotter
method's own resolution floor, and whose value therefore depended on which Pauli ordering
Python's hash randomization dealt to the Suzuki-Trotter synthesis that run.

THE LAW. A reference state |psi0> assigns each circuit eigenphase a population. Trotterization
leaks amplitude of order ||U_trot - U_exact||_2 from the DOMINANT eigenspaces into every other
one -- a spurious population of order dev^2 that interferes with the genuine signal. So an
eigenphase is *resolvable* only if the genuine population exceeds the leakage floor:

    resolvable  <=>  |<eig|psi0>|^2  >  ||U_trot - U_exact||_2^2

Below the floor, the computed population at that eigenphase is whatever interference leaves
(observed: anywhere from 1e-5 down to 2.5e-13 for a genuine 1.36e-5 signal), it can fall under
any sensible selection cut, and the extracted branch flips by the wrap quantum 2*pi/tau.

THE PROBE. Ordering is a deterministic perturbation of the same physics: re-synthesize under
different canonical term orderings and compare. A resolvable quantity moves by ordinary Trotter
bias (small, smooth in reps); an unresolvable one flips branch -- spread ~ 2*pi/tau. For Nb3F8:
reps=1 spread 3756.7 meV (the wrap quantum; below the floor), reps=2 spread 5.5 meV (above).

Numbers (Nb3F8 sector nelec=2, genuine population 1.36e-5): floor = 3.6e-5 at reps=1 (leak >
signal -> unresolvable), 1.3e-6 at reps=2, 7.5e-8 at reps=4. Exactly the one assertion that
coin-flipped was below the floor; every assertion that never flaked is above it.

HONEST SCOPE: dev is a global 2-norm, so dev^2 is a conservative leakage scale -- the criterion
may call "unresolvable" something that happens to resolve, never the reverse. Floor numbers here
are SuzukiTrotter(order=2); other orders rescale dev but the criterion is generic. The fix that
makes any of this reproducible is `canonical_term_order` inside `build_trotter_step`
(largest-|coeff| first) -- it also reproduces the historical green-run numbers exactly.

See specs/SPEC_trotter_resolution_floor.md.
"""
from __future__ import annotations

import numpy as np
from qiskit.quantum_info import Operator, SparsePauliOp
from scipy.linalg import expm

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.trotter_krylov import build_trotter_step


def _centered(mh: MolecularHamiltonian):
    """The centered frame of trotter_odmd.build_trotter_odmd_problem: (H_dense, psi0, mu, tau)."""
    H_dense = np.asarray(mh.qubit_hamiltonian.to_matrix())
    psi0 = np.asarray(mh.hf_state().data, dtype=complex)
    w_eig, V = np.linalg.eigh(H_dense)
    pops = np.abs(V.conj().T @ psi0) ** 2
    reach = w_eig[pops > 1e-8].real
    mu = float(0.5 * (reach.max() + reach.min()))
    tau = float(np.pi / (reach.max() - reach.min()))
    return H_dense, psi0, mu, tau, w_eig, V, pops


def reference_population(mh: MolecularHamiltonian) -> float:
    """|<lowest reachable eigenstate|psi0>|^2 against the EXACT Hamiltonian -- the genuine
    signal at the eigenphase the pipeline extracts (`min` over populations > 1e-8, the same
    reachable set that defines its tau/mu frame), uncorrupted by Trotter leakage. The GLOBAL
    ground state can be strictly unreachable (population 0 by symmetry) -- it never anchors
    the eigenphase, so it is not the relevant signal."""
    _, _, _, _, w_eig, _, pops = _centered(mh)
    reachable = pops > 1e-8
    i_low = int(np.argmin(np.where(reachable, w_eig.real, np.inf)))
    return float(pops[i_low])


def leakage_floor(mh: MolecularHamiltonian, reps: int, order: int = 2) -> float:
    """The spurious-population scale at any eigenphase: ||U_trot - U_exact||_2^2 in the same
    centered frame the ODMD pipeline uses."""
    H_dense, _, mu, tau, _, _, _ = _centered(mh)
    nq = mh.num_qubits
    H_s = (mh.qubit_hamiltonian - SparsePauliOp("I" * nq, coeffs=[mu])).simplify()
    U = Operator(build_trotter_step(H_s, tau, order=order, reps=reps)).data
    U_exact = expm(-1j * tau * (H_dense - mu * np.eye(H_dense.shape[0])))
    return float(np.linalg.norm(U - U_exact, 2)) ** 2


def is_resolvable(mh: MolecularHamiltonian, reps: int, order: int = 2) -> bool:
    """Can the ground eigenphase of this sector be trusted at this Trotter depth?"""
    return reference_population(mh) > leakage_floor(mh, reps, order=order)


# --- the ordering probe ---------------------------------------------------------------------

_ORDERINGS = {
    "coeff-desc": lambda label, coeff: (-abs(coeff), label),
    "coeff-asc": lambda label, coeff: (abs(coeff), label),
    "label-asc": lambda label, coeff: label,
}


def _reordered(op: SparsePauliOp, key) -> SparsePauliOp:
    labels = [p.to_label() for p in op.paulis]
    coeffs = [complex(c) for c in op.coeffs]
    idx = sorted(range(len(labels)), key=lambda i: key(labels[i], coeffs[i]))
    return SparsePauliOp.from_list([(labels[i], coeffs[i]) for i in idx])


def _eigenphase_energy(mh: MolecularHamiltonian, reps: int, key, order: int = 2) -> float:
    """circuit eigenphase ground energy under an explicit term ordering (bypasses the canonical
    sort by synthesizing the reordered operator directly -- the probe's whole point)."""
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import PauliEvolutionGate
    from qiskit.synthesis import SuzukiTrotter

    _, psi0, mu, tau, _, _, _ = _centered(mh)
    nq = mh.num_qubits
    H_s = _reordered((mh.qubit_hamiltonian - SparsePauliOp("I" * nq, coeffs=[mu])).simplify(), key)
    gate = PauliEvolutionGate(H_s, time=tau, synthesis=SuzukiTrotter(order=order, reps=reps))
    qc = QuantumCircuit(nq)
    qc.compose(gate.definition, inplace=True)
    U = Operator(qc).data
    lam, W = np.linalg.eig(U)
    pops_u = np.abs(W.conj().T @ psi0) ** 2
    e = float(np.min(-np.angle(lam[pops_u > 1e-8]) / tau))
    return e + mh.energy_offset + mu


def ordering_spread(params: dict, reps: int, order: int = 2) -> float:
    """The probe, on an Nb3X8 material's charge gap: max-min of the gap over the canonical
    orderings, in meV. Branch flips show up as ~ the 2*pi/tau wrap quantum; resolvable gaps
    spread only by ordinary Trotter bias.

    ``params`` is an NB3X8_LT_BULK-style dict (U0, t, Us).
    """
    from nb3x8_device_gap import SECTOR_WEIGHTS, sector_models

    models = sector_models(**params)
    gaps = []
    for key in _ORDERINGS.values():
        gap = 0.0
        for nelec, weight in SECTOR_WEIGHTS:
            gap += weight * _eigenphase_energy(models[nelec].to_hamiltonian(), reps, key,
                                               order=order)
        gaps.append(gap)
    return float(max(gaps) - min(gaps))


if __name__ == "__main__":
    from nb3x8_device_gap import sector_models
    from nb3x8_gaps import NB3X8_LT_BULK

    print("Trotter resolution floor -- Nb3F8 sector nelec=2 (the flake's origin)")
    mh = sector_models(**NB3X8_LT_BULK["Nb3F8"])[2].to_hamiltonian()
    pop = reference_population(mh)
    print(f"genuine ground population = {pop:.3e}")
    for reps in (1, 2, 4):
        floor = leakage_floor(mh, reps)
        print(f"reps={reps}: floor = {floor:.3e}  -> resolvable: {pop > floor}")
    for reps in (1, 2):
        print(f"ordering spread of the F8 gap at reps={reps}: "
              f"{ordering_spread(NB3X8_LT_BULK['Nb3F8'], reps):.1f} meV")
