"""
Acceptance gates G1-G4 for specs/SPEC_qksd_noise.md (excited-state QKSD under shot noise).

No new code: this validates the existing ``solve_excited`` (SPEC_qksd_excited) under the finite-
sampling shot-noise model already in ``QuantumKrylovSolver`` (Hermitian-symmetric perturbation +
noise-aware overlap cutoff). Claim: excited energies and the first gap degrade gracefully with the
shot budget, but the excited state is much more noise-fragile than the ground state. The noise is
random, so every gate averages over seeds.

Uses only pyscf/qiskit (no block2); `make gates` runs it in its own process.
"""
import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.noise import shot_noise_sigma
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

KDIM = 10
SEEDS = 12


def _h2():
    return build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")


def _reachable(mh, overlap_tol=1e-8):
    w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    hf = np.asarray(mh.hf_state().data, dtype=complex)
    overlaps = np.abs(V.conj().T @ hf) ** 2
    return np.sort(w[overlaps > overlap_tol].real) + mh.energy_offset


def _excited_errors(mh, ref0, ref1, shots, seeds=SEEDS):
    """Mean |dE0|, |dE1|, |dgap|, and the min excited energy over seeds at a shot budget."""
    sigma = shot_noise_sigma(shots)
    e0, e1, gap, emin = [], [], [], np.inf
    for s in range(seeds):
        en = np.sort(QuantumKrylovSolver(mh, noise_sigma=sigma, seed=s).solve_excited(KDIM, 2).energies)
        e0.append(abs(en[0] - ref0))
        if len(en) >= 2:                              # excited state resolved this seed
            e1.append(abs(en[1] - ref1))
            gap.append(abs((en[1] - en[0]) - (ref1 - ref0)))
            emin = min(emin, en[1])
    return np.mean(e0), np.mean(e1), np.mean(gap), emin


def test_G1_noiseless_excited_recovery():
    """sigma=0 anchor: the excited energy and the gap match the dense reachable reference."""
    mh = _h2()
    ref = _reachable(mh)
    en = np.sort(QuantumKrylovSolver(mh).solve_excited(KDIM, 2).energies)
    assert abs(en[1] - ref[1]) < 1e-6, (en[1], ref[1])
    assert abs((en[1] - en[0]) - (ref[1] - ref[0])) < 1e-6, (en[1] - en[0], ref[1] - ref[0])


def test_G2_bounded_under_noise():
    """At a modest budget the excited error is bounded and never blows up."""
    mh = _h2()
    ref = _reachable(mh)
    _, me1, _, emin = _excited_errors(mh, ref[0], ref[1], shots=4096)
    assert me1 < 0.1, me1                              # chemical-to-mHa scale, not -800 Ha pathology
    assert emin > ref[1] - 0.5, (emin, ref[1])         # controlled: no runaway below the target


def test_G3_gap_improves_with_shots():
    """DEFINITION OF DONE: more shots -> smaller mean first-excitation-gap error."""
    mh = _h2()
    ref = _reachable(mh)
    _, _, gap_low, _ = _excited_errors(mh, ref[0], ref[1], shots=4096)
    _, _, gap_high, _ = _excited_errors(mh, ref[0], ref[1], shots=262144)
    assert gap_high < gap_low, (gap_high, gap_low)


def test_G4_excited_more_fragile_than_ground():
    """THE FINDING: at a fixed budget the excited error dwarfs the ground error (measured ~24-36x)."""
    mh = _h2()
    ref = _reachable(mh)
    me0, me1, _, _ = _excited_errors(mh, ref[0], ref[1], shots=16384)
    assert me1 > 5.0 * me0, (me1, me0, me1 / me0)
