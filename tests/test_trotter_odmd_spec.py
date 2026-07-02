"""
Acceptance gates G1-G4 for specs/SPEC_trotter_odmd.md (circuit-real ODMD + Richardson).

Test-first: ``trotter_odmd`` does not exist yet, so this file is RED until the spec is
implemented -- and G1 is additionally RED against the UNFIXED ``build_trotter_step``: probing
this spec exposed that Operator()/Statevector.evolve() evaluate an opaque PauliEvolutionGate via
its exact matrix, silently ignoring the SuzukiTrotter synthesis, so TrotterKrylovSolver had been
doing exact evolution all along (trotter_order/reps were no-ops; the old "within Trotter error"
test passed vacuously). The fix materializes the synthesized definition; G1 is its regression
gate (pre-fix deviation ~1e-16 vs the gated > 0.05).

Claims: ODMD on the genuinely-Trotterized signal returns the ground eigenphase of the circuit
unitary to machine precision (DMD adds no approximation of its own); that eigenphase's bias vs
FCI follows the second-order dt^2 law (ratio ~4 per reps doubling); two-point Richardson removes
it below 0.1 mHa noiselessly and beats the plain estimate >3x under shot noise. References:
dense diagonalization of the circuit Operator (exact) and FCI. PySCF/qiskit, no block2;
`make gates` runs it in its own process.
"""
import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from odmd import odmd_energy, sample_odmd_energy
from trotter_odmd import build_trotter_odmd_problem, richardson_energy

SYSTEMS = {
    "h2": dict(atom="H 0 0 0; H 0 0 0.74"),
    "h4": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
}
REPS = (1, 2, 4)
_CACHE = {}


def _prob(key, reps):
    if (key, reps) not in _CACHE:
        mh = _CACHE.setdefault(key, build_molecular_hamiltonian(**SYSTEMS[key]))
        _CACHE[(key, reps)] = build_trotter_odmd_problem(mh, n=24, reps=reps)
    return _CACHE[(key, reps)]


def test_G1_step_is_genuinely_trotterized():
    """REGRESSION GATE for the build_trotter_step fix: the circuit unitary deviates from exact
    evolution at reps=1 (pre-fix: ~1e-16) and the deviation shrinks ~4x per reps doubling."""
    floor = {"h2": 0.05, "h4": 0.02}
    for key in SYSTEMS:
        dev = {r: _prob(key, r).unitary_deviation for r in REPS}
        assert dev[1] > floor[key], (key, dev[1])
        for coarse, fine in ((1, 2), (2, 4)):
            ratio = dev[coarse] / dev[fine]
            assert 3.5 < ratio < 5.5, (key, coarse, fine, ratio)


def test_G2_dmd_adds_nothing():
    """ODMD on the circuit signal == ground eigenphase of U_trot to < 1e-9 Ha (measured 5e-11):
    the ONLY systematic error left is the effective-Hamiltonian shift."""
    for key in SYSTEMS:
        for r in REPS:
            prob = _prob(key, r)
            e, _ = odmd_energy(prob.s, prob.tau)
            assert abs(e - prob.e_circuit) < 1e-9, (key, r, e - prob.e_circuit)


def test_G3_dt_squared_eigenphase_law():
    """Bias(reps) = E_U - FCI scales as dt_eff^2: ratio in [3.3, 5.5] per reps doubling."""
    for key in SYSTEMS:
        bias = {r: _prob(key, r).e_circuit - _prob(key, r).ref for r in REPS}
        assert abs(bias[1]) > 1e-3, (key, bias[1])          # visible bias to extrapolate away
        for coarse, fine in ((1, 2), (2, 4)):
            ratio = bias[coarse] / bias[fine]
            assert 3.3 < ratio < 5.5, (key, coarse, fine, ratio)


def test_G4_richardson_removes_the_bias():
    """DEFINITION OF DONE: reps-(2,4) Richardson residual < 0.1 mHa and >= 10x below the reps=4
    bias (noiseless); under shot noise (H4, 1e6 shots/element) the extrapolated median error
    beats the plain reps=2 estimate by > 3x."""
    for key in SYSTEMS:
        p2, p4 = _prob(key, 2), _prob(key, 4)
        e_rich = richardson_energy(p2.e_circuit, p4.e_circuit)
        resid, bias4 = abs(e_rich - p4.ref), abs(p4.e_circuit - p4.ref)
        assert resid < 1e-4, (key, resid)
        assert resid < bias4 / 10.0, (key, resid, bias4)
    p1, p2 = _prob("h4", 1), _prob("h4", 2)
    shots, plain, rich = 1_000_000, [], []
    for sd in range(60):
        e1 = sample_odmd_energy(p1, shots, sd)
        e2 = sample_odmd_energy(p2, shots, sd + 50_000)     # independent noise draws
        plain.append(abs(e2 - p2.ref))
        rich.append(abs(richardson_energy(e1, e2) - p2.ref))
    assert np.median(rich) < np.median(plain) / 3.0, (np.median(rich), np.median(plain))
