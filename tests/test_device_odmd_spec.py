"""
Acceptance gates G1-G4 for specs/SPEC_device_odmd.md (device-noise ODMD).

Test-first: ``device_odmd`` does not exist yet, so this file is RED until the spec is
implemented. Claim: a global depolarizing channel damps the survival signal s_k -> f^k s_k,
multiplying every DMD eigenvalue by f but leaving its PHASE untouched -- so the ODMD energy is
exactly damping-invariant while KQD's GEVP on the same damped data is not. Local gate-level
noise (Aer NoiseModel on real transpiled Hadamard-test circuits) is NOT a global channel, so the
residual phase bias is measured, not assumed zero; immunity ends when damping pushes the signal
under the shot-noise floor (the SPEC_odmd_excited noise-edge law).

References: exact centered ground energy (channel gates, H4) and the exact ground eigenphase of
the SAME Trotter step circuit (Aer gates, H2 -- isolating noise-induced bias from the already
spec'd Trotter bias). Medians over noise seeds. PySCF/qiskit-aer, no block2; `make gates` runs
it in its own process.
"""
import numpy as np
import scipy.sparse as sp
from scipy.linalg import toeplitz
from scipy.sparse.linalg import expm_multiply

from device_odmd import device_odmd_energy, measure_survival_signal
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.noise import build_depolarizing_noise_model
from hybrid_quantum_solver.quantum_krylov_solver import solve_generalized_eig
from odmd import build_odmd_problem, odmd_energy
from trotter_odmd import build_trotter_odmd_problem

_CACHE = {}


def _h4():
    """H4 channel-level fixture: (ODMDProblem, exact KQD H-row, per-element sigma at 1e5 shots)."""
    if "h4" not in _CACHE:
        mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7")
        prob = build_odmd_problem(mh, n=24)
        Hm = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
        Hs = (Hm - prob.mu * sp.identity(Hm.shape[0], format="csc")).tocsc()
        psi0 = np.asarray(mh.hf_state().data, dtype=complex)
        h_row = np.array([psi0.conj() @ (Hs @ expm_multiply(-1j * (k * prob.tau) * Hs, psi0))
                          for k in range(prob.n)])
        sigma = np.sqrt(2.0 * (2.0 - 1.0 / prob.dim) / 1e5)
        _CACHE["h4"] = (prob, h_row, sigma)
    return _CACHE["h4"]


def _kqd(h_vec, s_vec, floor):
    return solve_generalized_eig(toeplitz(np.conj(h_vec), h_vec),
                                 toeplitz(np.conj(s_vec), s_vec), 1e-12, floor)[0]


def _damped_noisy_medians(f, seeds=100):
    prob, h_row, sigma = _h4()
    damp = f ** np.arange(prob.n)
    e_dev, e_def, e_kqd = [], [], []
    for sd in range(seeds):
        rng = np.random.default_rng(sd)
        g = rng.normal(0, sigma / np.sqrt(2), prob.n) + 1j * rng.normal(0, sigma / np.sqrt(2), prob.n)
        g[0] = 0.0
        s_n = prob.s * damp + g
        gh = rng.normal(0, sigma / np.sqrt(2), prob.n) + 1j * rng.normal(0, sigma / np.sqrt(2), prob.n)
        h_n = h_row * damp + gh
        e_dev.append(abs(device_odmd_energy(s_n, prob.tau, sigma, amp_floor=0.02) - prob.ref))
        e_def.append(abs(odmd_energy(s_n, prob.tau, svd_threshold=5 * sigma)[0] - prob.ref))
        e_kqd.append(abs(_kqd(h_n, s_n, 5 * sigma) - prob.ref))
    return tuple(float(np.median(x)) for x in (e_dev, e_def, e_kqd))


def test_G1_exact_channel_immunity():
    """Uniform damping leaves ODMD eigenphases EXACT (< 1e-6 Ha) while KQD on identically damped
    rows drifts (> 0.5 / > 1 mHa at f = 0.9 / 0.7)."""
    prob, h_row, _ = _h4()
    for f, kqd_floor in ((0.9, 5e-4), (0.7, 1e-3)):
        damp = f ** np.arange(prob.n)
        e_dev = device_odmd_energy(prob.s * damp, prob.tau, 0.0, amp_floor=0.02)
        assert abs(e_dev - prob.ref) < 1e-6, (f, e_dev - prob.ref)
        assert abs(_kqd(h_row * damp, prob.s * damp, 0.0) - prob.ref) > kqd_floor, f


def test_G2_damped_plus_shot_noise_beats_kqd():
    """DEFINITION OF DONE: damping + 1e5-shot noise on matched data -- ODMD stays at the mHa
    level while KQD fails by orders of magnitude."""
    for f, odmd_tol, min_ratio in ((0.9, 1e-3, 100.0), (0.7, 5e-3, 50.0)):
        m_dev, _, m_kqd = _damped_noisy_medians(f)
        assert m_dev < odmd_tol, (f, m_dev)
        assert m_kqd / m_dev > min_ratio, (f, m_kqd / m_dev)


def test_G3_modulus_window_is_the_mechanism():
    """The unit-modulus filter of odmd_energy misidentifies damped modes under noise (>= 2x worse
    than the wide-window device estimator at f=0.9); at f=1 both are exact."""
    prob, _, _ = _h4()
    assert abs(device_odmd_energy(prob.s, prob.tau, 0.0, amp_floor=0.02) - prob.ref) < 1e-8
    assert abs(odmd_energy(prob.s, prob.tau)[0] - prob.ref) < 1e-8
    m_dev, m_def, _ = _damped_noisy_medians(0.9)
    assert m_def / m_dev > 2.0, (m_def, m_dev)


def test_G4_aer_phase_survival_and_boundary():
    """Real transpiled Hadamard-test circuits under an Aer depolarizing NoiseModel: (a) faithful
    at zero noise; (b) < 0.5 mHa phase error through >= 50% amplitude loss at cx=3e-4;
    (c) at cx=1e-3 the signal falls under the shot floor and immunity ends (> 1 mHa)."""
    mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")
    ref = build_trotter_odmd_problem(mh, n=8, reps=1)      # same frame, same step circuit
    K, shots, sigma = 8, 32768, 1.0 / np.sqrt(32768.0)

    def run(noise_model):
        errs, d7 = [], []
        for seed in range(5):
            s = measure_survival_signal(mh, K, shots=shots, noise_model=noise_model, seed=seed)
            errs.append(abs(device_odmd_energy(s, ref.tau, sigma) - ref.e_circuit))
            d7.append(abs(s[7]) / abs(ref.s[7]))
        return float(np.median(errs)), float(np.median(d7))

    err0, d0 = run(None)
    assert err0 < 1e-5, err0                               # (a) faithful plumbing
    assert d0 > 0.95, d0
    err3, d3 = run(build_depolarizing_noise_model(3e-5, 3e-4, 3e-4))
    assert d3 < 0.5, d3                                    # (b) >= 50% amplitude lost...
    assert err3 < 5e-4, err3                               #     ...yet the phase survives
    err10, d10 = run(build_depolarizing_noise_model(1e-4, 1e-3, 1e-3))
    assert d10 < 0.05, d10                                 # (c) signal under the shot floor...
    assert err10 > 1e-3, err10                             #     ...immunity ends (the boundary)
