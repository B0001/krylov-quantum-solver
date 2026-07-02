"""
Acceptance gates G1-G4 for specs/SPEC_odmd_excited.md (excited-state ODMD via noise-edge
thresholding).

Test-first: ``odmd_spectrum`` does not exist yet, so this file is RED until the spec is
implemented. Claim: the SAME survival-amplitude signal ODMD uses for the ground state carries the
low-lying excited spectrum in its higher DMD eigenphases -- with NO extra measurements -- provided
the SVD truncation is an absolute noise-edge cutoff c*sigma*(sqrt(d)+sqrt(m)) instead of the
ground-state spec's relative 5*sigma*sigma_max floor (which the dominant p0~0.95 ground mode
inflates until it swallows the excited singular value ~ p1*sqrt(dm)). Visibility law: mode n is
recoverable iff p_n*sqrt(dm) clears the noise edge, so depth K buys visibility as ~sqrt(K).

References: FCI spectrum restricted to HF-reachable states (as SPEC_qksd_excited.md); the noisy
comparator is the solver's own solve_excited under matched per-element noise (SPEC_qksd_noise.md
machinery, applied favourably to QKSD -- no lambda factor on its H elements). Median over noise
seeds. PySCF/qiskit, no block2; `make gates` runs it in its own process.
"""
import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from odmd import build_odmd_problem, noise_edge, odmd_spectrum, sample_odmd_spectrum

SYSTEMS = {
    "h4": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
    "n2": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
}
_CACHE = {}


def _case(key):
    """(mh, ODMDProblem(n=48), reachable FCI energies in the centered frame)."""
    if key not in _CACHE:
        mh = build_molecular_hamiltonian(**SYSTEMS[key])
        prob = build_odmd_problem(mh, n=48)
        w_eig, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
        psi0 = np.asarray(mh.hf_state().data, dtype=complex)
        pops = np.abs(V.conj().T @ psi0) ** 2
        _CACHE[key] = (mh, prob, w_eig[pops > 1e-8].real - prob.mu)
    return _CACHE[key]


def _noisy_gap_errors(prob, gap_fci, shots, n, seeds=100, c=1.2):
    errs = []
    for sd in range(seeds):
        e = sample_odmd_spectrum(prob, shots, sd, n=n, c=c)
        errs.append(abs((e[1] - e[0]) - gap_fci) if len(e) > 1 else np.inf)
    return np.array(errs)


def test_G1_same_signal_recovers_the_spectrum():
    """Noiseless: E1 and the gap match reachable-FCI < 1e-5 Ha on H4 (K=24) and N2 (K=48)."""
    for key, K in (("h4", 24), ("n2", 48)):
        _, prob, ex = _case(key)
        E, _, _ = odmd_spectrum(prob.s[:K], prob.tau, cutoff=0.0)
        assert abs(E[1] - ex[1]) < 1e-5, (key, E[1] - ex[1])
        assert abs((E[1] - E[0]) - (ex[1] - ex[0])) < 1e-5, key


def test_G2_depth_is_the_excited_resource():
    """Noiseless boundary: on H4 the K=16 signal misses the gap by > 1 mHa while K=24 nails it --
    excited eigenphases need deeper K than the ground state."""
    _, prob, ex = _case("h4")
    gap_fci = ex[1] - ex[0]
    E16, _, _ = odmd_spectrum(prob.s[:16], prob.tau, cutoff=0.0)
    E24, _, _ = odmd_spectrum(prob.s[:24], prob.tau, cutoff=0.0)
    assert abs((E16[1] - E16[0]) - gap_fci) > 1e-3, (E16[1] - E16[0]) - gap_fci
    assert abs((E24[1] - E24[0]) - gap_fci) < 1e-5, (E24[1] - E24[0]) - gap_fci


def test_G3_noisy_gap_beats_qksd_at_matched_noise():
    """DEFINITION OF DONE: H4, K=48, 1e5 shots/element: median gap error < 10 mHa, resolved in
    >= 95% of seeds, and > 10x below noisy QKSD solve_excited at BOTH M=16 and M=24 (QKSD is
    noiselessly converged at M=24, so noise -- not depth -- is its limit)."""
    mh, prob, ex = _case("h4")
    gap_fci = float(ex[1] - ex[0])
    shots = 100_000
    errs = _noisy_gap_errors(prob, gap_fci, shots, n=48)
    m_odmd = float(np.median(errs))
    assert np.mean(np.isinf(errs)) <= 0.05, np.mean(np.isinf(errs))
    assert m_odmd < 10e-3, m_odmd
    sigma = np.sqrt(2.0 * (2.0 - 1.0 / prob.dim) / shots)
    for M in (16, 24):
        solver = QuantumKrylovSolver(mh, noise_sigma=sigma, seed=1)
        kerrs = []
        for _ in range(100):
            es = solver.solve_excited(M, n_states=2).energies
            kerrs.append(abs((es[1] - es[0]) - gap_fci) if len(es) > 1 else np.inf)
        m_kqd = float(np.median(kerrs))
        assert m_kqd / m_odmd > 10.0, (M, m_kqd, m_odmd)


def test_G4_threshold_mechanism_and_visibility_law():
    """(a) The ground-state spec's relative 5*sigma*sigma_max floor loses the excited mode in
    >= 90% of seeds at K=48; (b) at K=16 even the noise edge cannot see it (>= 50% unresolved:
    p1*sqrt(dm) below the edge) while K=48 resolves >= 95% -- the sqrt(K) visibility onset."""
    _, prob, ex = _case("h4")
    gap_fci = float(ex[1] - ex[0])
    shots = 100_000
    sigma = np.sqrt(2.0 * (2.0 - 1.0 / prob.dim) / shots)
    # (a) relative floor: rebuild the SPEC_odmd cutoff (5*sigma relative to sigma_max) by hand
    from scipy.linalg import hankel, svd
    sig_max = svd(hankel(prob.s[:24], prob.s[23:47]), compute_uv=False)[0]
    rel_unresolved = 0
    for sd in range(100):
        rng = np.random.default_rng(sd)
        g = rng.normal(0, sigma / np.sqrt(2), 48) + 1j * rng.normal(0, sigma / np.sqrt(2), 48)
        g[0] = 0.0
        E, _, _ = odmd_spectrum(prob.s + g, prob.tau, cutoff=5.0 * sigma * sig_max)
        rel_unresolved += int(len(E) < 2)
    assert rel_unresolved >= 90, rel_unresolved
    # (b) visibility onset with depth at the noise edge
    assert noise_edge(sigma, 8, 8) > noise_edge(sigma, 8, 8, c=0.0)  # sanity: positive edge
    unresolved_16 = float(np.mean(np.isinf(_noisy_gap_errors(prob, gap_fci, shots, n=16))))
    unresolved_48 = float(np.mean(np.isinf(_noisy_gap_errors(prob, gap_fci, shots, n=48))))
    assert unresolved_16 >= 0.50, unresolved_16
    assert unresolved_48 <= 0.05, unresolved_48
