"""
Acceptance gates for specs/SPEC_excited_state_certification.md (certified excited-state bounds).

The spec's v1.0 headline -- a *self-certified* excited-state solver, rigorous with no classical
oracle -- was FALSIFIED during implementation, and so were two of its four gates. The gates below
encode the falsifications rather than the claims (specs/README.md step 5):

  G1   oracle-mode containment of E_0 AND E_1 on LiH / H4 / N2 CAS(6,6), every M in 4..24.
  G1b  what is actually implemented: Lehmann's 2-vector pencil, which is provably tighter than the
       per-vector Temple inequality the spec's own sec 3.2 wrote down (and mislabelled Lehmann).
  G1c  the bounds are BOUNDS, not estimates: 300 seeded random Hermitian matrices, zero escapes --
       and escapes appear the moment the eigenvalue COUNT below the separator is off by one, which
       is the exact mechanism that kills G4.
  G2   THE FIRST KILL (replaces v1.0's "micro-Hartree at M=12 on STO-3G H2"): STO-3G H2's
       HF-reachable sector holds exactly TWO levels, so E_2 does not exist, no separator exists,
       and E_1 cannot be bounded from below at any M. The Krylov space saturates it at rank 2 from
       M=2 and the Ritz values are exact to 1e-15 -- there is nothing to converge. Convergence is
       re-gated on LiH/H4, where it is 3+ orders but NOT monotone.
  G3   v1.0's "< 5x overhead" survives on LiH (1.23-1.30x) and FAILS on H4 (5.42x) and N2
       (11.68x). Re-gated as the measured ceiling.
  G4   THE SECOND KILL: self-certified mode is NOT rigorous. beta = theta_2 - sigma_2 is an
       estimate; it escapes on real molecules (LiH M=4,6,8 by 2.0 mHa; N2 M=16 by 17.4 mHa).
  G4b  WHY it can never be fixed: a reachable level with tiny HF overlap is invisible to the
       subspace. A constructed witness hides a level of amplitude 1e-4 (population 1e-8, above
       this repo's 1e-10 reachability cut); the self-certified bound overshoots E_1 by 0.500 Ha
       while sigma_1 = 7e-5 says "converged". Also pins the sound condition and its conditioning
       limit.

All RNG seeded -> deterministic. PySCF/qiskit, no block2; `make gates` isolates it.
"""
import numpy as np
import pytest

from excited_bounds import (
    bracket_from_states,
    bracket_ladder,
    lehmann_excited_brackets,
    lehmann_lower_bounds,
    reachable_separator,
)
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver, ritz_pairs
from reachability import reachable_eigenpairs

DIMS = (4, 6, 8, 10, 12, 16, 20, 24)
SYSTEMS = {
    "LiH": dict(atom="Li 0 0 0; H 0 0 1.6", active_electrons=2, active_orbitals=5),
    "H4": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
    "N2": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
}
_CACHE = {}


def _case(name):
    """(E_0, E_1, E_2, oracle ladder, self ladder) for one system -- built once, reused by every
    gate. Reference energies are the reachable-sector eigenvalues (dense; validation oracle)."""
    if name not in _CACHE:
        mh = build_molecular_hamiltonian(**SYSTEMS[name])
        w, _ = reachable_eigenpairs(mh)
        off = mh.energy_offset
        beta = reachable_separator(mh)
        solver = QuantumKrylovSolver(mh)
        _CACHE[name] = (float(w[0]) + off, float(w[1]) + off, beta,
                        bracket_ladder(mh, DIMS, beta_oracle=beta, solver=solver),
                        bracket_ladder(mh, DIMS, beta_oracle=None, solver=solver))
    return _CACHE[name]


def _hidden_level_witness(m=4, amplitude=1e-4, dt=np.pi / 3.0):
    """A spectrum with a REACHABLE level the Krylov space cannot see, and its Ritz vectors.

    Levels {0, 0.5, 1, 2, 3}; the 0.5 level carries amplitude 1e-4, i.e. population 1e-8 -- ABOVE
    the 1e-10 HF-reachability cut this repo's certified arc uses (reachability.py), so it is part
    of the sector the bound is about. H is diagonal, so the real-time Krylov vectors are exact
    phase sums and no propagator error enters. Returns (H, ritz states, ritz values).
    """
    w = np.array([0.0, 0.5, 1.0, 2.0, 3.0])
    c = np.array([0.6, amplitude, 0.6, 0.4, 0.34])
    c = c / np.linalg.norm(c)
    H = np.diag(w).astype(complex)
    B = np.array([c * np.exp(-1j * k * dt * w) for k in range(m)])       # Krylov basis (m, 5)
    S = B.conj() @ B.T
    Hs = B.conj() @ (H @ B.T)
    vals, C, _ = ritz_pairs(0.5 * (Hs + Hs.conj().T), 0.5 * (S + S.conj().T))
    return H, C.T @ B, vals


@pytest.mark.parametrize("name", list(SYSTEMS))
def test_G1_simultaneous_containment_with_an_oracle_separator(name):
    """Both exact eigenvalues sit inside their brackets, at every gated depth, in oracle mode.

    Measured: 24 brackets over LiH / H4 / N2 CAS(6,6) at M = 4..24, ZERO escapes of either E_0 or
    E_1. The gate also requires a minimum number of FINITE lower bounds per system, so it cannot
    be satisfied by refusing everywhere (-inf is a valid but vacuous bound): LiH 8/8, H4 7/8,
    N2 5/8 -- the refusals are exactly the depths where theta_1 has not yet dropped below E_2.

    Containment is decided well above the arithmetic floor: the tightest margins over the whole
    ladder are E_0 - L_0 = 1.1e-11 Ha and E_1 - L_1 = 1.4e-6 Ha (both H4, M = 20), against a
    run-to-run scatter of 9e-15 and 9e-9 respectively (four repeat builds), so these are
    3+ orders of separation rather than luck.
    """
    e0, e1, _, ladder, _ = _case(name)
    finite = 0
    for br in ladder:
        assert br.lower_0 <= e0 <= br.upper_0, (name, br.m, br.lower_0 - e0, br.upper_0 - e0)
        assert br.lower_1 <= e1 <= br.upper_1, (name, br.m, br.lower_1 - e1, br.upper_1 - e1)
        finite += int(np.isfinite(br.lower_1))
    assert finite >= 5, (name, finite)


def test_G1b_lehmann_is_tighter_than_the_per_vector_temple_bound():
    """What is implemented is Lehmann's pencil, not the spec's sec 3.2 formula -- and it is better.

    sec 3.2 wrote E_i >= theta_i - sigma_i^2/(beta - theta_i), the single-vector Temple/Kato
    inequality, under the name of Lehmann's multi-dimensional theorem. Restricting the pencil to
    span{u_1} reproduces that formula algebraically, and Cauchy interlacing of the pencil makes the
    d-dimensional root the larger one. Measured on LiH at M = 4: Lehmann -0.964511 vs Temple
    -0.964578, i.e. +0.067 mHa of bound.

    The comparison tolerance is the pipeline's arithmetic floor, not slack: the Ritz vectors of a
    thresholded canonical orthogonalisation are orthonormal only to ~1e-9 (measured on H4), and the
    pencil divides that defect by beta - theta_1 ~ 0.3 Ha. The inequality therefore inverts by up
    to 4.5e-8 Ha on H4 at M = 10-12, and the number moves run to run because the default dt comes
    from ARPACK with a random start vector. Nothing in this module is a certificate below ~1e-8 Ha.
    """
    gains = []
    for name in SYSTEMS:
        for br in _case(name)[3]:
            if np.isfinite(br.lower_1) and np.isfinite(br.temple_lower_1):
                assert br.lower_1 >= br.temple_lower_1 - 1e-6, (name, br.m)
                gains.append(br.lower_1 - br.temple_lower_1)
    assert max(gains) > 5e-5, max(gains)          # LiH M=4: 6.7e-5 Ha


def test_G1c_the_pencil_is_a_bound_not_an_estimate():
    """300 seeded random Hermitian matrices: zero escapes -- and escapes the instant the count is off.

    Lehmann's hypothesis is that exactly nu eigenvalues lie below beta; the j-th largest root then
    bounds E_{nu-j}. The subspaces here are perturbed exact eigenvectors (eta = 0.02), so the
    bounds are TIGHT -- the closest comes within 1.6 mHa of its eigenvalue -- which is what gives
    the containment assertion teeth; a loose bound would satisfy it for free. Measured over 300
    8x8 matrices: worst escape exactly 0.0. With nu misstated by one -- the G4 failure mode, in
    miniature and at the matrix level -- 510 violations appear. A test that only checked the first
    half would pass on an estimate too.
    """
    rng = np.random.default_rng(0)
    worst_true, bad_shifted, tightest = 0.0, 0, -np.inf
    for _ in range(300):
        a = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
        H = 0.5 * (a + a.conj().T)
        w, V = np.linalg.eigh(H)
        nu = int(rng.integers(2, 8))                         # eigenvalues below the separator
        d = int(rng.integers(1, min(nu, 3) + 1))
        W = V[:, [nu - 1 - j for j in range(d)]] + 0.02 * (
            rng.normal(size=(8, d)) + 1j * rng.normal(size=(8, d)))
        q, _ = np.linalg.qr(W)
        _, y = np.linalg.eigh(q.conj().T @ H @ q)
        states = (q @ y).T                                   # rows are Ritz vectors
        beta = 0.5 * (w[nu - 1] + w[nu])                     # exactly nu eigenvalues below beta
        for j, b in enumerate(lehmann_lower_bounds(H, states, beta)):
            if nu - 1 - j >= 0:
                worst_true = max(worst_true, b - w[nu - 1 - j])
                tightest = max(tightest, b - w[nu - 1 - j])
            if nu - 2 - j >= 0:                              # pretend nu was one smaller
                bad_shifted += int(b - w[nu - 2 - j] > 1e-9)
    assert worst_true <= 1e-12, worst_true
    assert tightest > -2e-3, tightest                        # the bounds are tight, not vacuous
    assert bad_shifted > 100, bad_shifted


def test_G2_sto3g_h2_has_no_separator_at_all():
    """THE KILL: v1.0's G2 system cannot be certified at any M, and has nothing to converge.

    v1.0 required the excited bracket to close "monotonically to micro-Hartree at M = 12 on STO-3G
    H2". Measured: the HF-reachable sector of STO-3G H2 holds exactly 2 levels, so E_2 does not
    exist and NO separator theta_1 < beta <= E_2 can be chosen -- the lower bound on E_1 is -inf at
    every M by construction, not by lack of depth. The Krylov space saturates that sector at rank 2
    from M = 2 onward with theta_i - E_i = +-1e-15 and sigma_1 = 0, so M = 12 carries exactly the
    same information as M = 2 (at that roundoff level the sign of theta_1 - E_1 is not even fixed,
    which is why the Poincare upper bound is quoted with a tolerance here).
    """
    mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")
    w, _ = reachable_eigenpairs(mh)
    assert len(w) == 2, len(w)                               # no E_2 exists
    assert not np.isfinite(reachable_separator(mh))
    solver = QuantumKrylovSolver(mh)
    e1 = float(w[1]) + mh.energy_offset
    for m in (2, 4, 8, 12, 24):
        br = lehmann_excited_brackets(mh, m, beta_oracle=reachable_separator(mh), solver=solver)
        assert br.rank == 2, (m, br.rank)                    # saturated from M=2
        assert br.lower_1 == -np.inf                         # refusal, at every depth
        assert abs(br.upper_1 - e1) < 1e-13, (m, br.upper_1 - e1)


def test_G2b_the_bracket_closes_by_orders_but_not_monotonically():
    """Where convergence IS measurable, it is 3+ orders of magnitude -- and non-monotone.

    Replaces the "closes cleanly and monotonically" half of v1.0's G2. Measured U_1 - L_1 (mHa):
    LiH 88.26 -> 89.01 -> 90.06 -> 3.89 -> 3.98 -> 4.19 -> 0.0256 -> 0.0277 over M = 4..24;
    H4 981 -> 215.7 -> 113.3 -> 15.77 -> 6.98 -> 0.752 -> 0.0017 -> 0.0033. Both close by >= 3
    orders from M=6 to M=20 and both WIDEN again at least once -- LiH at five of the seven steps,
    H4 at M = 20 -> 24. The Krylov basis grows ill-conditioned faster than it grows informative,
    so a gate on monotonicity would fail.
    """
    for name, floor in (("LiH", 3e-2), ("H4", 3e-3)):
        widths = {br.m: br.width_1 * 1e3 for br in _case(name)[3]}
        assert widths[20] < floor, (name, widths[20])
        assert widths[6] / widths[20] > 1e3, (name, widths[6] / widths[20])
        increases = [m for a, m in zip(DIMS, DIMS[1:]) if widths[m] > widths[a]]
        assert increases, (name, widths)                     # NOT monotone: LiH [6,8,12,16,24], H4 [24]


def test_G3_certification_overhead_is_order_unity_but_not_under_5x():
    """v1.0's "< 5x the raw uncertified Ritz gap error" holds on LiH only. Re-gated as a ceiling.

    Overhead = (certified gap-bracket width) / |Ritz gap - exact gap| -- how much wider the
    certificate is than the error it is covering, which is the honest way to price certification.
    Measured over all finite brackets: LiH 1.23-1.30x, H4 2.30-5.42x, N2 CAS(6,6) 2.58-11.68x.
    So the spec's 5x is FALSE on two of three systems; what is true is that the overhead is O(1)
    at all -- a variance-based interval is usually orders looser, and this one is not, because
    bracket width and Ritz error both scale as sigma^2/gap.
    """
    measured = {}
    for name in SYSTEMS:
        e0, e1, _, ladder, _ = _case(name)
        ratios = [br.gap_width / abs((br.upper_1 - br.upper_0) - (e1 - e0))
                  for br in ladder if np.isfinite(br.gap_width)]
        measured[name] = (min(ratios), max(ratios))
    assert measured["LiH"][1] < 1.5, measured["LiH"]              # spec's claim survives here
    assert measured["H4"][1] > 5.0, measured["H4"]                # ... and dies here: 5.42x
    assert measured["N2"][1] > 5.0, measured["N2"]                # ... and here: 11.68x
    assert max(hi for _, hi in measured.values()) < 12.0, measured   # the honest ceiling
    assert min(lo for lo, _ in measured.values()) > 1.0, measured    # never cheaper than the error


def test_G4_self_certified_mode_is_not_rigorous():
    """THE KILL: v1.0's G4 ("must be mathematically rigorous ... for all M >= 8") is false.

    beta = theta_2 - sigma_2 is an ESTIMATE of E_2. Kato's interval guarantees only that SOME
    eigenvalue lies within sigma_2 of theta_2, never that it is E_2, so beta can sit above E_2;
    three eigenvalues are then below it, and the pencil's top root -- a true bound on E_2 -- is
    returned as a bound on E_1. Measured escapes: LiH M = 4, 6, 8 (E_1 under-bounded by 2.02, 1.97,
    1.89 mHa) and N2 CAS(6,6) M = 16 (by 17.4 mHa, i.e. at M well past the spec's M >= 8 claim).
    The sufficient condition is one-directional and that is the point: every escape has
    beta_self > E_2, but beta_self > E_2 does NOT always escape (N2 M = 4 survives by luck), so a
    self-mode run that passes is not evidence of anything.
    """
    escapes, lucky = [], 0
    for name in SYSTEMS:
        _, e1, beta_true, _, ladder = _case(name)
        for br in ladder:
            if np.isfinite(br.lower_1) and br.lower_1 > e1:
                escapes.append((name, br.m, (br.lower_1 - e1) * 1e3))
                assert br.beta > beta_true, (name, br.m)      # the mechanism, every time
            elif np.isfinite(br.lower_1) and br.beta > beta_true:
                lucky += 1                                    # unsound premise, survived anyway
    assert len(escapes) >= 4, escapes
    assert max(e[2] for e in escapes) > 15.0, escapes         # N2 M=16: 17.4 mHa
    assert lucky > 0, "unsound separators that happen to pass are what make self mode untestable"


def test_G4b_a_hidden_level_defeats_self_certification_by_construction():
    """WHY G4 can never be repaired: the data a self-certified beta is built from cannot see it.

    The witness spectrum is {0, 0.5, 1, 2, 3} with amplitude 1e-4 on the 0.5 level -- population
    1e-8, ABOVE the 1e-10 reachability cut the certified arc uses, so E_1 = 0.5 genuinely belongs
    to the sector being certified. At M = 4 the Krylov space returns Ritz values [0, 1, 2, 3]: the
    level is simply absent, and sigma_1 = 6.8e-5 reports a converged subspace. Self mode then
    picks beta = 1.9999 and certifies E_1 >= 1.000000 -- an overshoot of 0.500 Ha on a quantity
    whose true value is 0.5. Same failure class as SPEC_subspace_floor_resolvability's
    ~1e-4-amplitude level and SPEC_reachability_tolerance's threshold divergence.

    The sound condition, and its limit: with the true separator beta = E_2 = 1.0 the precondition
    theta_1 < beta survives only by 3e-9 (theta_1 IS the level the subspace mistook for E_1), the
    pencil is correspondingly ill-conditioned, and the bound collapses onto E_1 +- 1e-5 -- it stops
    lying, and stops saying anything. A separator must be strictly separated from theta_1 to carry
    digits, not merely valid.
    """
    H, states, vals = _hidden_level_witness(m=4)
    assert np.allclose(vals, [0.0, 1.0, 2.0, 3.0], atol=1e-6), vals   # 0.5 is invisible

    self_br = bracket_from_states(H, states, beta=None, m=4)
    assert self_br.sigma1 < 1e-4, self_br.sigma1                      # "converged"
    assert self_br.beta > 1.0, self_br.beta                           # beta overshoots E_2 = 1.0
    assert self_br.lower_1 - 0.5 > 0.4, self_br.lower_1               # 0.500 Ha escape

    oracle = bracket_from_states(H, states, beta=1.0, m=4)            # beta = E_2, the true one
    assert 0.0 < 1.0 - oracle.upper_1 < 1e-7, oracle.upper_1          # beta - theta_1 ~ 3e-9
    assert abs(oracle.lower_1 - 0.5) < 1e-5, oracle.lower_1           # collapsed to the boundary
    assert oracle.lower_0 <= 0.0 + 1e-9, oracle.lower_0               # E_0 = 0 still bounded
