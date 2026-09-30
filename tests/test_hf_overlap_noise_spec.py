"""
Acceptance gates for specs/SPEC_hf_overlap_noise.md (chem-mqu): the HF overlap certificate's
noise sibling. Unlike every other noisy certificate in this repo, a DETERMINANT guiding state's
lambda_u and r are exact/classical -- all shot noise enters through the single gap floor beta.

Prediction under test: coverage does not degrade smoothly with margin (delta_exact - r); it
bifurcates -- and the z-inflation knob that buys coverage everywhere else in the noise arc
(`certified_noise`, `gap_selfcheck_noise`) can instead COST coverage here, because the dominant
failure mode flips from "invalid" (overclaim, fixed by pushing beta down) at wide margin to
"vacuous" (fixed by nothing -- pushing beta down only makes it worse) at thin margin.
"""
import pytest

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

M = 8
SHOTS = (1e4, 1e5, 1e6)
Z_GRID = (0.0, 1.0, 2.0, 3.0)
TRIALS = 6000

# Four systems spanning the margin (delta_exact - r) from wide to vacuous, per the bead's own
# framing: H2 eq (wide, gamma=0.9936) -> H4 moderate -> H4 thin -> square H4 a=1.2 (vacuous, the
# self-mode M=8 direct certificate is already vacuous at ZERO noise).
CASES = {
    "H2_wide": dict(atom="H 0 0 0; H 0 0 0.74"),
    "H4_moderate": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
    "H4_thin": dict(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"),
    "sqH4_vacuous": dict(atom="H 0 0 0; H 1.2 0 0; H 1.2 1.2 0; H 0 1.2 0"),
}


@pytest.fixture(scope="module")
def grid():
    """The full {system x shots x z} coverage grid, >= 6000 trials/cell, seed=0."""
    from hf_overlap_noise import hf_overlap_noise_coverage

    out = {}
    for name, spec in CASES.items():
        mh = build_molecular_hamiltonian(**spec)
        solver = QuantumKrylovSolver(mh)
        for shots in SHOTS:
            for z in Z_GRID:
                out[(name, shots, z)] = hf_overlap_noise_coverage(
                    mh, M, shots, z=z, trials=TRIALS, seed=0, solver=solver,
                )
    return out


def test_G0_margin_ordering_matches_the_intended_span(grid):
    """Sanity: the four systems really do span wide -> vacuous, and sqH4 is vacuous already at
    zero noise (margin_exact <= 0), which is the point of including it."""
    margins = {name: grid[(name, 1e5, 0.0)]["margin"] for name in CASES}
    assert margins["H2_wide"] > margins["H4_moderate"] > margins["H4_thin"] > 0 > margins["sqH4_vacuous"], margins


def test_G1_padding_is_one_sided_and_monotone(grid):
    """Structural/killable: the pad only ever pushes beta DOWN, so for a fixed system+shots,
    frac_vacuous is non-decreasing in z and frac_invalid is non-increasing in z. A violation means
    the one-sided-pad implementation is wrong."""
    for name in CASES:
        for shots in SHOTS:
            vac = [grid[(name, shots, z)]["frac_vacuous"] for z in Z_GRID]
            inv = [grid[(name, shots, z)]["frac_invalid"] for z in Z_GRID]
            assert vac == sorted(vac), (name, shots, "frac_vacuous not non-decreasing in z", vac)
            assert inv == sorted(inv, reverse=True), (
                name, shots, "frac_invalid not non-increasing in z", inv
            )


def test_G2_wide_margin_is_near_saturated_by_z2(grid):
    """At wide margin (H2, gamma=0.9936) the certificate's only failure mode is INVALID
    (overclaim), which z-padding fixes directly: coverage >= 0.995 at z=2 for every shot count."""
    for shots in SHOTS:
        cov = grid[("H2_wide", shots, 2.0)]["coverage"]
        assert cov >= 0.995, (shots, cov)
        assert grid[("H2_wide", shots, 2.0)]["frac_invalid"] <= 0.005, shots


@pytest.mark.parametrize("shots", SHOTS)
def test_G3_kill_criterion_z_le_2_does_NOT_restore_coverage_everywhere(grid, shots):
    """THE KILL CRITERION (pre-registered in the bead): dies if z <= 2 restores coverage >= 0.9
    on every system. It does not -- H4_moderate and H4_thin stay well below 0.9 at every z <= 2
    and every shot count, and sqH4 is pinned at 0 (see G6). The claim survives."""
    for name in ("H4_moderate", "H4_thin"):
        for z in (0.0, 1.0, 2.0):
            cov = grid[(name, shots, z)]["coverage"]
            assert cov < 0.9, (name, shots, z, cov)


def test_G4_bifurcation_sensitivity_to_z_spikes_near_the_vacuous_boundary(grid):
    """THE FINDING: coverage's sensitivity to the SAME z-sweep (0 -> 3) is not smooth across the
    margin spectrum -- it is small (saturating) at wide margin, and spikes by 3-4x right before
    the certificate goes structurally vacuous. This is the observable signature of gamma_min's
    unbounded derivative in delta as delta -> r."""
    shots = 1e4
    swing = {
        name: grid[(name, shots, 3.0)]["coverage"] - grid[(name, shots, 0.0)]["coverage"]
        for name in CASES
    }
    assert abs(swing["H4_thin"]) > 3 * abs(swing["H2_wide"]), swing
    assert abs(swing["H4_thin"]) > 2 * abs(swing["H4_moderate"]), swing


def test_G5_inflation_reverses_sign_between_regimes(grid):
    """THE NOVEL PART (not a repeat of "inflation buys coverage"): at shots=1e4, z-inflation
    INCREASES coverage at wide margin (fixes overclaim) but DECREASES it at moderate/thin margin
    (worsens vacuousness) -- the opposite direction from every other certificate in the noise arc,
    where z uniformly buys coverage (certified_noise, gap_selfcheck_noise)."""
    shots = 1e4
    swing = {
        name: grid[(name, shots, 3.0)]["coverage"] - grid[(name, shots, 0.0)]["coverage"]
        for name in ("H2_wide", "H4_moderate", "H4_thin")
    }
    assert swing["H2_wide"] > 0, swing
    assert swing["H4_moderate"] < 0, swing
    assert swing["H4_thin"] < 0, swing


def test_G6_already_vacuous_system_cannot_be_rescued_by_inflation(grid):
    """Honesty diagnostic: a system that is vacuous at ZERO noise (sqH4 a=1.2, self-mode M=8)
    stays vacuous at every shots/z tested -- the one-sided pad can only push beta further down, so
    it can never manufacture a certificate where the noiseless one does not exist."""
    for shots in SHOTS:
        for z in Z_GRID:
            r = grid[("sqH4_vacuous", shots, z)]
            assert r["frac_vacuous"] > 0.999, (shots, z, r)
            assert r["coverage"] < 0.001, (shots, z, r)


def test_G7_shots_always_help_unlike_z(grid):
    """Contrast worth recording: MORE SHOTS always increase coverage (shrinks the noise around the
    true, non-negative-at-zero-noise margin), for every system with a positive exact margin --
    unlike z, which can hurt. Checked at z=0 across the shot ladder."""
    for name in ("H2_wide", "H4_moderate", "H4_thin"):
        covs = [grid[(name, s, 0.0)]["coverage"] for s in SHOTS]
        assert covs == sorted(covs), (name, "coverage should be non-decreasing in shots", covs)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
