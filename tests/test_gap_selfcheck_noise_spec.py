"""
Acceptance gates G1-G4 for specs/SPEC_gap_selfcheck_noise.md (gap self-check under shot noise --
intersection concentrates noise where composition, certified_thermochem_noise, diluted it).

`gap_selfcheck_noise.self_check_noise_coverage`/`minimal_z_for_selfcheck_coverage` do not exist
yet: this file is RED until they are implemented.
"""
import pytest

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

DIMS = (6, 8, 12, 16, 20, 24)
SHOTS = (1e4, 1e5, 1e6)

CASES = {
    "H4": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
    "LiH": dict(atom="Li 0 0 0; H 0 0 1.6", active_electrons=2, active_orbitals=5),
}


@pytest.fixture(scope="module")
def mhs():
    return {name: build_molecular_hamiltonian(**spec) for name, spec in CASES.items()}


@pytest.mark.parametrize("name", CASES)
def test_G1_selfcheck_inherits_collapse_and_is_shot_independent(mhs, name):
    """Raw (z=0) coverage of the self-checked interval is broken (< 0.3) at every shot count, and
    (mirroring certified_noise) roughly N-independent."""
    from gap_selfcheck_noise import self_check_noise_coverage

    mh = mhs[name]
    covs = [self_check_noise_coverage(mh, DIMS, s, z=0.0)["coverage"] for s in SHOTS]
    assert all(c < 0.3 for c in covs), covs
    assert max(covs) - min(covs) < 0.1, covs


@pytest.mark.parametrize("name", CASES)
def test_G2_single_bracket_z2_rule_does_not_fix_it(mhs, name):
    """THE CONTRAST: z=2 -- the rule that restores certified_noise's single bracket to ~0.98 and
    is already more than enough for certified_thermochem_noise's composed bracket -- leaves the
    self-checked interval's coverage below 0.85 at shots=1e5."""
    from gap_selfcheck_noise import self_check_noise_coverage

    mh = mhs[name]
    cov = self_check_noise_coverage(mh, DIMS, 1e5, z=2.0)["coverage"]
    assert cov < 0.85, cov


@pytest.mark.parametrize("name", CASES)
def test_G3_selfcheck_needs_more_inflation_than_single_bracket(mhs, name):
    """THE FINDING / definition of done: the minimal z restoring 90% coverage on the self-checked
    interval is GREATER than 2.0 (the single-bracket rule) at shots=1e5 -- the opposite direction
    from certified_thermochem_noise's composed (difference) bracket, which needed less."""
    from gap_selfcheck_noise import minimal_z_for_selfcheck_coverage

    mh = mhs[name]
    z_star = minimal_z_for_selfcheck_coverage(mh, DIMS, 1e5, target=0.9, resolution=0.25)
    assert z_star > 2.0, z_star


@pytest.mark.parametrize("name", CASES)
def test_G4_padding_fixes_inconclusive_before_it_fixes_correct(mhs, name):
    """Honesty diagnostic: the self-checked interval is empty (no bracket corroborates) in a
    non-trivial fraction of raw (z=0) trials, but that fraction nearly vanishes once z=2.0 padding
    is applied -- padding restores a conclusive answer well before it restores a correct one."""
    from gap_selfcheck_noise import self_check_noise_coverage

    mh = mhs[name]
    raw = self_check_noise_coverage(mh, DIMS, 1e5, z=0.0)
    padded = self_check_noise_coverage(mh, DIMS, 1e5, z=2.0)
    assert raw["frac_empty"] > 0.05, raw
    assert padded["frac_empty"] < 0.02, padded
