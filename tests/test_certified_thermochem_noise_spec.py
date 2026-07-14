"""
Acceptance gates G1-G4 for specs/SPEC_certified_thermochem_noise.md (the certified relative-energy
bracket under shot noise -- composition beats either endpoint alone).

`certified_thermochem_noise.thermochem_noise_coverage`/`minimal_z_for_coverage` do not exist yet:
this file is RED until they are implemented.
"""
import numpy as np
import pytest

from certified_noise import shot_noise_coverage
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

M = 16
SHOTS = (1e4, 1e5, 1e6)


def chain(r):
    return f"H 0 0 0; H 0 0 {r}; H 0 0 {2 * r}; H 0 0 {3 * r}"


@pytest.fixture(scope="module")
def geometries():
    return build_molecular_hamiltonian(atom=chain(0.9)), build_molecular_hamiltonian(atom=chain(2.3))


def test_G1_composition_inherits_collapse_and_is_shot_independent(geometries):
    """Raw coverage of the composed relative-energy bracket is broken (< 0.85) at every shot count,
    and (mirroring `certified_noise`) roughly N-independent -- shots do not buy this back."""
    from certified_thermochem_noise import thermochem_noise_coverage

    mh_a, mh_b = geometries
    covs = [thermochem_noise_coverage(mh_a, mh_b, M, s, z=0.0)["cov_raw"] for s in SHOTS]
    assert all(c < 0.85 for c in covs), covs
    assert max(covs) - min(covs) < 0.05, covs


def test_G2_existing_single_bracket_rule_still_restores_coverage(geometries):
    """Reusing the certified_noise z=2 inflation rule unchanged (applied per-endpoint before
    composing) still restores >= 0.9 coverage on the composed bracket -- a sanity/regression check
    that composition does not break the existing rule."""
    from certified_thermochem_noise import thermochem_noise_coverage

    mh_a, mh_b = geometries
    covs = [thermochem_noise_coverage(mh_a, mh_b, M, s, z=2.0)["cov_inflated"] for s in SHOTS]
    assert all(c >= 0.9 for c in covs), covs


def test_G3_composition_needs_less_inflation_than_single_bracket(geometries):
    """THE FINDING / definition of done: the minimal z restoring 90% coverage on the COMPOSED
    bracket is strictly less than the single-bracket z=2 rule -- composition partially
    self-corrects the coin-flip collapse. Also must be > 0: G1 established raw (z=0) coverage is
    already broken, so some inflation is still required."""
    from certified_thermochem_noise import minimal_z_for_coverage

    mh_a, mh_b = geometries
    z_star = minimal_z_for_coverage(mh_a, mh_b, M, 1e5, target=0.9, resolution=0.05)
    assert 0.0 < z_star < 2.0, z_star


def test_G4_composed_coverage_beats_either_single_endpoint(geometries):
    """Mechanism check: at z=0, the composed relative-energy bracket's raw coverage exceeds EITHER
    single-geometry bracket's raw coverage (certified_noise, same M/shots/seed) by a wide margin,
    not a rounding-level one -- the two independent coin-flips do not simply compound."""
    from certified_thermochem_noise import thermochem_noise_coverage

    mh_a, mh_b = geometries
    shots = 1e5
    composed = thermochem_noise_coverage(mh_a, mh_b, M, shots, z=0.0)["cov_raw"]
    single_a = shot_noise_coverage(mh_a, M, shots, z=0.0)["cov_raw"]
    single_b = shot_noise_coverage(mh_b, M, shots, z=0.0)["cov_raw"]
    assert composed >= single_a + 0.2, (composed, single_a)
    assert composed >= single_b + 0.2, (composed, single_b)


def test_G_reference_matches_certified_thermochem_exact_delta(geometries):
    """Sanity: the exact relative energy used as ground truth here is the same quantity
    `certified_thermochem` validates its noiseless bracket against (dense diagonalization,
    per-geometry energy_offset included -- offsets differ across geometries and do not cancel)."""
    mh_a, mh_b = geometries
    delta_exact = mh_b.ground_state_energy() - mh_a.ground_state_energy()
    assert np.isclose(delta_exact, 0.3022795693589475, atol=1e-6)
