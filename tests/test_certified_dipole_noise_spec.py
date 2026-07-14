"""
Acceptance gates G1-G4 for specs/SPEC_certified_dipole_noise.md (certified dipole under shot noise
-- inflation cannibalizes the gap margin it needs to stay finite).

`certified_dipole_noise.dipole_noise_coverage` does not exist yet: this file is RED until it is
implemented.
"""
import pytest

from hybrid_quantum_solver.molecular_hamiltonian import build_dipole_operators, build_molecular_hamiltonian

M = 16
SHOTS = (1e4, 1e5, 1e6)

CASES = {
    "HeH+": dict(atom="He 0 0 0; H 0 0 0.772", charge=1),
    "LiH": dict(atom="Li 0 0 0; H 0 0 1.6"),
}


@pytest.fixture(scope="module")
def systems():
    out = {}
    for name, spec in CASES.items():
        mh = build_molecular_hamiltonian(**spec)
        az = build_dipole_operators(**spec)[2]
        out[name] = (mh, az)
    return out


def test_G1_raw_coverage_broken_on_both_systems(systems):
    """Raw (z=0) coverage of the certified dipole bracket under shot noise is broken on both
    systems -- healthy-margin HeH+ less catastrophically than fragile-margin LiH."""
    from certified_dipole_noise import dipole_noise_coverage

    mh, az = systems["HeH+"]
    covs = [dipole_noise_coverage(mh, az, M, s, z=0.0)["coverage"] for s in SHOTS]
    assert all(c < 0.6 for c in covs), covs

    mh, az = systems["LiH"]
    covs = [dipole_noise_coverage(mh, az, M, s, z=0.0)["coverage"] for s in SHOTS]
    assert all(c < 0.05 for c in covs), covs


def test_G2_moderate_inflation_restores_coverage_on_healthy_margin_system(systems):
    """On HeH+ (healthy Delta_lo margin), z=1.0 restores >= 0.9 coverage at shots=1e5/1e6 --
    inflation still works here, at moderate z."""
    from certified_dipole_noise import dipole_noise_coverage

    mh, az = systems["HeH+"]
    for shots in (1e5, 1e6):
        cov = dipole_noise_coverage(mh, az, M, shots, z=1.0)["coverage"]
        assert cov >= 0.9, (shots, cov)


def test_G3_inflation_ceiling_at_tight_shot_budget(systems):
    """THE FINDING / definition of done: at shots=1e4, HeH+'s finite-bracket rate is strictly
    LOWER at z=3 than at z=1 -- more inflation shrinks the fraction of trials where a bracket can
    even be constructed, the opposite of certified_thermochem_noise / gap_selfcheck_noise (both
    monotonically non-decreasing in z)."""
    from certified_dipole_noise import dipole_noise_coverage

    mh, az = systems["HeH+"]
    r1 = dipole_noise_coverage(mh, az, M, 1e4, z=1.0)
    r3 = dipole_noise_coverage(mh, az, M, 1e4, z=3.0)
    assert r3["finite_frac"] < r1["finite_frac"], (r1, r3)
    assert r3["coverage"] < r1["coverage"], (r1, r3)


def test_G4_lih_margin_too_thin_for_inflation_to_rescue(systems):
    """Boundary, recorded not fixed: across every tested z and shot count, LiH's finite-bracket
    rate stays low -- inflation cannot rescue a system whose Delta_lo margin is already this thin
    at M=16."""
    from certified_dipole_noise import dipole_noise_coverage

    mh, az = systems["LiH"]
    fracs = [
        dipole_noise_coverage(mh, az, M, shots, z=z)["finite_frac"]
        for shots in SHOTS for z in (0.0, 1.0, 2.0, 3.0)
    ]
    assert max(fracs) < 0.3, fracs
