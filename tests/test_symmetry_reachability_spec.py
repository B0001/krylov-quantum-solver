"""Gates for specs/SPEC_symmetry_reachability.md.

specs/SPEC_reachability_tolerance.md proved no fixed amplitude threshold separates a physical
HF overlap from an SCF convergence residue on a symmetry-forbidden state. The proposed replacement
is a spatial-irrep filter, whose obvious objection is that RHF breaks symmetry on exactly the
strongly-correlated systems that matter -- leaving no irrep to match.

THE FINDING: that objection is self-consistent, not fatal. The filter is available EXACTLY when the
artifact is present (9/9 on the square-H4 family). Where RHF breaks symmetry the reference genuinely
overlaps the low state (~0.45, not ~1e-10), so there is nothing to remove.

chem-8ls (2026-11-05): the "artifact" half was the free SCF stopping on the D2h saddle point, at
sides that depended on the platform's rounding. The builder now follows that instability to the
stable broken-symmetry minimum (specs/SPEC_scf_stability_rule.md). So on this family the artifact is
gone, and the filter unavailable, at every side. G1 still holds, as one regime; G2 and G3 assert the
branch the rule selects (SPEC section 10).
"""
import numpy as np
import pytest

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from reachability import (
    SYMMETRY_BREAK_TOL,
    TIGHT_SCF_CONV_TOL,
    hf_orbital_irreps,
    hf_population_spectrum,
    scf_symmetry_status,
    symmetry_filter_available,
)

# atom ORDER is load-bearing -- it fixes which determinant is HF (SPEC_reachability_tolerance R2).
def _square_h4(a):
    return f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0"


SWEEP = (1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40)
ARTIFACT_P0 = 1e-6          # below this, p0 is residue on a forbidden level, not an overlap


# --- G1: the correlation (DEFINITION OF DONE) -----------------------------------------------------

def test_G1_filter_is_available_exactly_when_the_artifact_is_present():
    """The whole finding, in one gate. Killed by a single mismatched geometry."""
    mismatches = []
    for a in SWEEP:
        geom = _square_h4(a)
        p0 = float(hf_population_spectrum(build_molecular_hamiltonian(atom=geom))[0])
        artifact = p0 < ARTIFACT_P0
        if artifact != symmetry_filter_available(geom):
            mismatches.append((a, p0, artifact))
    assert not mismatches, mismatches


# --- G2: one regime -- the stability rule leaves no artifact geometry (chem-8ls) -------------------

def test_G2_every_side_carries_a_physical_overlap():
    """Was: artifact geometries carry p0 < 1e-7, symmetry-broken ones p0 > 0.4, nothing between.
    The artifact geometries were the SCF stopping on the saddle (macOS: 1.05/1.10/1.35 in this
    sweep; Linux: 1.15/1.20/1.25/1.30). With the stability rule every side is on the minimum, and p0
    falls smoothly from 0.476 to 0.440 (SPEC section 10)."""
    p0 = [float(hf_population_spectrum(build_molecular_hamiltonian(atom=_square_h4(a)))[0])
          for a in SWEEP]
    assert min(p0) > 0.4 and max(p0) < 0.5, p0
    assert all(x > y for x, y in zip(p0, p0[1:])), p0              # monotone in a: no branch jump


# --- G3: the mechanism -- unavailability tracks a genuinely lower broken-symmetry solution --------

@pytest.mark.parametrize("a", SWEEP)
def test_G3_unavailable_means_RHF_found_a_lower_solution(a):
    """The filter refuses precisely because the reference is not a symmetry eigenfunction -- and
    that reference is variationally BETTER, so refusing is correct, not a shortcoming. Gated at 1.20
    and 1.40 only until chem-8ls; with the stability rule it holds at every side."""
    geom = _square_h4(a)
    broken, d_e = scf_symmetry_status(geom)
    assert broken, (a, d_e)
    assert d_e < -0.05, (a, d_e)                    # ~0.08 Ha lower, not a convergence wobble
    assert hf_orbital_irreps(geom) is None
    assert not symmetry_filter_available(geom)



# --- G4: the filter is not vacuous -- it must accept ordinary systems ------------------------------

@pytest.mark.parametrize("geom", [
    "H 0 0 0; H 0 0 0.74",                                  # H2 equilibrium
    "H 0 0 0; H 0 0 2.0",                                   # H2 stretched
    "H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0",             # linear H4
    "Li 0 0 0; H 0 0 1.6",                                  # LiH
])
def test_G4_filter_available_on_the_repo_s_ordinary_gated_systems(geom):
    """Killed if the filter refuses the systems the specs actually gate on -- it would then be
    unusable in practice regardless of how well it behaves on square H4."""
    assert symmetry_filter_available(geom), geom
    assert not scf_symmetry_status(geom)[0], geom


# --- G5: those ordinary systems carry no artifact in the first place -------------------------------

@pytest.mark.parametrize("geom", [
    "H 0 0 0; H 0 0 0.74",
    "H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0",
])
def test_G5_ordinary_systems_have_no_near_threshold_contamination(geom):
    """Bounds the scope of the whole reachability problem: on linear chains and H2 every population
    below the physical ground state is at the machine-zero floor, so no threshold choice matters.
    Square H4 is where this lives.
    """
    mh = build_molecular_hamiltonian(atom=geom)
    p = hf_population_spectrum(mh)
    physical = np.where(p > 1e-6)[0]
    below = p[: physical[0]] if physical[0] > 0 else np.array([0.0])
    assert below.max() < 1e-20, (geom, below.max())


# --- G6: the tight-SCF constant is the one that collapses the residue ------------------------------

def test_G6_tight_conv_tol_constant_is_pinned():
    assert TIGHT_SCF_CONV_TOL == 1e-13
    assert SYMMETRY_BREAK_TOL == 1e-6
