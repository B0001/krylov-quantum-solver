"""
Spec gates for SPEC_hf_overlap_subspace: certified HF SUBSPACE overlap on molecules (SPEC-21b).

Feeds the repo's premise-gated self-mode Krylov E_d floor (theta_d - sigma_d) into the SPEC-21b
block Davis-Kahan machinery and certifies gamma_min <= ||P_S u|| for the lowest-d reachable
eigenspace of square H4 -- non-vacuously where the d=1 single-vector certificate is vacuous, and
with no oracle.

  G1 validity (killable, zero-tol): gamma_min <= exact reachable ||P_S u||, self M in {6,8} + oracle.
  G2 THE FINDING: across the square-H4 sweep, d=1 certificate VACUOUS while d=2 non-vacuous (>=0.45).
  G3 no oracle needed: self-mode gamma matches oracle-mode within 0.05 at M >= 8.
  G4 premise boundary: self mode at M < 6 raises; oracle mode does not.
  G5 usefulness trend: d=2 self-mode floor non-decreasing as the square tightens (a: 1.4->1.2->1.0).
"""

import pytest

from hf_overlap_certificate import certify_hf_overlap
from hf_overlap_subspace import (
    _reachable_e_d_total,
    certify_hf_subspace_overlap,
    exact_hf_subspace_overlap,
)
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

_SIDES = (1.4, 1.2, 1.0)  # square H4 side lengths (Angstrom), rectangle -> tight square


def _square_h4(a: float):
    return build_molecular_hamiltonian(atom=f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0")


@pytest.fixture(scope="module")
def bundles():
    out = {}
    for a in _SIDES:
        mh = _square_h4(a)
        out[a] = {
            "mh": mh,
            "solver": QuantumKrylovSolver(mh),
            "exact_d2": exact_hf_subspace_overlap(mh, 2),
            "e_d2": _reachable_e_d_total(mh, 2),
        }
    return out


@pytest.mark.parametrize("a", _SIDES)
@pytest.mark.parametrize("m", [6, 8])
def test_G1_validity_self_and_oracle(bundles, a, m):
    """G1: certified d=2 floor never exceeds the exact reachable subspace overlap. Zero tol."""
    b = bundles[a]
    c_self = certify_hf_subspace_overlap(b["mh"], 2, m=m, solver=b["solver"])
    c_orac = certify_hf_subspace_overlap(b["mh"], 2, m=m, e_d=b["e_d2"], solver=b["solver"])
    for c, mode in ((c_self, "self"), (c_orac, "oracle")):
        floor = 0.0 if c.vacuous else c.gamma_min
        assert floor <= b["exact_d2"] + 1e-12, (
            f"G1 VIOLATION square-H4 a={a} M={m} {mode}: "
            f"gamma_min={floor} > exact ||P_S u||={b['exact_d2']}"
        )
        assert c.cluster_size == 2


@pytest.mark.parametrize("a", _SIDES)
def test_G2_d1_vacuous_d2_useful(bundles, a):
    """G2 (the finding): the single-vector certificate is vacuous where the d=2 block is useful."""
    b = bundles[a]
    c1 = certify_hf_overlap(b["mh"], m=8, solver=b["solver"])
    c2 = certify_hf_subspace_overlap(b["mh"], 2, m=8, solver=b["solver"])
    assert c1.vacuous, f"G2: expected d=1 VACUOUS on square-H4 a={a}, got gamma={c1.gamma_min}"
    assert not c2.vacuous, f"G2: expected d=2 non-vacuous on square-H4 a={a}"
    assert c2.gamma_min >= 0.45, f"G2: d=2 floor {c2.gamma_min} unexpectedly weak at a={a}"


@pytest.mark.parametrize("a", _SIDES)
def test_G3_self_matches_oracle(bundles, a):
    """G3: at M >= 8 the self-mode floor matches oracle within 0.05 -- certified with no oracle."""
    b = bundles[a]
    g_self = certify_hf_subspace_overlap(b["mh"], 2, m=8, solver=b["solver"]).gamma_min
    g_orac = certify_hf_subspace_overlap(b["mh"], 2, m=8, e_d=b["e_d2"], solver=b["solver"]).gamma_min
    assert g_self <= g_orac + 1e-12, f"G3: self {g_self} exceeds oracle {g_orac} at a={a}"
    assert abs(g_self - g_orac) <= 0.05, f"G3: self {g_self} vs oracle {g_orac} gap > 0.05 at a={a}"


def test_G4_self_mode_below_M6_raises(bundles):
    """G4: the self-mode premise boundary (M >= 6) is inherited as a loud raise."""
    b = bundles[1.2]
    with pytest.raises(ValueError, match="m >= 6"):
        certify_hf_subspace_overlap(b["mh"], 2, m=4, solver=b["solver"])
    # oracle mode has no premise -- must not raise below M=6
    c = certify_hf_subspace_overlap(b["mh"], 2, m=4, e_d=b["e_d2"], solver=b["solver"])
    assert c.gamma_min <= b["exact_d2"] + 1e-12


def test_G5_usefulness_trend(bundles):
    """G5: the d=2 self-mode floor is non-decreasing as the square tightens (1.4 -> 1.2 -> 1.0)."""
    floors = [
        certify_hf_subspace_overlap(bundles[a]["mh"], 2, m=8, solver=bundles[a]["solver"]).gamma_min
        for a in (1.4, 1.2, 1.0)
    ]
    assert floors[0] <= floors[1] + 1e-9 <= floors[2] + 1e-9, (
        f"G5: d=2 self-mode floors not non-decreasing as square tightens: {floors}"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
