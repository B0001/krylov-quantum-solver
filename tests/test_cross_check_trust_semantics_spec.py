"""
Acceptance gates G1-G4 for specs/SPEC_cross_check_trust_semantics.md (cross_check's own trust
semantics -- what "reference" means when CASCI can't run).

Deliberately no new library code: every gate is driven entirely through `cross_check.py`'s existing
public cost-cap keyword arguments (`fci_max_dim`, `krylov_max_dim`, `qubit_dense_max_orb`), never
touching internals -- a genuine external characterization of the harness's own fallback logic.
"""
import pytest
from pyscf import ao2mo, gto, mcscf, scf

from cross_check import cross_check

SYSTEMS = {
    "H2": ("H 0 0 0; H 0 0 0.74", 2, 2),
    "H4": ("H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3", 4, 4),
}


def _reference(label):
    atom, norb, ne = SYSTEMS[label]
    mol = gto.M(atom=atom, basis="sto-3g")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    cas = mcscf.CASCI(mf, norb, ne)
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), norb)
    return h1, eri, float(e_core), (ne // 2, ne // 2), norb


@pytest.mark.parametrize("label", SYSTEMS)
def test_G1_baseline_agreement_pins_the_informal_assertion(label):
    """With all four solvers reachable, they agree within the harness's own 5 mHa tolerance --
    pins the informal __main__ `assert out["agree"]` into a real gate."""
    h1, eri, e_core, nelec, norb = _reference(label)
    out = cross_check(h1, eri, e_core, nelec, norb, tol_mHa=5.0)
    assert out["skipped"] == [], out["skipped"]
    assert out["agree"], out
    assert out["max_dev_mHa"] <= 5.0, out["max_dev_mHa"]


def test_G2_reference_fallback_is_exact_insertion_order_not_a_ranking():
    """THE FINDING / definition of done: with CASCI forced unreachable, the "reference" becomes
    EXACTLY the Krylov value (not a recomputed best estimate) -- the fallback is a literal
    first-available pick in source insertion order, not a documented trust ranking."""
    h1, eri, e_core, nelec, norb = _reference("H4")
    out = cross_check(h1, eri, e_core, nelec, norb, tol_mHa=5.0, fci_max_dim=0)
    assert "CASCI" in out["skipped"], out["skipped"]
    assert "Krylov" in out["results"], out["results"]
    assert out["reference"] == out["results"]["Krylov"][0], out


def test_G3_fallback_priority_never_promotes_sqd_ahead_of_adapt():
    """With CASCI and Krylov both forced unreachable, the reference becomes ADAPT, not SQD -- the
    insertion order in the source (CASCI, Krylov, ADAPT, SQD) happens to keep the
    configuration-sampling method last in line, even though nothing in the code documents this as
    a trust ranking."""
    h1, eri, e_core, nelec, norb = _reference("H4")
    out = cross_check(h1, eri, e_core, nelec, norb, tol_mHa=5.0, fci_max_dim=0, krylov_max_dim=0)
    assert set(out["skipped"]) == {"CASCI", "Krylov"}, out["skipped"]
    assert "ADAPT" in out["results"] and "SQD" in out["results"], out["results"]
    assert out["reference"] == out["results"]["ADAPT"][0], out


def test_G4_sqd_cannot_be_suppressed_through_the_public_cost_caps():
    """Boundary, recorded not smoothed over: with ALL THREE other caps forced to zero, CASCI,
    Krylov, and ADAPT are all skipped, but SQD still ran unconditionally and became the sole
    reference -- there is no way to drive cross_check into a fully-empty "no reference" state
    through its public cost caps alone; SQD is not a symmetric fourth knob."""
    h1, eri, e_core, nelec, norb = _reference("H4")
    out = cross_check(h1, eri, e_core, nelec, norb, tol_mHa=5.0,
                      fci_max_dim=0, krylov_max_dim=0, qubit_dense_max_orb=0)
    assert set(out["skipped"]) == {"CASCI", "Krylov", "ADAPT"}, out["skipped"]
    assert "SQD" in out["results"], out["results"]
    assert out["reference"] == out["results"]["SQD"][0], out
    assert out["reference"] is not None
