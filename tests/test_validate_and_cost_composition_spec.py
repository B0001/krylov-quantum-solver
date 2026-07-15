"""
Acceptance gates G1-G4 for specs/SPEC_validate_and_cost_composition.md (validate_and_cost -- one
threshold quietly governs two independent regime boundaries).

Deliberately no new library code: every gate is driven through `validate_and_cost.py`'s existing
public `validate_and_cost`/`print_report`, never touching internals.
"""
import pytest
from pyscf import ao2mo, gto, mcscf, scf

from validate_and_cost import _HAVE_FT, print_report, validate_and_cost


def _reference(atom, norb, ne, basis="sto-3g"):
    mol = gto.M(atom=atom, basis=basis)
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    cas = mcscf.CASCI(mf, norb, ne)
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), norb)
    return h1, eri, float(e_core), (ne // 2, ne // 2)


@pytest.fixture(scope="module")
def small_h4_report():
    h1, eri, e_core, nelec = _reference("H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3", 4, 4)
    return validate_and_cost(h1, eri, e_core, nelec, 4, target_mHa=1.0)


@pytest.fixture(scope="module")
def large_n2_report():
    # 8 orbitals, 10 active electrons (na=nb=5): exceeds the default qubit_dense_max_orb=7, but
    # CI dim (3136) stays under cross_check's own krylov_max_dim/fci_max_dim caps, so CASCI/
    # Krylov/SQD still run while only taper and ADAPT (both qubit_dense_max_orb-gated) drop out.
    h1, eri, e_core, nelec = _reference("N 0 0 0; N 0 0 1.0977", 8, 10)
    return h1, eri, e_core, nelec, validate_and_cost(h1, eri, e_core, nelec, 8, target_mHa=1.0)


def test_G1_small_cas_full_pipeline_composes(small_h4_report):
    """On H4 CAS(4,4): taper is not skipped and reduces qubit count; cross-check's four methods
    agree; FT cost gracefully reports None in the standard chem env (no openfermion)."""
    rep = small_h4_report
    assert _HAVE_FT is False, "expected no openfermion in the standard chem env"
    assert "skipped" not in rep["taper"], rep["taper"]
    assert rep["taper"]["n_qubits_tapered"] < rep["taper"]["n_qubits_original"], rep["taper"]
    assert rep["cross_check"]["agree"] is True, rep["cross_check"]
    assert rep["ft_cost"] is None, rep["ft_cost"]


def test_G2_one_threshold_governs_two_regime_boundaries(small_h4_report):
    """THE FINDING / definition of done: on the SAME small system (H4 CAS(4,4), norb=4), taper
    AND ADAPT (inside cross_check) both flip off together when qubit_dense_max_orb is lowered
    below norb, and both are on together at the default -- one forwarded threshold governs two
    conceptually independent regime boundaries.

    Deliberately does NOT test this by RAISING the threshold on the large (norb=8) N2 system: doing
    so would let taper actually attempt `pauli_decompose` at 16 qubits, which is exponential in
    qubit count (SPEC_taper_spectrum.md's own R1 measured ~1000s at 8 qubits; 16 is many orders of
    magnitude worse) -- this was tried while probing this spec and OOM-killed the process. Proving
    the SAME coupling by lowering the threshold on an already-small, already-tractable system is
    the safe direction."""
    h1, eri, e_core, nelec = _reference("H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3", 4, 4)

    rep_default = small_h4_report  # qubit_dense_max_orb=7 (default): both taper and ADAPT run
    assert "skipped" not in rep_default["taper"], rep_default["taper"]
    assert "ADAPT" not in rep_default["cross_check"]["skipped"], rep_default["cross_check"]["skipped"]

    rep_lowered = validate_and_cost(h1, eri, e_core, nelec, 4, target_mHa=1.0,
                                    qubit_dense_max_orb=3)  # below norb=4: both should skip
    assert "skipped" in rep_lowered["taper"], rep_lowered["taper"]
    assert "ADAPT" in rep_lowered["cross_check"]["skipped"], rep_lowered["cross_check"]["skipped"]


def test_G3_graceful_degradation_surviving_methods_still_agree(large_n2_report):
    """Despite taper and ADAPT both being out of regime, the surviving methods (CASCI, Krylov,
    SQD) still reach a valid, non-vacuous agreement verdict -- the pipeline degrades gracefully,
    not catastrophically."""
    _h1, _eri, _e_core, _nelec, rep = large_n2_report
    xc = rep["cross_check"]
    assert xc["reference"] is not None
    assert xc["agree"] is True, xc
    assert {"CASCI", "Krylov", "SQD"} <= set(xc["results"]), xc["results"]


def test_G4_print_report_does_not_crash_in_either_regime(small_h4_report, large_n2_report):
    """print_report is exercised directly (not just the structural dict checks) on both the
    full-pipeline (H4) and taper/ADAPT-skipped (N2) reports, with no exception."""
    _h1, _eri, _e_core, _nelec, large_rep = large_n2_report
    print_report(small_h4_report, title="H4 CAS(4,4) small")
    print_report(large_rep, title="N2 CASCI(8,10) large")
