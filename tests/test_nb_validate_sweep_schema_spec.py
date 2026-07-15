"""
Acceptance gates G1-G4 for specs/SPEC_nb_validate_sweep_schema.md (nb_validate_sweep -- the
ragged-schema CSV union, checked not just asserted in a docstring).

Deliberately no new library code: every gate calls `nb_validate_sweep.py`'s existing
`validate_sweep`/`_to_frame` directly. Stays at the module's own small LiH CAS(2,2) __main__
example throughout (never raises the active space toward validate_and_cost's qubit_dense_max_orb
threshold -- see SPEC_validate_and_cost_composition.md's R2 OOM near-miss for why that matters).
"""
import polars as pl
import pytest
from pyscf import gto, mcscf, scf

from nb_validate_sweep import _to_frame, validate_sweep
from run_nbn_sqd_sweep import valid_spin_sectors
from validate_and_cost import _HAVE_FT

ATOM = "Li 0 0 0; H 0 0 1.6"
BASIS = "sto-3g"


def _independent_casci(mol_spin):
    mol = gto.M(atom=ATOM, basis=BASIS, spin=mol_spin)
    mf = scf.UHF(mol) if mol_spin != 0 else scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    cas = mcscf.CASCI(mf, 2, 2)
    cas.verbose = 0
    cas.kernel()
    return float(cas.e_tot)


@pytest.fixture(scope="module")
def sweep_rows(tmp_path_factory):
    out = tmp_path_factory.mktemp("nb_sweep") / "sector_verdicts.csv"
    rows = validate_sweep(ATOM, BASIS, {}, cas_electrons=2, cas_orbitals=2, output_csv=str(out))
    return rows, out


def test_G1_sweep_produces_correct_sectors_with_correct_energies(sweep_rows):
    """The sweep produces exactly the two physically reachable spin sectors (singlet, triplet),
    both OK, with reference_energy matching independently-computed CASCI."""
    assert valid_spin_sectors(2, 2) == [0, 2]
    rows, _out = sweep_rows
    assert len(rows) == 2, rows
    assert all(r["status"] == "OK" for r in rows), rows

    by_spin = {r["mol_spin"]: r for r in rows}
    assert abs(by_spin[0]["reference_energy"] - _independent_casci(0)) < 1e-6
    assert abs(by_spin[2]["reference_energy"] - _independent_casci(2)) < 1e-6


def test_G2_graceful_ft_degradation_through_the_csv_flattening_layer(sweep_rows):
    """Every OK row has ft_status == "no_openfermion" and every ft_* numeric field is None in the
    standard chem env -- the graceful-degradation path through the CSV-flattening layer, not just
    validate_and_cost's own return dict."""
    assert _HAVE_FT is False, "expected no openfermion in the standard chem env"
    rows, _out = sweep_rows
    for r in rows:
        assert r["ft_status"] == "no_openfermion", r
        for key in ("ft_threshold", "ft_lambda_DF", "ft_ccsd_t_err_mHa", "ft_toffoli",
                    "ft_logical_qubits"):
            assert r[key] is None, (key, r)


def test_G3_ragged_schema_union_matches_the_docstring_claim():
    """THE FINDING / definition of done: _to_frame on a synthetic [OK-shaped row, FAILED-shaped
    row] list produces a DataFrame whose columns are the union of both rows' keys, in FIRST-SEEN
    order, with the FAILED row's missing fields null -- not absent, not crashed, not misaligned."""
    ok_row = {
        "mol_spin": 0, "multiplicity": "singlet", "n_alpha": 1, "n_beta": 1, "norb": 2,
        "reference_energy": -7.86, "cross_check_agree": True, "max_dev_mHa": 0.001,
        "methods_run": "CASCI+Krylov", "methods_skipped": "", "n_qubits_original": 4,
        "n_qubits_tapered": 2, "ft_threshold": None, "ft_lambda_DF": None,
        "ft_ccsd_t_err_mHa": None, "ft_toffoli": None, "ft_logical_qubits": None,
        "ft_status": "no_openfermion", "status": "OK",
    }
    failed_row = {"mol_spin": 2, "multiplicity": "triplet", "status": "FAILED: SCF did not converge"}

    df = _to_frame([ok_row, failed_row])

    expected_order = list(ok_row.keys())  # "status" already last in ok_row; failed_row adds none new
    assert df.columns == expected_order, df.columns

    assert df.row(0, named=True)["mol_spin"] == 0
    assert df.row(1, named=True)["mol_spin"] == 2
    failed = df.row(1, named=True)
    assert failed["status"] == "FAILED: SCF did not converge"
    for key in ("n_alpha", "n_beta", "norb", "reference_energy", "cross_check_agree",
                "methods_run", "n_qubits_original", "ft_status"):
        assert failed[key] is None, (key, failed[key])


def test_G4_sweep_writes_a_valid_rereadable_csv(sweep_rows):
    """The sweep actually writes a valid, re-readable CSV to disk -- the on-disk artifact the
    module's own docstring promises, not just the in-memory return value."""
    rows, out = sweep_rows
    assert out.exists()
    reread = pl.read_csv(str(out))
    assert reread.height == 2
    assert reread.columns == _to_frame(rows).columns
