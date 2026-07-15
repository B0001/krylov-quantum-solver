"""
Acceptance gates G1-G4 for specs/SPEC_sqd_plumbing_verdict_coverage.md (SQD plumbing -- promoting
the smoke test, and exercising the verdicts it never triggers).

Deliberately no new library code: G1/G2/G4 mirror `smoke_test_sqd_plumbing.py`'s own calls into
`run_nbn_sqd_sweep.py`'s existing public functions (tightened where the original's assertions were
loose); G3 calls `validate_row`/`analyze_sector_trend` directly on constructed inputs, no SQD
execution needed for those branches.
"""
import numpy as np
from pyscf import gto, scf

import run_nbn_sqd_sweep as R


def _build_mf(atom, basis="sto-3g", spin=0):
    mol = gto.M(atom=atom, basis=basis, spin=spin)
    mf = scf.RHF(mol) if spin == 0 else scf.ROHF(mol)
    mf.verbose = 0
    mf.kernel()
    return mf


def test_G1_h2_full_space_reproduces_casci_pass_verdict():
    """H2 CAS(2,2) = full space: SQD over all determinants reproduces CASCI exactly, stays
    variational, and validate_row says PASS exactly."""
    rng = np.random.default_rng(0)
    mf = _build_mf("H 0 0 0; H 0 0 0.74")
    hcore, eri, e_core, nelec, casci = R.integrals_for_spin(mf, 2, 2, 0)
    norb = hcore.shape[0]
    ba = R.generate_bit_array_uniform(20_000, norb * 2, rand_seed=rng)
    best, _max_dim, _iters, _ = R.run_sqd_for_sector(hcore, eri, e_core, nelec, norb, ba,
                                                      samples_per_batch=100, spin_sq=0.0, rng=rng)
    err_mha = abs(best - casci) * 1e3
    assert err_mha < 1e-3, err_mha
    assert best >= casci - 1e-6, (best, casci)
    assert R.validate_row(best, casci) == "PASS"


def test_G2_frame_break_simulation_is_caught_exactly_as_frame_error():
    """THE (tightened) SMOKE-TEST PROMOTION: dropping e_core from a correct energy simulates the
    "507 Ha class of bug" -- must be caught. Tightened from the original's looser "one of
    FRAME_ERROR/ABOVE_TOL" disjunction: on this system e_core > 0, so dropping it makes the energy
    MORE NEGATIVE (below CASCI beyond VARIATIONAL_TOL_MHA), which validate_row's own logic maps to
    FRAME_ERROR specifically, not ABOVE_TOL."""
    rng = np.random.default_rng(0)
    mf = _build_mf("H 0 0 0; H 0 0 0.74")
    hcore, eri, e_core, nelec, casci = R.integrals_for_spin(mf, 2, 2, 0)
    norb = hcore.shape[0]
    ba = R.generate_bit_array_uniform(20_000, norb * 2, rand_seed=rng)
    best, _max_dim, _iters, _ = R.run_sqd_for_sector(hcore, eri, e_core, nelec, norb, ba,
                                                      samples_per_batch=100, spin_sq=0.0, rng=rng)
    assert e_core > 0, "test assumes a positive e_core on this system (checked, not assumed blindly)"
    broken = best - e_core
    assert R.validate_row(broken, casci) == "FRAME_ERROR"


def test_G3_every_validate_row_verdict_is_reachable():
    """THE FINDING / definition of done, part 1: validate_row's SKIP path -- untested by the
    original smoke test, which never passes a non-finite energy -- fires correctly for both a
    missing and a NaN energy."""
    assert R.validate_row(None, -1.0) == "SKIP"
    assert R.validate_row(float("nan"), -1.0) == "SKIP"
    assert R.validate_row(-1.0, None) == "SKIP"


def test_G3_every_analyze_sector_trend_verdict_is_reachable():
    """THE FINDING / definition of done, part 2: FLAT_SUBSPACE, STALLED, CONVERGING, and
    INSUFFICIENT -- four of analyze_sector_trend's five verdicts -- are never produced by any real
    run in this repo (only CONVERGED is, per G4). Constructed directly against the module's own
    published tolerance constants, all four fire correctly."""
    conv_tol = R.CONVERGENCE_TOL_MHA          # 10.0
    trend_min = R.TREND_MIN_IMPROVEMENT       # 0.20

    # FLAT_SUBSPACE: subspace_dim never grows, and the endpoint stays above CONVERGENCE_TOL_MHA.
    flat_rows = [
        {"delta_mHa": conv_tol + 10.0, "subspace_dim": 10},
        {"delta_mHa": conv_tol + 5.0, "subspace_dim": 10},
    ]
    assert R.analyze_sector_trend(flat_rows) == "FLAT_SUBSPACE"

    # STALLED: subspace grows, but improvement is below TREND_MIN_IMPROVEMENT, endpoint still
    # above CONVERGENCE_TOL_MHA.
    first_delta = conv_tol + 10.0
    stalled_last = first_delta * (1.0 - trend_min / 2)  # half the required improvement
    stalled_rows = [
        {"delta_mHa": first_delta, "subspace_dim": 10},
        {"delta_mHa": stalled_last, "subspace_dim": 50},
    ]
    assert stalled_last > conv_tol, "test construction sanity: endpoint must stay above tolerance"
    assert R.analyze_sector_trend(stalled_rows) == "STALLED"

    # CONVERGING: subspace grows, improvement clears TREND_MIN_IMPROVEMENT, but the endpoint is
    # still (just) above CONVERGENCE_TOL_MHA -- not yet CONVERGED.
    converging_last = conv_tol + 1.0
    assert converging_last <= (1.0 - trend_min) * first_delta, (
        "test construction sanity: must clear the improvement bar"
    )
    converging_rows = [
        {"delta_mHa": first_delta, "subspace_dim": 10},
        {"delta_mHa": converging_last, "subspace_dim": 50},
    ]
    assert R.analyze_sector_trend(converging_rows) == "CONVERGING"

    # INSUFFICIENT: fewer than 2 finite points.
    assert R.analyze_sector_trend([]) == "INSUFFICIENT"
    assert R.analyze_sector_trend([{"delta_mHa": 5.0, "subspace_dim": 10}]) == "INSUFFICIENT"


def test_G4_h4_sweep_variational_and_trend_asserted():
    """H4 CAS(4,4) samples_per_batch sweep: every point stays variational, and
    analyze_sector_trend on the REAL swept data is asserted CONVERGED -- the original smoke test
    computed and printed this value but never asserted on it."""
    rng = np.random.default_rng(0)
    mf = _build_mf("H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0")
    hcore, eri, e_core, nelec, casci = R.integrals_for_spin(mf, 4, 4, 0)
    norb = hcore.shape[0]
    ba = R.generate_bit_array_uniform(40_000, norb * 2, rand_seed=rng)
    rows = []
    for spb in (10, 30, 80, 200):
        best, max_dim, _iters, _ = R.run_sqd_for_sector(hcore, eri, e_core, nelec, norb, ba,
                                                         samples_per_batch=spb, spin_sq=0.0, rng=rng)
        assert best >= casci - 1e-6, (spb, best, casci)
        rows.append({"delta_mHa": abs(best - casci) * 1e3, "subspace_dim": max_dim})
    assert R.analyze_sector_trend(rows) == "CONVERGED", rows
