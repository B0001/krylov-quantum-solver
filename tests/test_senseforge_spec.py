"""
Acceptance gates for specs/tasks/04-senseforge.md (SenseForge M1-M2: sweep harness + certified
derivatives).

Two deviations from the literal task spec are gated as findings, not hidden -- see
senseforge/hamiltonian.py's module docstring for the full physics reasoning:

  (1) No ab-initio CIF-strain geometry path exists; strain eps is defined as the fractional
      hopping perturbation t(eps) = t0*(1+eps), reusing nb3x8_strain.py's own "|t| is the sole
      strain proxy" convention.
  (2) certified_gaps.gap_bracket silently returns the WRONG gap (bright ionic ~1117 meV, not the
      spin gap J~66 meV) on this system -- the singlet-triplet gap is dark to the HF reference.
      Every SenseForge gap is exact diagonalization / closed form instead (the system is
      FCI-trivial; certified_gaps.py's own docstring says certification is "pointless" here).

PySCF-free (pure numpy/qiskit-nature on a 4-qubit model); fast. `make gates` runs this in its
own process.
"""
from pathlib import Path

import numpy as np
import pytest

from certchem.contract import Bracket, Certificate, CertifiedResult
from nb3x8_gaps import NB3X8_LT_BULK
from odmd_spin import dimer_exchange_analytic
from senseforge.candidates import (
    RANKING_DEGENERATE,
    RANKING_INTERIOR,
    RANKING_MONOTONE,
    build_candidates,
    figure_of_merit,
    rank_candidates,
    ranking_verdict,
)
from senseforge.config import ConfigError, SweepConfig, load_config
from senseforge.hamiltonian import (
    G_FACTOR,
    MU_B_MEV_PER_T,
    certified_strain_gap,
    check_zeeman_hermitian,
    strained_params,
    zeeman_split_gap,
)
from senseforge.headers import ADR_0003_NOTE, has_adr0003_note
from senseforge.pipeline import run_pipeline
from senseforge.sensitivity import certified_central_differences
from senseforge.sweep import csv_path, read_sweep_csv, run_sweep
from senseforge.validation import cross_check_closed_form_vs_exact


# --- G1: config schema --------------------------------------------------------------------


def test_G1_config_loads_and_resolves_defaults(tmp_path):
    p = tmp_path / "cfg.yaml"
    p.write_text("halide: Cl\naxis: strain\n")
    cfg = load_config(str(p))
    assert cfg.system == "Nb3Cl8"
    assert cfg.grid_min == pytest.approx(-0.02)
    assert cfg.grid_max == pytest.approx(0.02)
    assert len(cfg.grid()) == 17  # [-2%, +2%] step 0.25%, PRD sec 3


def test_G1_bad_config_fails_with_named_field(tmp_path):
    p = tmp_path / "bad.yaml"
    p.write_text("halide: Xx\naxis: strain\n")
    with pytest.raises(ConfigError) as exc:
        load_config(str(p))
    assert exc.value.field == "halide"


def test_G1_config_hash_is_deterministic_and_sensitive():
    a = SweepConfig(halide="Cl", axis="strain", grid_min=-0.02, grid_max=0.02, grid_step=0.0025)
    b = SweepConfig(halide="Cl", axis="strain", grid_min=-0.02, grid_max=0.02, grid_step=0.0025)
    c = SweepConfig(halide="Br", axis="strain", grid_min=-0.02, grid_max=0.02, grid_step=0.0025)
    assert a.content_hash() == b.content_hash()
    assert a.content_hash() != c.content_hash()


# --- G2: strain + Zeeman physics ------------------------------------------------------------


def test_G2_strain_eps_zero_reproduces_validated_J():
    """t(eps=0) must be the unperturbed literature parameters, and the strain-gap machinery must
    reproduce the already-gated dimer_exchange_analytic exactly."""
    p = strained_params("Nb3Cl8", 0.0)
    assert p == NB3X8_LT_BULK["Nb3Cl8"]
    J = dimer_exchange_analytic(**NB3X8_LT_BULK["Nb3Cl8"])
    result = certified_strain_gap("Nb3Cl8", 0.0)
    assert result.bracket.best_estimate_hartree == pytest.approx(J, abs=1e-9)
    assert result.bracket.width == 0.0  # exact, by construction


def test_G2_strain_gap_direction_matches_gruneisen():
    """Positive eps increases |t|; nb3x8_strain.py's finding is gamma_J > 0 (spin gap stiffens
    under compression) for every halide -- so J(eps>0) > J(0) here too."""
    for system in NB3X8_LT_BULK:
        j0 = certified_strain_gap(system, 0.0).bracket.best_estimate_hartree
        j_plus = certified_strain_gap(system, 0.01).bracket.best_estimate_hartree
        assert j_plus > j0, system


def test_G2_zeeman_field_is_hermitian():
    check = check_zeeman_hermitian(B=7.3)
    assert check.is_hermitian
    assert check.max_asymmetry < 1e-10


def test_G2_zeeman_zero_field_matches_exact_J():
    J = dimer_exchange_analytic(**NB3X8_LT_BULK["Nb3Cl8"])
    result = zeeman_split_gap("Nb3Cl8", 0.0)
    assert result.bracket.best_estimate_hartree == pytest.approx(J, abs=1e-9)


def test_G2_zeeman_field_response_is_linear_below_crossing():
    """The finding (see hamiltonian.py docstring): a uniform field splits the triplet Sz=+/-1
    sublevels EXACTLY linearly (constant slope -g*mu_B) until the level-crossing field -- and
    that crossing sits at ~572 T for Nb3Cl8, far outside the PRD's [0, 10] T sweep, so the whole
    default sweep range is safely linear."""
    J = dimer_exchange_analytic(**NB3X8_LT_BULK["Nb3Cl8"])
    for B in (1.0, 5.0, 10.0):
        gap = zeeman_split_gap("Nb3Cl8", B).bracket.best_estimate_hartree
        expected = J - G_FACTOR * MU_B_MEV_PER_T * B
        assert gap == pytest.approx(expected, abs=1e-9)


def test_G2_zeeman_uniform_field_has_no_effect_within_hf_sector():
    """The physics finding that ruled out routing this through certified_gaps.gap_bracket: the
    HF-reachable (1,1)/Sz=0 sector's OWN excitation (E1-E0 within that sector alone) is exactly
    field-independent, because Sz_tot = 0 identically for every state in it. This is exactly
    what would make ``certified_gaps.gap_bracket`` (which only ever explores this sector, seeded
    from the closed-shell HF reference) blind to the field entirely."""
    from nb3x8_gaps import dimer_cluster_integrals
    from senseforge.hamiltonian import _N_OP, _SZ_OP

    p = NB3X8_LT_BULK["Nb3Cl8"]
    mh = dimer_cluster_integrals(**p).to_hamiltonian()
    H0 = mh.qubit_hamiltonian.to_matrix()
    N_op = _N_OP
    Sz = _SZ_OP.to_matrix()

    def two_lowest_in_n2_sz0_sector(H):
        w, V = np.linalg.eigh(H)
        n = np.real(np.einsum("ji,jk,ki->i", V.conj(), N_op, V))
        sz = np.real(np.einsum("ji,jk,ki->i", V.conj(), Sz, V))
        keep = (np.abs(n - 2.0) < 1e-8) & (np.abs(sz) < 1e-6)
        return np.sort(w[keep].real)[:2]

    e_plain = two_lowest_in_n2_sz0_sector(H0)
    e_field = two_lowest_in_n2_sz0_sector(H0 + 3.0 * Sz)
    assert np.allclose(e_plain, e_field, atol=1e-8)


# --- G3: sweep driver + resumability ---------------------------------------------------------


def test_G3_sweep_produces_one_row_per_grid_point(tmp_path):
    cfg = SweepConfig(halide="Cl", axis="strain", grid_min=-0.01, grid_max=0.01, grid_step=0.005,
                      output_dir=str(tmp_path / "out"))
    run_sweep(cfg)
    rows = read_sweep_csv(cfg)
    assert len(rows) == len(cfg.grid()) == 5


def test_G3_kill_mid_sweep_and_resume_gives_identical_csv(tmp_path):
    out_a = tmp_path / "full"
    out_b = tmp_path / "resumed"
    cfg_full = SweepConfig(halide="Cl", axis="strain", grid_min=-0.01, grid_max=0.01,
                           grid_step=0.005, output_dir=str(out_a))
    cfg_resumed = SweepConfig(halide="Cl", axis="strain", grid_min=-0.01, grid_max=0.01,
                              grid_step=0.005, output_dir=str(out_b))

    run_sweep(cfg_full)
    full_rows = read_sweep_csv(cfg_full)

    run_sweep(cfg_resumed, _limit=2)  # simulate a kill after 2 points
    with pytest.raises(ValueError):
        read_sweep_csv(cfg_resumed)  # incomplete
    run_sweep(cfg_resumed)  # resume
    resumed_rows = read_sweep_csv(cfg_resumed)

    assert [r["gap"] for r in full_rows] == [r["gap"] for r in resumed_rows]
    assert [r["x"] for r in full_rows] == [r["x"] for r in resumed_rows]


def test_G3_resume_does_not_recompute_existing_points(tmp_path, monkeypatch):
    cfg = SweepConfig(halide="Cl", axis="strain", grid_min=-0.01, grid_max=0.01, grid_step=0.005,
                      output_dir=str(tmp_path / "out"))
    run_sweep(cfg)

    calls = []
    import senseforge.sweep as sweep_mod
    real_compute = sweep_mod.compute_gap

    def spy(config, x):
        calls.append(x)
        return real_compute(config, x)

    monkeypatch.setattr(sweep_mod, "compute_gap", spy)
    run_sweep(cfg)  # fully cached -- should compute nothing new
    assert calls == []


# --- G4: ADR-0003 header automation -----------------------------------------------------------


def test_G4_every_artifact_carries_the_adr0003_note(tmp_path):
    cfg = SweepConfig(halide="Cl", axis="strain", grid_min=-0.01, grid_max=0.01, grid_step=0.005,
                      output_dir=str(tmp_path / "out"))
    report_path = run_pipeline(cfg)
    artifacts = [csv_path(cfg), report_path] + list(Path(cfg.output_dir).glob("design_card_*.md"))
    assert len(artifacts) >= 3
    for path in artifacts:
        assert has_adr0003_note(path.read_text()), path


def test_G4_header_note_is_the_exact_adr_text():
    assert ADR_0003_NOTE == "cluster-model prediction; not validated for the periodic solid."


# --- G5: certified finite differences ---------------------------------------------------------


def test_G5_synthetic_quadratic_slope_is_contained_in_the_propagated_bracket():
    def make_result(x, width=0.5):
        true = 2 * x**2 + 3 * x + 10
        return CertifiedResult(
            bracket=Bracket(lower_hartree=true - width, upper_hartree=true + width,
                            best_estimate_hartree=true),
            certificate=Certificate(method="synthetic", floor_check="n/a", krylov_dim=0,
                                    convergence="exact", solver_version="0"),
        )

    xs = [round(-1.0 + 0.1 * i, 6) for i in range(21)]
    results = [make_result(x) for x in xs]
    sens = certified_central_differences(xs, results)
    assert len(sens) == len(xs) - 2
    for s in sens:
        analytic_slope = 4 * s.x + 3
        assert s.slope_lower <= analytic_slope <= s.slope_upper
        assert s.second_derivative == pytest.approx(4.0, abs=1e-8)


def test_G5_requires_uniform_grid():
    xs = [0.0, 0.1, 0.5]  # non-uniform
    results = [CertifiedResult(bracket=Bracket(0, 0, 0),
                               certificate=Certificate("x", "n/a", 0, "exact", "0"))
              for _ in xs]
    with pytest.raises(ValueError):
        certified_central_differences(xs, results)


# --- G6: FoM ranking + candidates ---------------------------------------------------------------


def test_G6_fom_falls_back_to_raw_slope_at_zero_width():
    assert figure_of_merit(slope=2.5, bracket_width=0.0) == pytest.approx(2.5)
    assert figure_of_merit(slope=-2.5, bracket_width=0.0) == pytest.approx(2.5)
    assert figure_of_merit(slope=4.0, bracket_width=2.0) == pytest.approx(2.0)


def test_G6_ranking_is_descending_by_fom():
    from senseforge.candidates import Candidate

    cands = [
        Candidate(halide="Cl", axis="strain", x=0.0, gap=66.0, bracket_width=0.0, slope=1.0,
                  slope_lower=1.0, slope_upper=1.0, second_derivative=0.0, fom=1.0),
        Candidate(halide="Cl", axis="strain", x=0.01, gap=67.0, bracket_width=0.0, slope=5.0,
                  slope_lower=5.0, slope_upper=5.0, second_derivative=0.0, fom=5.0),
    ]
    ranked = rank_candidates(cands)
    assert [c.fom for c in ranked] == [5.0, 1.0]


def test_G6_every_candidate_row_carries_its_bracket(tmp_path):
    cfg = SweepConfig(halide="Cl", axis="strain", grid_min=-0.01, grid_max=0.01, grid_step=0.005,
                      output_dir=str(tmp_path / "out"))
    report_path = run_pipeline(cfg)
    text = report_path.read_text()
    assert "| rank |" in text
    # every data row (skip header/separator) must have the bracket column populated
    rows = [ln for ln in text.splitlines() if ln.startswith("| ") and "rank" not in ln
           and "---" not in ln]
    assert len(rows) > 0
    for row in rows:
        cols = [c.strip() for c in row.strip("|").split("|")]
        assert cols[6] != ""  # bracket column


# --- G7: first real sweep artifacts exist -----------------------------------------------------


def test_G7_nb3cl8_strain_sweep_config_is_committed_and_runs():
    cfg = load_config("senseforge/configs/nb3cl8_strain.yaml")
    assert cfg.system == "Nb3Cl8"
    assert cfg.axis == "strain"


def test_G7_first_real_sweep_artifacts_present_in_results():
    # Produced by `python run_senseforge_sweep.py senseforge/configs/nb3cl8_strain.yaml`.
    out = Path("results/senseforge/Nb3Cl8_strain")
    assert (out / "candidates.md").exists()
    assert (out / "gap_vs_strain.csv").exists()


# --- G8: Gate 1 (adapted) -- the substitute cross-check, and its honest verdict ---------------


@pytest.mark.parametrize("system", list(NB3X8_LT_BULK))
def test_G8_closed_form_matches_exact_diagonalization(system):
    """The substitute for literal cluster-size convergence (which does not apply -- see
    senseforge/validation.py): two independent exact computations of the same B=0 gap must
    agree to machine precision, on every halide in the validated family."""
    result = cross_check_closed_form_vs_exact(system)
    assert result.agrees, result
    assert result.discrepancy < 1e-8


# --- G9 (THE FINDING): the FoM ranking does not discriminate on this model ---------------------


def _candidates_for(axis, grid_min, grid_max, grid_step, tmp_path):
    cfg = SweepConfig(halide="Cl", axis=axis, grid_min=grid_min, grid_max=grid_max,
                      grid_step=grid_step, output_dir=str(tmp_path))
    run_sweep(cfg)
    rows = read_sweep_csv(cfg)  # csv.DictReader -> every value is a str
    xs = [float(r["x"]) for r in rows]
    gaps = [float(r["gap"]) for r in rows]
    widths = [float(r["upper"]) - float(r["lower"]) for r in rows]
    results = [
        CertifiedResult(
            bracket=Bracket(lower_hartree=float(r["lower"]), upper_hartree=float(r["upper"]),
                            best_estimate_hartree=float(r["gap"])),
            certificate=Certificate(method="exact", floor_check="n/a", krylov_dim=0,
                                    convergence="exact", solver_version="0"),
        )
        for r in rows
    ]
    sens = certified_central_differences(xs, results)
    return build_candidates("Cl", axis, xs, gaps, widths, sens)


def test_G9_field_axis_ranking_is_degenerate(tmp_path):
    """The Zeeman response is EXACTLY linear (G2), so d(gap)/dB is the same constant at every B:
    every operating point ties on FoM and the rank order is pure sort noise. The pre-fix
    design_card_1.md published '+2 T' as the top operating point -- a recommendation that was
    provably no better than any other point in the sweep."""
    cands = _candidates_for("field", 0.0, 10.0, 0.5, tmp_path)
    assert ranking_verdict(cands) == RANKING_DEGENERATE
    foms = [c.fom for c in cands]
    assert max(foms) - min(foms) < 1e-9 * max(foms), (min(foms), max(foms))


def test_G9_strain_axis_ranking_is_a_window_edge_not_an_optimum(tmp_path):
    """|S| is monotone across the strain window, so rank 1 is the WINDOW EDGE. Falsifiable
    directly: widen the window and the 'best' operating point moves with it -- a discovered
    optimum would not."""
    narrow = _candidates_for("strain", -0.02, 0.02, 0.0025, tmp_path / "narrow")
    assert ranking_verdict(narrow) == RANKING_MONOTONE
    best_narrow = rank_candidates(narrow)[0].x

    wide = _candidates_for("strain", -0.05, 0.05, 0.0025, tmp_path / "wide")
    best_wide = rank_candidates(wide)[0].x

    assert best_wide > best_narrow, (best_narrow, best_wide)   # the "optimum" chased the window
    assert best_narrow == pytest.approx(0.0175, abs=1e-6)
    assert best_wide == pytest.approx(0.0475, abs=1e-6)


def test_G9_no_axis_yields_an_interior_optimum(tmp_path):
    """The headline: on the only model this repo can reach, NEITHER axis produces a genuine
    interior optimum -- so SenseForge's screening premise (rank operating points by FoM) is
    vacuous here. This is the boundary, recorded as a gate so it cannot be quietly forgotten."""
    for axis, lo, hi, step in (("strain", -0.02, 0.02, 0.0025), ("field", 0.0, 10.0, 0.5)):
        cands = _candidates_for(axis, lo, hi, step, tmp_path / axis)
        assert ranking_verdict(cands) != RANKING_INTERIOR, axis


def test_G9_artifacts_disclose_the_non_discrimination(tmp_path):
    """Every published artifact must SAY so -- a reader must not mistake a tie or a window edge
    for a result. The design cards are recommendations; they may not imply one was made."""
    for axis, marker in (("field", "NO DISCRIMINATION"), ("strain", "NO INTERIOR OPTIMUM")):
        cfg = SweepConfig(halide="Cl", axis=axis,
                          **({"grid_min": 0.0, "grid_max": 10.0, "grid_step": 0.5} if axis == "field"
                             else {"grid_min": -0.02, "grid_max": 0.02, "grid_step": 0.0025}),
                          output_dir=str(tmp_path / axis))
        report = run_pipeline(cfg)
        assert marker in report.read_text(), axis
        for card in Path(cfg.output_dir).glob("design_card_*.md"):
            assert marker in card.read_text(), (axis, card)


def test_G9_krylov_dim_is_gone_from_config_and_artifacts(tmp_path):
    """The dead field that wrote FALSE PROVENANCE: nothing read krylov_dim (every Certificate
    here is built with krylov_dim=0 -- the gaps are exact diagonalization, no Krylov subspace is
    ever constructed), yet it was hashed into content_hash and stamped on every design card as
    'krylov_dim=12'. Removed; a config that still sets it now fails with the named field."""
    assert "krylov_dim" not in SweepConfig.__dataclass_fields__

    p = tmp_path / "stale.yaml"
    p.write_text("halide: Cl\naxis: strain\nkrylov_dim: 12\n")
    with pytest.raises(ConfigError) as exc:
        load_config(str(p))
    assert exc.value.field == "krylov_dim"

    # NB: assert on the stamped "krylov_dim=" key, not the bare word -- pytest names tmp_path
    # after the test, so the word appears in the (legitimately stamped) output_dir path.
    cfg = SweepConfig(halide="Cl", axis="strain", grid_min=-0.01, grid_max=0.01, grid_step=0.005,
                      output_dir=str(tmp_path / "out"))
    report = run_pipeline(cfg)
    assert "krylov_dim=" not in report.read_text()
