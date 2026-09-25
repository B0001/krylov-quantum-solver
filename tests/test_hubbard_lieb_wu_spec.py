"""
Acceptance gates for specs/SPEC_hubbard_bethe.md §10 (bead chem-tjr): DMRG Hubbard chains at large
L against the exact Lieb-Wu thermodynamic-limit energy.

  G5  every vendored ladder passes the pre-registered point check (regime, undershoot vs spread,
      stage convergence), and the recorded headline is reproduced by hubbard_tdl_analysis.
  G6  acceptance (1) |e_a - lieb_wu| < 1e-3 Ha/site and (2) route agreement, pinned to the table.
  G7  pure: the §10.1 stalled ladder, and why `regime` alone cannot be trusted on it.
  G8  block2 anchors: DMRG with init_occs reproduces the analytic free-fermion and dimer energies.

G5-G7 are pure (csv + numpy). G8 needs block2 and no pyscf -- this file never imports pyscf, so it
is safe in its own process (CLAUDE.md segfault note).
"""
import numpy as np
import pytest

from hybrid_quantum_solver.dmrg_reference import dmrg_available, truncation_regime
from hybrid_quantum_solver.model_hamiltonians import (
    hubbard_chain_integrals,
    hubbard_dimer_energy,
    lieb_wu_energy,
)

requires_dmrg = pytest.mark.skipif(not dmrg_available(), reason="block2 not installed")

# Filled in from the vendored table; if the table or the analysis changes, G6 fails until the spec
# is updated with it. Keys: U -> (e_a, bar_a, e_b, bar_b).
RECORDED = {}
TDL_GATE = 1e-3


def _analysis():
    import hubbard_tdl_analysis as an
    rows = an.load_table()
    return an, rows, an.analyse(rows)


def test_G5_every_vendored_point_passes_the_preregistered_check():
    _, rows, res = _analysis()
    assert rows, "empty table"
    bad = [p for p in res["points"] if not p[3]]
    assert not bad, bad
    for r in rows:
        assert r["regime"] in ("converged", "truncation"), (r["route"], r["U"], r["L"])


def test_G5_table_covers_the_preregistered_grid():
    _, rows, _ = _analysis()
    have = {(r["route"], float(r["U"]), int(r["L"])) for r in rows}
    for U in (2.0, 4.0, 8.0):
        for L in (20, 40, 60, 80, 100):
            assert ("open", U, L) in have, ("open", U, L)
        for L in (16, 20, 24, 28, 32):
            assert ("ring", U, L) in have, ("ring", U, L)


def test_G6_recorded_tdl_numbers_and_acceptance():
    _, _, res = _analysis()
    for U, (e_a, bar_a, e_b, bar_b) in RECORDED.items():
        r = res["U"][U]
        assert abs(r["a"][0] - e_a) < 1e-9 and abs(r["a"][1] - bar_a) < 1e-9, (U, r["a"])
        assert abs(r["b"][0] - e_b) < 1e-9 and abs(r["b"][1] - bar_b) < 1e-9, (U, r["b"])
    assert set(RECORDED) == {2.0, 4.0, 8.0}
    for U in (2.0, 4.0, 8.0):
        r = res["U"][U]
        # (1) headline: route (a) vs the exact Bethe-ansatz integral
        assert abs(r["a"][0] - lieb_wu_energy(U)) < TDL_GATE, (U, r["resid_a"])
        # (2) the two routes agree within their own bars
        assert r["diff_ab"] < r["a"][1] + r["b"][1], (U, r["diff_ab"], r["a"][1], r["b"][1])


def test_G7_stalled_ladder_is_labelled_converged_but_fails_the_spread_free_checks():
    """§10.1: open L=60, U=4 from block2's default random MPS. Weights all <= 4e-9, so `regime`
    reads "converged" -- the stall that chem-4e9 targets, missed on the 1/D branch. Pinned so a
    classifier fix has to update this spec. The ladder's energies are 13 Ha apart; the spec's
    stage_dE check (not reproducible from these triples) is what rejects such a run."""
    stalled = [(100, 3.524227915646634e-09, -20.780402942374906),
               (200, 5.688756095253755e-10, -30.024624553972146),
               (400, 1.2913782637611208e-12, -34.04430453150846)]
    assert truncation_regime(stalled) == "converged"


@requires_dmrg
def test_G8_dmrg_free_fermion_and_dimer_anchors():
    from hybrid_quantum_solver.dmrg_reference import dmrg_energy_extrapolated
    for open_chain in (True, False):
        m = hubbard_chain_integrals(20, 0.0, open_chain=open_chain)
        e_free = 2.0 * float(np.sort(np.linalg.eigvalsh(m.h1))[:10].sum())
        r = dmrg_energy_extrapolated(m.h1, m.eri, m.nelec, m.e_core, bond_dims=(100, 200, 400),
                                     init_occs=np.ones(20), scratch="./.dmrg_tmp/g8_free")
        assert abs(r.per_D[-1][2] - e_free) < 1e-6, (open_chain, r.per_D, e_free)
        assert max(r.stage_dE) < 1e-6, r.stage_dE
    m = hubbard_chain_integrals(2, 4.0)
    r = dmrg_energy_extrapolated(m.h1, m.eri, m.nelec, m.e_core, bond_dims=(4, 4),
                                 init_occs=np.ones(2), scratch="./.dmrg_tmp/g8_dimer")
    assert abs(r.per_D[-1][2] - hubbard_dimer_energy(1.0, 4.0)) < 1e-8, r.per_D
