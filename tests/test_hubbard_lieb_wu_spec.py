"""
Acceptance gates for specs/SPEC_hubbard_bethe.md §10 (bead chem-tjr): DMRG Hubbard chains at large
L against the exact Lieb-Wu thermodynamic-limit energy.

  G5  every vendored ladder passes the pre-registered point check (regime, undershoot vs spread,
      stage convergence) except the one recorded failure (open U=8 L=100), and the grid is complete.
  G6  acceptance (1) |e_a - lieb_wu| < 1e-3 Ha/site (asserted) and (2) route agreement (its
      measured outcome pinned -- it fails at U=2), both from the vendored table.
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
RECORDED = {
    2.0: (-0.8443500712762756, 1.587491516055195e-05, -0.8439654191321566, 8.536640116134714e-05),
    4.0: (-0.5737118083289061, 1.12570729054573e-05, -0.5735708668966306, 0.00013715000089707297),
    8.0: (-0.3275340436769117, 3.422257810632434e-05, -0.3275268404050134, 2.4467769166009023e-06),
}
# The one point that fails the pre-registered check, after its one allowed re-run (§10.4): its
# first ladder stage stalled. Pinned so the failure stays visible and nothing else may join it.
FAILED_POINTS = {("open", 8.0, 100)}
ROUTES_AGREE = {2.0: False, 4.0: True, 8.0: True}   # acceptance (2), as measured
TDL_GATE = 1e-3


def _analysis():
    import hubbard_tdl_analysis as an
    rows = an.load_table()
    return an, rows, an.analyse(rows)


def test_G5_every_vendored_point_passes_the_preregistered_check():
    _, rows, res = _analysis()
    assert rows, "empty table"
    bad = {(p[0], p[1], p[2]) for p in res["points"] if not p[3]}
    assert bad == FAILED_POINTS, bad
    for r in rows:
        if (r["route"], float(r["U"]), int(r["L"])) not in FAILED_POINTS:
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
        # (2) the two routes agree within their own bars -- recorded outcome pinned, not asserted
        # true: at U=2 it FAILS (rings L<=32 sit inside the charge crossover, 1/L^2 fit biased;
        # §10.4). A future table that closes the gap must update the spec to pass.
        assert r["pass_2"] is ROUTES_AGREE[U], (U, r["diff_ab"], r["a"][1], r["b"][1])


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
    """Through the production path (the driver's integrals(): folded ring order) at U=0, where the
    chain is gapless in both sectors -- the hardest case for a fixed D."""
    from benchmark_hubbard_lieb_wu import integrals
    from hybrid_quantum_solver.dmrg_reference import dmrg_energy_extrapolated
    for route in ("open", "ring"):
        h1, eri, ne, ec = integrals(route, 20, 0.0)
        e_free = 2.0 * float(np.sort(np.linalg.eigvalsh(h1))[:10].sum())
        r = dmrg_energy_extrapolated(h1, eri, ne, ec, bond_dims=(100, 200, 400),
                                     init_occs=np.ones(20), scratch="./.dmrg_tmp/g8_free")
        # The U=0 ring keeps a real truncation error at D=400 (6e-5 Ha, both sectors gapless), so
        # the anchor tests what the analysis relies on: E(D_max) is variational and the exact
        # energy lies within the pre-registered per-point sigma = max(stderr, |E_x - E(D_max)|).
        sigma = max(r.stderr, abs(r.energy - r.per_D[-1][2]))
        assert r.per_D[-1][2] >= e_free - 1e-9, (route, r.per_D, e_free)
        assert abs(r.energy - e_free) <= sigma, (route, r.energy, e_free, sigma)
        assert max(r.stage_dE) < 1e-6, r.stage_dE
    m = hubbard_chain_integrals(2, 4.0)
    r = dmrg_energy_extrapolated(m.h1, m.eri, m.nelec, m.e_core, bond_dims=(4, 4),
                                 init_occs=np.ones(2), scratch="./.dmrg_tmp/g8_dimer")
    assert abs(r.per_D[-1][2] - hubbard_dimer_energy(1.0, 4.0)) < 1e-8, r.per_D
