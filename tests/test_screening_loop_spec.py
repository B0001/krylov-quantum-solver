"""
Acceptance gates for specs/SPEC_interval_dominance_screening.md (certified-bracket pruning).

The dominance arithmetic is sound and the loop never picks the wrong winner. What v1.0 got wrong
is the ECONOMICS, and the gates below encode the falsification (specs/README.md step 5):

  G1   winner correctness -- the soundness gate, and it holds on both libraries in both modes.
  G1b  THE INHERITED KILL: Temple's SELF mode (the only certificate a screening loop actually
       has, since an oracle E_1 means the candidate is already solved) is NOT rigorous -- measured
       up to 0.695 mHa ABOVE the true energy. Oracle mode never violates. Same failure class as
       SPEC_excited_state_certification's G4, one rung up.
  G2   FALSIFIED: >= 75% saving is reachable ONLY with the unsound bound (79.2%). With rigorous
       oracle brackets the saving is 62.5-66.7%. The headline number was bought with an invalid
       certificate.
  G3   v1.0's "strictly non-increasing staircase" is VACUOUS -- an active pool cannot grow by
       construction. Replaced by the ratio that actually governs pruning: spread/width at the
       cheap dimension.
  G4   near-degenerate candidates resolve without a false prune.

All deterministic (exact statevector, no RNG). PySCF/qiskit, no block2; `make gates` isolates it.
"""
import numpy as np

from screening_loop import (
    brute_force_sweep,
    h4_library,
    interval_dominance_sweep,
)
from temple_bounds import krylov_bracket

WIDE = np.linspace(0.8, 2.2, 8)
TIGHT = np.linspace(1.00, 1.07, 8)
_CACHE = {}


def _run(spacings, oracle, max_m=12):
    key = (tuple(np.round(spacings, 5)), oracle, max_m)
    if key not in _CACHE:
        lib = h4_library(spacings, oracle=oracle)
        win, met = interval_dominance_sweep(lib, max_m=max_m)
        ref = h4_library(spacings, oracle=oracle)
        bwin, bmet = brute_force_sweep(ref, target_m=max_m)
        _CACHE[key] = (win, met, bwin, bmet)
    return _CACHE[key]


def test_G1_pruning_never_loses_the_true_winner():
    """The soundness gate: the adaptive screen returns the brute-force winner, both modes.

    Measured on H4 chains: WIDE (a = 0.8..2.2, spread 285 mHa) and TIGHT (a = 1.00..1.07, spread
    18.9 mHa), oracle and self mode -- 4/4 agreements, and the winner is never among the pruned.
    """
    for spacings in (WIDE, TIGHT):
        for oracle in (True, False):
            win, met, bwin, _ = _run(spacings, oracle)
            assert win.name == bwin.name, (spacings[0], oracle, win.name, bwin.name)
            assert win.active and win.pruned_at is None, win


def test_G1b_the_cheap_certificate_is_not_a_certificate():
    """THE INHERITED KILL: Temple self mode is unsound, so 'certified' pruning is not certified.

    Across the WIDE library at M = 2..16, the self-mode lower bound exceeds the exact FCI energy
    on 5 (candidate, M) pairs, worst by 0.695 mHa -- it is licence to prune a candidate that is
    still viable. The oracle bound (eps = exact E_1) violates on none. An E_1 oracle, however,
    presupposes the candidate is already solved, which is the thing screening exists to avoid.

    On this HOMOGENEOUS library no false prune results: every candidate's bound is inflated
    together, so the winner's margin L_winner - U_best stays negative (measured <= -0.43 mHa).
    That is a measured property of same-molecule libraries, NOT a guarantee -- a heterogeneous
    library whose members converge at different rates has no such common-mode protection.
    """
    lib = h4_library(WIDE, oracle=False)
    orc = h4_library(WIDE, oracle=True)
    self_viol, oracle_viol, worst = 0, 0, 0.0
    for cs, co in zip(lib, orc):
        exact = cs.mh.ground_state_energy()
        for m in (2, 4, 6, 8, 12, 16):
            bs = krylov_bracket(cs.mh, m, eps=None)
            bo = krylov_bracket(co.mh, m, eps=co.eps)
            if np.isfinite(bs.lower) and bs.lower - exact > 1e-9:
                self_viol += 1
                worst = max(worst, bs.lower - exact)
            if np.isfinite(bo.lower) and bo.lower - exact > 1e-9:
                oracle_viol += 1
    assert self_viol >= 3, self_viol                     # the bound is not a bound
    assert worst > 5e-4, worst                           # 0.695 mHa: far above arithmetic noise
    assert oracle_viol == 0, oracle_viol                 # the rigorous mode is rigorous

    win, met, _, _ = _run(WIDE, False)
    assert win.pruned_at is None                         # common-mode, so no false prune HERE


def test_G2_the_75_percent_saving_needs_the_unsound_bound():
    """FALSIFIES v1.0's economics: rigorous pruning saves 62.5-66.7%, not >= 75%.

    Basis-vector cost, 8 candidates, max_m = 12 (brute force = 96). Oracle (sound): WIDE 36 = 62.5%,
    TIGHT 32 = 66.7% -- both real, both under the specced 75%. Self mode reaches 79.2% on WIDE, but
    G1b shows that bound is invalid, so the only schedule that clears v1.0's bar is the one that is
    not certified. Gated as a ceiling on the SOUND saving so the kill cannot decay into a pass.
    """
    sound = [_run(s, True)[1]["saving"] for s in (WIDE, TIGHT)]
    assert all(0.50 < s < 0.75 for s in sound), sound    # real, but short of the target
    unsound = _run(WIDE, False)[1]["saving"]
    assert unsound > 0.75, unsound                       # the target is met only unsoundly
    assert unsound > max(sound), (unsound, sound)


def test_G3_pruning_power_is_the_spread_over_width_ratio():
    """Replaces v1.0's vacuous staircase gate with the quantity that actually decides savings.

    An active pool cannot grow, so "strictly non-increasing" is true by construction and gates
    nothing. What governs pruning is the library's energy spread over the certified bracket width
    at the CHEAP dimension. Measured at M = 2 on H4: WIDE spread/width = 4.45 prunes 6/8
    immediately; TIGHT spread/width = 1.00 prunes 0/8 and must reach M = 4 before anything moves.
    Same loop, same budget -- the saving is a property of the LIBRARY, not of the method.
    """
    ratios, pruned = [], []
    for spacings in (WIDE, TIGHT):
        lib = h4_library(spacings, oracle=False)
        for c in lib:
            c.refine(2)
        uppers = np.array([c.upper for c in lib])
        widths = np.array([c.upper - c.lower for c in lib])
        fin = widths[np.isfinite(widths)]
        ratios.append(float(uppers.max() - uppers.min()) / float(np.median(fin)))
        u_best = uppers.min()
        pruned.append(sum(1 for c in lib if np.isfinite(c.lower) and c.lower > u_best))

    assert ratios[0] > 3.0 > 1.5 > ratios[1], ratios     # WIDE ~4.45, TIGHT ~1.00
    assert pruned[0] >= 5, pruned                        # WIDE prunes at once
    assert pruned[1] == 0, pruned                        # TIGHT cannot prune at all at M=2
    for spacings in (WIDE, TIGHT):                       # and the pool never grows
        curve = [n for _, n in _run(spacings, True)[1]["pool_curve"]]
        assert curve == sorted(curve, reverse=True), curve


def test_G4_exactly_degenerate_candidates_resolve_without_a_false_prune():
    """A true tie must run both twins to max_m; a RESOLVABLE near-tie may legitimately prune one.

    Found while gating: twins 1e-5 A apart (0.02 mHa) are NOT a tie for this loop -- by M = 10 the
    oracle bracket is 0.0012 mHa wide, dominance is proven honestly, and the loop stops early at
    M = 10. That is correct behaviour, so the gate uses EXACTLY degenerate geometries, where no
    bracket can ever separate the pair. Both then reach max_m, nothing is pruned between them, and
    the winner is still a brute-force winner. The finding: near-degeneracy is resolved by
    CONVERGING, and only an exact tie exhausts the budget.
    """
    a = float(WIDE[0])
    geoms = [a, a, 1.4, 1.8, 2.2]
    lib = h4_library(geoms, oracle=True)
    twins = lib[:2]
    assert abs(twins[0].mh.ground_state_energy() - twins[1].mh.ground_state_energy()) < 1e-12

    win, met = interval_dominance_sweep(lib, max_m=12)
    ref = h4_library(geoms, oracle=True)
    bwin, _ = brute_force_sweep(ref, target_m=12)
    assert win.name == bwin.name, (win.name, bwin.name)
    for c in twins:
        assert c.pruned_at is None, (c.name, c.pruned_at)   # an exact tie never prunes
        assert c.current_m == 12, (c.name, c.current_m)     # both scaled to max_m
    assert met["pruned"] == 3, met["pruned"]                # the three non-twins do prune


def test_G5_a_vacuous_bound_never_prunes():
    """Temple goes vacuous (-inf) at some dimensions; that must block pruning, not enable it.

    In the WIDE sweep the Temple bound is -inf for 2-3 of 8 candidates at M = 6-8 (eps <= theta_0).
    A -inf lower bound can never exceed any finite upper bound, so those candidates survive -- the
    conservative direction. This pins that the loop reads a vacuous bound as "unknown", which is
    why pruning power is non-monotonic in M while the pool itself still never grows.
    """
    lib = h4_library(WIDE, oracle=False)
    seen_vacuous = False
    for m in (2, 4, 6, 8):
        for c in lib:
            c.refine(m)
        u_best = min(c.upper for c in lib)
        for c in lib:
            if not np.isfinite(c.lower):
                seen_vacuous = True
                assert not (c.lower > u_best)            # -inf never prunes
    assert seen_vacuous, "expected at least one vacuous Temple bound in the WIDE sweep"
