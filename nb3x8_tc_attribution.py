#!/usr/bin/env python3
"""
bead chem-phh/chem-b0y -- what explains the Nb3X8 5.3x (Cl) / 2.3x (Br) Tc overprediction, now that
the coordination/mean-field attribution is KILLED (chem-g78, nb3x8_gaps.py /
tests/test_nb3x8_gaps_spec.py::test_G8_*)?

``SPEC_nb3x8_magnetometry.md`` section 7 names three candidates the isolated-dimer model cannot see
and chem-g78 did not test:
  (1) the cooperative/first-order structural transition itself,
  (2) in-plane kagome exchange,
  (3) a direct lattice renormalization of t_s_perp.

This module pre-registers a killable check for each. (2) and (3) are computable from data vendored
or sourced into this repo; only (1) is still genuinely blocked (needs elastic/lattice free energies
nobody has measured or computed for this family).

CANDIDATE 3 -- THE EXECUTED CHECK: hold U0, Us fixed (the on-site terms, set by the local Nb3 cluster
chemistry) and ask how much the interlayer hopping |t_s_perp| alone would have to shrink for the
isolated-dimer T_max to land on the measured Tc. Compare that required shrinkage to an *independent,
already-tabulated* empirical anchor: how much t_s_perp actually differs between the bulk and bilayer
stacking variants of the SAME dimerized (LT) structural motif, across the whole halide family (Table
I already carries both). That spread (<=2.04% family-wide) is a real, measured bound on how much this
hopping moves between closely related variants of the same phase -- not a parameter nobody has fit,
a second independent slice of the same vendored ab-initio table.

FINDING: the needed shrinkage is 57.7% (Cl) / 35.8% (Br) -- 28-600x LARGER than the 0.09%-2.04%
spread that same family's own data shows between bulk and bilayer stacking of the SAME dimerized
motif. A same-motif renormalization of t_s_perp alone, at the scale set by this repo's own data for
how much that hopping actually varies between closely related dimerized structures, CANNOT explain
the overprediction -- this candidate is KILLED in the narrow form stated. (It is NOT the same claim as
candidate (1): the HT/LT hopping ratio -- 0.126 (Cl) / 0.121 (Br) -- shows |t_s_perp| CAN differ by a
large factor when the structure actually changes PHASE, not just stacking; a cooperative transition
that partially relaxes the lattice toward the HT structure is a distinct, untested mechanism that this
check does not rule out.)

CANDIDATE 2 -- THE EXECUTED CHECK (chem-b0y): SPEC_nb3x8_magnetometry.md G3(b) shows -J_interlayer/4
overshoots the measured (HT) Curie-Weiss theta_W by 3.75x (Cl), so SOME weaker coupling sets theta_W.
``in_plane_kagome_exchange`` conjectured this is a frustrated in-plane kagome exchange, and that
standard kagome frustration physics could suppress Tc below the naive single-dimer scale by a factor
growing with J_inplane/J_interlayer.

Sourced in-plane coupling: Zhang, Dai, Wang, Cheng & Ji, "Robust Mottness and tunable interlayer
magnetism in Nb3X8 (X = F, Cl, Br, I) bilayers" (arXiv:2508.17352), Table 1, DFT+U total-energy
mapping of in-plane NN/NNN exchange between Nb3 trimers (meV/Nb3, "positive = AFM" convention,
verified against the downloaded PDF text, not a search summary): Nb3Cl8 J1=+0.12, J2=-0.24;
Nb3Br8 J1=-2.74, J2=-0.37. Two findings already correct the premise before any suppression-factor
arithmetic: (a) the in-plane net is a TRIANGULAR lattice of trimers in this source, not a kagome net
-- "kagome exchange" in the conjecture's name is this repo's geometric shorthand for the trimer
lattice, not the edge-sharing-triangle kagome net frustration formulas assume; (b) the coupling is
not robustly antiferromagnetic -- Nb3Br8's J1 is FM by the paper's own sign convention, and Nb3Cl8's
J1 is a tiny near-degenerate AFM/FM difference (the paper: "Nb3Cl8 has nearly degenerate FM, ZZ-AFM
and Stripe-AFM configurations"). Frustration-driven Tc suppression requires a dominant, genuinely
antiferromagnetic coupling on a frustrated net; neither halide's sourced coupling satisfies that.

The quantitative check run anyway (steelmanning the conjecture past its own premise problem): the
standard RPA/mean-field self-consistency closure for weakly-coupled quantum dimers, J'_eff *
chi_dimer(T_c) = 1 (same structure that produces the Curie-Weiss law; chi_dimer from the already-
validated ``nb3x8_susceptibility.susceptibility``), gives the EFFECTIVE in-plane coupling that WOULD
be needed to pull the isolated-dimer crossover down to the measured Tc. Compared against the most
generous possible reading of the sourced J1/J2 -- full triangular-lattice coordination (z1=z2=6) on
the absolute values, i.e. ``sourced_inplane_bound`` = 6*|J1| + 6*|J2| -- the needed J'_eff exceeds
even that generous bound by ~9160x (Cl8) / ~36x (Br8). KILLED, well outside the margin the conjecture
asked for ("requires J' unphysically close to or exceeding J_interlayer" -- here it is 2-3 ORDERS OF
MAGNITUDE beyond J_interlayer itself, let alone the sourced in-plane scale).
"""
from __future__ import annotations

from scipy.optimize import brentq

from nb3x8_gaps import NB3X8_CLUSTERS, NB3X8_LT_BULK
from nb3x8_magnetometry import EXPERIMENT, MEV_PER_K, predicted_transition_K
from nb3x8_susceptibility import susceptibility
from odmd_spin import dimer_exchange_analytic

# In-plane NN (J1) / NNN (J2) exchange between Nb3 trimers, monolayer DFT+U total-energy mapping,
# meV per Nb3 trimer, "positive = AFM" sign convention. Zhang, Dai, Wang, Cheng & Ji, "Robust
# Mottness and tunable interlayer magnetism in Nb3X8 (X = F, Cl, Br, I) bilayers", arXiv:2508.17352,
# Table 1 -- verified against the downloaded PDF text (chem-jiy methodology), not a search summary.
# The paper's own in-plane lattice for the Nb3 trimers is TRIANGULAR, not kagome (see module
# docstring); only Cl/Br are used here (the two halides with a measured Tc to attribute).
NB3X8_INPLANE_EXCHANGE = {
    "Nb3Cl8": dict(J1=0.12, J2=-0.24, source="arXiv:2508.17352 Table 1"),
    "Nb3Br8": dict(J1=-2.74, J2=-0.37, source="arXiv:2508.17352 Table 1"),
}

# Same dimerized (LT) structural motif, two stacking variants, per halide -- both already in
# NB3X8_CLUSTERS (arXiv:2501.10320 Table I). Used ONLY as an empirical anchor for "how much can
# t_s_perp plausibly move without a phase change"; no new parameters.
_BULK_BIL_KEYS = {
    "Nb3F8": ("F  LT-bulk", "F  LT-bil"),
    "Nb3Cl8": ("Cl LT-bulk", "Cl LT-bil"),
    "Nb3Br8": ("Br LT-bulk", "Br LT-bil"),
    "Nb3I8": ("I  LT-bulk", "I  LT-bil"),
}

# chem-g78's open follow-up candidates this module does NOT test, with the pre-registered check
# that WOULD kill or support each once the data exists, and why it cannot be run now.
UNTESTED_CANDIDATES = {
    "cooperative_structural_transition": dict(
        description=(
            "Nb3Cl8/Nb3Br8's singlet transition is first-order (SPEC_nb3x8_susceptibility.md, "
            "SPEC_nb3x8_thermo.md: 'real Nb3Cl8 has a first-order one at ~90 K'), not the smooth "
            "Bleaney-Bowers crossover the isolated-dimer model assumes. A first-order transition is "
            "set by a free-energy crossing between the HT (undimerized) and LT (dimerized) phases, "
            "not by chi(T)'s maximum -- the lattice only ACQUIRES the large LT-phase J once ordered, "
            "so the naive 0.625*J(LT) estimate can overshoot the true crossing temperature."
        ),
        pre_registered_check=(
            "Build F_phase(T) = E_mag,phase(T) - T*S_mag,phase(T) + E_lattice,phase for phase in "
            "{HT, LT}, reusing the already-validated nb3x8_thermo Boltzmann-trace machinery for "
            "E_mag/S_mag with NB3X8_HT_BULK vs NB3X8_LT_BULK; solve F_HT(Tc)=F_LT(Tc) for Tc. KILL "
            "if, for every E_lattice value in the physically plausible range (van der Waals "
            "interlayer relaxation energies), Tc cannot be pulled down into the measured range "
            "(within the same factor-of-2 band G2 of SPEC_nb3x8_magnetometry.md already uses); "
            "SURVIVES if a plausible E_lattice reproduces measured Tc without extreme tuning."
        ),
        blocked_on=(
            "No elastic/lattice free-energy data for Nb3X8 is vendored anywhere in this repo -- "
            "arXiv:2501.10320 (the only ab-initio source used by nb3x8_gaps.py/nb3x8_susceptibility.py) "
            "reports electronic cRPA parameters only, no phonon/elastic constants. This sandbox has no "
            "working internet access to source E_lattice from a different paper (reconfirmed this "
            "session; see specs/BACKLOG.md note under chem-phh)."
        ),
        follow_up_bead="chem-cfl",
    ),
}


def same_motif_t_spread(name: str) -> float:
    """Relative |t_s_perp| difference between the bulk and bilayer stacking variants of the SAME
    dimerized (LT) structural motif -- an empirical, already-vendored bound on how much this hopping
    can move without a phase change (zero new parameters; both keys are in ``NB3X8_CLUSTERS``)."""
    bulk_key, bil_key = _BULK_BIL_KEYS[name]
    t_bulk = NB3X8_CLUSTERS[bulk_key]["t"]
    t_bil = NB3X8_CLUSTERS[bil_key]["t"]
    return abs(t_bil - t_bulk) / abs(t_bulk)


def max_family_same_motif_t_spread() -> float:
    """The largest same-motif |t_s_perp| spread across the whole F/Cl/Br/I family."""
    return max(same_motif_t_spread(name) for name in _BULK_BIL_KEYS)


def ht_lt_t_ratio(name: str) -> float:
    """|t_s_perp| ratio between the HT (undimerized) and LT (dimerized) PHASES -- context only, not
    part of the candidate-3 kill test (a phase change is candidate 1's territory, not candidate 3's).
    Only Cl/Br have a tabulated HT phase (Table IV)."""
    ht_key = {"Nb3Cl8": "Cl HT-bulk", "Nb3Br8": "Br HT-bulk"}[name]
    return NB3X8_CLUSTERS[ht_key]["t"] / NB3X8_LT_BULK[name]["t"]


def required_t_rescale(name: str, *, xtol: float = 1e-10) -> float:
    """The multiplier eta on t_s_perp (U0, Us held fixed) such that the isolated-dimer predicted
    transition scale ``predicted_transition_K(U0, eta*t, Us)`` equals the measured Tc. Solved by
    bisection on eta in (0, 1] (the model overpredicts, so eta < 1 is expected)."""
    p = NB3X8_LT_BULK[name]
    tc_measured = EXPERIMENT[name]["Tc_K"]

    def residual(eta: float) -> float:
        return predicted_transition_K(p["U0"], eta * p["t"], p["Us"]) - tc_measured

    return brentq(residual, 1e-4, 1.0, xtol=xtol)


def candidate3_verdict(name: str, *, margin: float = 10.0) -> dict:
    """Candidate 3 ('a direct lattice renormalization of t_s_perp') kill test for one halide.

    KILL bar: the rescaling needed to match measured Tc must stay within ``margin`` times the
    largest same-motif |t_s_perp| spread this family's own data shows (bulk vs bilayer, both
    dimerized) -- a 10x margin on top of an already-measured, independent quantity, not a number
    picked to force a result. If the needed rescaling blows past even that generous bar, a
    same-structural-motif renormalization of t_s_perp alone cannot be the explanation.
    """
    eta = required_t_rescale(name)
    needed_frac = abs(1.0 - eta)
    anchor = max_family_same_motif_t_spread()
    bound = margin * anchor
    return dict(
        name=name,
        eta=eta,
        needed_pct=100.0 * needed_frac,
        anchor_pct=100.0 * anchor,
        bound_pct=100.0 * bound,
        killed=needed_frac > bound,
    )


def required_inplane_coupling(name: str) -> float:
    """The effective in-plane coupling J'_eff (meV) that would need to satisfy the RPA/mean-field
    self-consistency closure ``J'_eff * chi_dimer(T_c) = 1`` for the isolated interlayer dimer's
    susceptibility (``nb3x8_susceptibility.susceptibility``, evaluated at the measured Tc in meV) to
    be pulled down to order the measured transition -- i.e. how strong an in-plane coupling candidate
    2 needs to do the suppression job entirely by itself."""
    p = NB3X8_LT_BULK[name]
    tc_meV = EXPERIMENT[name]["Tc_K"] * MEV_PER_K
    chi = susceptibility(p["U0"], p["t"], p["Us"], tc_meV)
    return 1.0 / chi


def sourced_inplane_bound(name: str, *, z1: float = 6.0, z2: float = 6.0) -> float:
    """The most generous possible reading of the sourced in-plane J1 (NN) / J2 (NNN) exchange
    (``NB3X8_INPLANE_EXCHANGE``, arXiv:2508.17352 Table 1) as an effective coupling scale: full
    triangular-lattice coordination (z1=z2=6) on the absolute values, summed with no cancellation.
    This deliberately overstates the physically expected in-plane coupling scale (ignores sign,
    ignores frustration, ignores that z=6 is every neighbor at once) to give candidate 2 the
    strongest possible case before the kill bar is applied."""
    c = NB3X8_INPLANE_EXCHANGE[name]
    return z1 * abs(c["J1"]) + z2 * abs(c["J2"])


def candidate2_verdict(name: str) -> dict:
    """Candidate 2 ('in-plane kagome exchange') kill test for one halide.

    KILL bar: the RPA-required effective in-plane coupling (``required_inplane_coupling``) must
    exceed the most generous sourced bound (``sourced_inplane_bound``) -- the pre-registered
    condition in ``UNTESTED_CANDIDATES`` (now resolved) asked for this once a real J' was sourced:
    "KILL if ... requires J' unphysically close to or exceeding J_interlayer". Here the needed J'
    exceeds even J_interlayer itself (not just the sourced in-plane scale) by 5.6x-299x, so the
    sourced-bound comparison is a fortiori.
    """
    needed = required_inplane_coupling(name)
    bound = sourced_inplane_bound(name)
    j_interlayer = dimer_exchange_analytic(**NB3X8_LT_BULK[name])
    return dict(
        name=name,
        needed_meV=needed,
        bound_meV=bound,
        j_interlayer_meV=j_interlayer,
        needed_over_bound=needed / bound,
        needed_over_j_interlayer=needed / j_interlayer,
        killed=needed > bound,
    )


if __name__ == "__main__":
    print("chem-phh/chem-b0y: Nb3X8 Tc-overprediction attribution -- three candidates, two tested\n")

    print("Candidate 2 (in-plane kagome/triangular exchange, sourced arXiv:2508.17352 Table 1):")
    for name in ("Nb3Cl8", "Nb3Br8"):
        c = NB3X8_INPLANE_EXCHANGE[name]
        v = candidate2_verdict(name)
        print(f"  {name}: sourced J1={c['J1']} meV, J2={c['J2']} meV ({c['source']})")
        print(f"    RPA-required J'_eff = {v['needed_meV']:.2f} meV "
              f"({v['needed_over_j_interlayer']:.2f}x J_interlayer={v['j_interlayer_meV']:.2f} meV)")
        print(f"    most generous sourced bound (z=6 each, |J1|+|J2|) = {v['bound_meV']:.2f} meV "
              f"-> needed exceeds bound by {v['needed_over_bound']:.1f}x "
              f"-> {'KILLED' if v['killed'] else 'survives'}")

    print("\nCandidate 3 (direct t_s_perp renormalization, same dimerized motif):")
    anchor = max_family_same_motif_t_spread()
    print(f"  empirical anchor (max bulk-vs-bilayer |t_s_perp| spread, family-wide): {anchor*100:.2f}%")
    for name in ("Nb3F8", "Nb3Cl8", "Nb3Br8", "Nb3I8"):
        print(f"    {name}: {same_motif_t_spread(name)*100:.3f}%")
    for name in ("Nb3Cl8", "Nb3Br8"):
        v = candidate3_verdict(name)
        print(f"\n  {name}: needs t_s_perp cut by {v['needed_pct']:.1f}% to match measured Tc "
              f"(eta={v['eta']:.4f}); kill bar = {v['bound_pct']:.1f}% (10x the {anchor*100:.2f}% "
              f"anchor) -> {'KILLED' if v['killed'] else 'survives'}")
        print(f"    context only (not part of the kill test): HT/LT |t_s_perp| ratio = "
              f"{ht_lt_t_ratio(name):.3f} (a full phase change, candidate 1's territory, moves t_s_perp "
              "by far more than a same-motif renormalization could)")

    print("\nCandidate 1: still BLOCKED on data this sandbox cannot source; see UNTESTED_CANDIDATES "
          "for the pre-registered check it would need to pass.")
    for key, v in UNTESTED_CANDIDATES.items():
        print(f"\n  [{key}] (follow-up: {v['follow_up_bead']})\n    blocked_on: {v['blocked_on']}")
