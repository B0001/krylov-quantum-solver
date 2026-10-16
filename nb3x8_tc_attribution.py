#!/usr/bin/env python3
"""
bead chem-phh -- what explains the Nb3X8 5.3x (Cl) / 2.3x (Br) Tc overprediction, now that the
coordination/mean-field attribution is KILLED (chem-g78, nb3x8_gaps.py / tests/test_nb3x8_gaps_spec.py
::test_G8_*)?

``SPEC_nb3x8_magnetometry.md`` section 7 names three candidates the isolated-dimer model cannot see
and chem-g78 did not test:
  (1) the cooperative/first-order structural transition itself,
  (2) in-plane kagome exchange,
  (3) a direct lattice renormalization of t_s_perp.

This module pre-registers a killable check for each. Only (3) is actually computable from data
already vendored in this repo (arXiv:2501.10320's Table I/IV cRPA parameters, via the already-
validated ``nb3x8_magnetometry``/``odmd_spin`` primitives) -- so it is the only one executed here.
(1) and (2) need ab-initio inputs (elastic/lattice free energies; in-plane kagome exchange) that
arXiv:2501.10320 does not report and that this sandbox cannot source (no working internet access --
reconfirmed this session, same finding as chem-a04's ``nb3x8_optical_phase_collapse_check.py``). They
are recorded as BLOCKED, not tested, with the check that would kill or support them written down now
so a future session with the data can run it without re-deriving the physics.

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
"""
from __future__ import annotations

from scipy.optimize import brentq

from nb3x8_gaps import NB3X8_CLUSTERS, NB3X8_LT_BULK
from nb3x8_magnetometry import EXPERIMENT, predicted_transition_K

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
    "in_plane_kagome_exchange": dict(
        description=(
            "The isolated bilayer dimer carries only the interlayer J; SPEC_nb3x8_magnetometry.md "
            "G3(b) already shows -J/4 overshoots the measured (HT) Curie-Weiss theta_W by 3.75x "
            "(Cl), so SOME weaker coupling (the in-plane kagome exchange) sets theta_W. If that "
            "in-plane coupling is antiferromagnetic and frustrated on the kagome net, standard "
            "frustration physics suppresses the 3-D ordering temperature below the naive "
            "single-dimer scale by a factor that grows with J_inplane/J_interlayer."
        ),
        pre_registered_check=(
            "Once an ab-initio or experimental (INS/EPR) in-plane J' is available: compute the "
            "frustration ratio f = J'/J_interlayer and compare to the needed suppression (5.3x Cl, "
            "2.3x Br, already measured by this module's sibling nb3x8_magnetometry.overprediction_"
            "factor). KILL if the frustration-suppression formulas standard for a kagome net cannot "
            "reach a 5.3x/2.3x suppression for a J' at the scale cRPA would plausibly give (i.e. "
            "requires J' unphysically close to or exceeding J_interlayer); SURVIVES if a modest "
            "J'/J_interlayer (<~0.3, typical of a weak superexchange path) reproduces it."
        ),
        blocked_on=(
            "No in-plane kagome exchange constant is tabulated anywhere in this repo or in "
            "arXiv:2501.10320's vendored tables (Table I/IV cover the interlayer dimer only); this "
            "sandbox has no working internet access to source one from a different ab-initio or "
            "experimental paper (reconfirmed this session; see specs/BACKLOG.md note under chem-phh)."
        ),
        follow_up_bead="chem-b0y",
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


if __name__ == "__main__":
    print("chem-phh: Nb3X8 Tc-overprediction attribution -- three candidates, one testable with "
          "vendored data\n")

    print("Candidate 3 (direct t_s_perp renormalization, same dimerized motif):")
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

    print("\nCandidates 1 & 2: BLOCKED on data this sandbox cannot source (no internet access, "
          "reconfirmed this session); see UNTESTED_CANDIDATES for the pre-registered check each "
          "would need to pass.")
    for key, v in UNTESTED_CANDIDATES.items():
        print(f"\n  [{key}] (follow-up: {v['follow_up_bead']})\n    blocked_on: {v['blocked_on']}")
