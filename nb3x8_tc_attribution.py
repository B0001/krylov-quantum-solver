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

This module pre-registers a killable check for each. (3) was executed first (chem-phh), from data
already vendored in this repo (arXiv:2501.10320's Table I/IV cRPA parameters, via the already-
validated ``nb3x8_magnetometry``/``odmd_spin`` primitives). (1) (chem-cfl) is now also executed,
using elastic/lattice data sourced from a SECOND primary reference this sandbox could fetch over the
network this session (arXiv:2412.03360 -- see below); the "no working internet access" finding
recorded by chem-phh/chem-a04 (``WebSearch`` only) did not hold for ``WebFetch`` on direct URLs,
including the arXiv API and PDF URLs. (2) still needs an in-plane kagome exchange constant that
neither arXiv:2501.10320 nor arXiv:2412.03360 reports, so it remains BLOCKED, with the check that
would kill or support it written down now so a future session with the data can run it without
re-deriving the physics.

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

CANDIDATE 1 -- THE EXECUTED CHECK (chem-cfl): build F_phase(T) = E_mag,phase(T) - T*S_mag,phase(T) +
E_lattice,phase for phase in {HT, LT}, reusing the already-validated N=2 Boltzmann trace
(``nb3x8_susceptibility.n2_spectrum``, the same primitive ``nb3x8_thermo`` uses) for E_mag/S_mag with
``NB3X8_HT_BULK`` vs ``NB3X8_LT_BULK``. Because E_lattice enters F_HT(T)-F_LT(T) linearly and only
once, the E_lattice that makes the crossing land EXACTLY on the measured Tc is solved directly (no
root-finding): ``E_lattice = Tc*dS_mag(Tc) - dE_mag(Tc)``. The physically-plausible range for
E_lattice is anchored to a SECOND primary reference this sandbox fetched this session -- Mortazavi,
Zhuang & Rabczuk, "A first-principles study on the physical properties of two-dimensional Nb3Cl8,
Nb3Br8 and Nb3I8," arXiv:2412.03360v1 (DFT exfoliation energies and monolayer lattice constants,
verified against the primary PDF text directly, not a search-engine summary, per the chem-jiy
lesson) -- using the full monolayer exfoliation energy per formula unit as a loose UPPER bound (a
partial HT<->LT structural relaxation within the bound interlayer well cannot plausibly cost more
than separating the layers to infinity).

FINDING: the needed E_lattice is 9.94 meV (Cl, 1.7% of the 596.9 meV exfoliation bound) / 41.6 meV
(Br, 5.6% of the 740.1 meV bound) -- both small, positive (correct-sign: LT net favored at T=0), and
far from "extreme tuning" against the bound. The candidate SURVIVES this check: a modest, physically
plausible lattice energy (a small fraction of the family's own measured interlayer binding scale)
reproduces the measured Tc without extreme tuning. This is NOT confirmation that a cooperative
structural transition IS the mechanism -- the bound is generous (full exfoliation, not a direct HT/LT
structural-energy difference the source paper does not resolve) and the magnetic parameters
(NB3X8_HT_BULK) come from a single-dimer electronic model with no feedback from the lattice -- it
only means this candidate is not ruled out by the data now available, unlike candidate 3.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from nb3x8_gaps import NB3X8_CLUSTERS, NB3X8_HT_BULK, NB3X8_LT_BULK
from nb3x8_magnetometry import EXPERIMENT, MEV_PER_K, predicted_transition_K
from nb3x8_susceptibility import n2_spectrum
from nb3x8_thermo import entropy

# Same dimerized (LT) structural motif, two stacking variants, per halide -- both already in
# NB3X8_CLUSTERS (arXiv:2501.10320 Table I). Used ONLY as an empirical anchor for "how much can
# t_s_perp plausibly move without a phase change"; no new parameters.
_BULK_BIL_KEYS = {
    "Nb3F8": ("F  LT-bulk", "F  LT-bil"),
    "Nb3Cl8": ("Cl LT-bulk", "Cl LT-bil"),
    "Nb3Br8": ("Br LT-bulk", "Br LT-bil"),
    "Nb3I8": ("I  LT-bulk", "I  LT-bil"),
}

# chem-cfl sourced elastic/lattice data for candidate 1 (below) and it is now tested -- see
# `candidate1_verdict`. Candidate 2 remains open, with the pre-registered check that WOULD kill or
# support it once an in-plane exchange constant exists, and why it cannot be run now.
UNTESTED_CANDIDATES = {
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
            "No in-plane kagome exchange constant is tabulated anywhere in this repo, in "
            "arXiv:2501.10320's vendored tables (Table I/IV cover the interlayer dimer only), or in "
            "arXiv:2412.03360 (the elastic/lattice reference chem-cfl sourced for candidate 1 -- it "
            "reports exfoliation energies and elastic constants, not a magnetic exchange constant). "
            "WebFetch (unlike WebSearch) DOES work in this sandbox (chem-cfl, this session) -- a "
            "future session should retry sourcing with WebFetch on direct URLs/the arXiv API before "
            "assuming this is unreachable, rather than relying on chem-phh/a04's WebSearch-only "
            "finding."
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


# ---------------------------------------------------------------------------------------------
# CANDIDATE 1 (chem-cfl): cooperative/first-order structural-transition free-energy crossing check
# ---------------------------------------------------------------------------------------------

# Mortazavi, Zhuang & Rabczuk, "A first-principles study on the physical properties of
# two-dimensional Nb3Cl8, Nb3Br8 and Nb3I8," arXiv:2412.03360v1 -- DFT (VASP, PBE-GGA+D3)
# exfoliation energies and monolayer hexagonal lattice constants, read directly off the primary PDF
# text this session (not a search-engine summary, per the chem-jiy lesson -- see
# specs/SPEC_nb3x8_tc_attribution.md S3b). Independent of arXiv:2501.10320 (used by nb3x8_gaps.py):
# this paper reports no electronic/magnetic parameters, only structural/elastic ones, and it
# resolves only the single (LT, dimerized) ground-state structure -- no HT phase.
NB3X8_EXFOLIATION = {
    "Nb3Cl8": dict(gamma_J_per_m2=0.24, a_monolayer_A=6.783),
    "Nb3Br8": dict(gamma_J_per_m2=0.27, a_monolayer_A=7.121),
}
# 1 J/m^2 acting over 1 Angstrom^2 of area, in meV: 1e-20 J / 1.602176634e-19 J/eV * 1000 meV/eV.
_J_PER_M2_TIMES_A2_TO_MEV = 1e-20 / 1.602176634e-19 * 1000.0


def mean_energy(U0: float, t: float, Us: float, T):
    """Exact reduced magnetic mean energy <E>(T) over the N=2 sector (k_B=1, own-ground referenced,
    <E>(0)=0) -- the same Boltzmann trace ``nb3x8_thermo.entropy``/``heat_capacity`` use internally,
    exposed here because the F_phase(T) check below needs E_mag(T) directly."""
    e, _ = n2_spectrum(U0, t, Us)
    T = np.atleast_1d(np.asarray(T, dtype=float))
    b = np.exp(-(e[:, None] - e.min()) / T)
    Z = b.sum(0)
    e_avg = (e[:, None] * b).sum(0) / Z
    return float(e_avg[0]) if e_avg.size == 1 else e_avg


def exfoliation_energy_per_fu(name: str) -> float:
    """Upper bound (meV per Nb3X8 formula unit) on a physically plausible lattice/structural energy
    scale: the full monolayer exfoliation energy from arXiv:2412.03360, assuming the standard
    1-formula-unit hexagonal monolayer cell (area = (sqrt(3)/2)*a^2, a = the paper's reported
    monolayer lattice constant). A HT<->LT structural relaxation is a PARTIAL distortion within the
    bound interlayer well, so it cannot plausibly cost more than separating the layers to infinity --
    a loose, defensible upper bound, not a measured HT/LT energy difference. Robust to the
    1-vs-2-formula-units-per-cell convention: halving this bound (298/370 meV) still leaves it far
    above the needed E_lattice computed below (9.9/41.6 meV)."""
    d = NB3X8_EXFOLIATION[name]
    area = (np.sqrt(3.0) / 2.0) * d["a_monolayer_A"] ** 2
    return d["gamma_J_per_m2"] * area * _J_PER_M2_TIMES_A2_TO_MEV


def required_e_lattice(name: str) -> float:
    """E_lattice (meV per dimer, LT minus HT) such that F_HT(Tc) = F_LT(Tc) lands exactly on the
    MEASURED Tc (``EXPERIMENT[name]["Tc_K"]``), using NB3X8_HT_BULK vs NB3X8_LT_BULK for the magnetic
    (own-ground-referenced) energy/entropy. F_phase(T) = E_mag,phase(T) - T*S_mag,phase(T) +
    E_lattice,phase (LT carries -E_lattice relative to HT, i.e. LT is lattice-stabilized). Since
    E_lattice enters F_HT(T)-F_LT(T) linearly and only once, solving for the value that zeroes it at
    T=Tc_measured is direct algebra, not root-finding: E_lattice = Tc*dS_mag(Tc) - dE_mag(Tc), with
    dE_mag = E_mag,HT - E_mag,LT, dS_mag = S_mag,HT - S_mag,LT (both at T=Tc_measured)."""
    Tc_meV = EXPERIMENT[name]["Tc_K"] * MEV_PER_K
    ht, lt = NB3X8_HT_BULK[name], NB3X8_LT_BULK[name]
    dE = mean_energy(**ht, T=Tc_meV) - mean_energy(**lt, T=Tc_meV)
    dS = entropy(**ht, T=Tc_meV) - entropy(**lt, T=Tc_meV)
    return float(Tc_meV * dS - dE)


def candidate1_verdict(name: str) -> dict:
    """Candidate 1 ('cooperative/first-order structural transition') kill test, chem-cfl.

    KILL bar: the E_lattice needed to pull the HT/LT free-energy crossing down onto the measured Tc
    must be nonpositive (wrong sign -- would mean HT, not LT, needs to be lattice-favored at T=0, the
    opposite of the observed LT ground state) or must exceed the physically plausible upper bound
    (``exfoliation_energy_per_fu``) -- a bound generous enough that only a candidate needing an
    unphysically large lattice energy fails it. SURVIVES otherwise, with ``needed_frac_of_bound``
    reported so "without extreme tuning" is an auditable number, not an adjective.
    """
    needed = required_e_lattice(name)
    bound = exfoliation_energy_per_fu(name)
    return dict(
        name=name,
        needed_meV=needed,
        bound_meV=bound,
        needed_frac_of_bound=needed / bound,
        killed=not (0.0 < needed < bound),
    )


if __name__ == "__main__":
    print("chem-phh/chem-cfl: Nb3X8 Tc-overprediction attribution -- three candidates, two tested\n")

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

    print("\nCandidate 1 (cooperative/first-order structural transition, F_phase(T) crossing):")
    for name in ("Nb3Cl8", "Nb3Br8"):
        v = candidate1_verdict(name)
        print(f"  {name}: needs E_lattice = {v['needed_meV']:.2f} meV to match measured Tc "
              f"({100*v['needed_frac_of_bound']:.1f}% of the {v['bound_meV']:.1f} meV exfoliation-"
              f"energy bound, arXiv:2412.03360) -> {'KILLED' if v['killed'] else 'survives'}")

    print("\nCandidate 2: BLOCKED on data this sandbox has not sourced; see UNTESTED_CANDIDATES for "
          "the pre-registered check it would need to pass.")
    for key, v in UNTESTED_CANDIDATES.items():
        print(f"\n  [{key}] (follow-up: {v['follow_up_bead']})\n    blocked_on: {v['blocked_on']}")
