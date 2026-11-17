#!/usr/bin/env python3
"""
bead chem-phh/chem-cfl/chem-b0y -- what explains the Nb3X8 5.3x (Cl) / 2.3x (Br) Tc
overprediction, now that the coordination/mean-field attribution is KILLED (chem-g78,
nb3x8_gaps.py / tests/test_nb3x8_gaps_spec.py::test_G8_*)?

``SPEC_nb3x8_magnetometry.md`` section 7 names three candidates the isolated-dimer model cannot see
and chem-g78 did not test:
  (1) the cooperative/first-order structural transition itself,
  (2) in-plane kagome exchange,
  (3) a direct lattice renormalization of t_s_perp.

This module pre-registers a killable check for each, and all three are now executed. (3) first
(chem-phh), from data already vendored in this repo (arXiv:2501.10320's Table I/IV cRPA parameters,
via the already-validated ``nb3x8_magnetometry``/``odmd_spin`` primitives): KILLED. (1) (chem-cfl),
with elastic/lattice data from a second primary reference (arXiv:2412.03360): SURVIVES. (2)
(chem-b0y), with an in-plane exchange constant from a third (arXiv:2508.17352): KILLED. The two later
references were fetched over the network from inside the sandbox; the "no working internet access"
finding chem-phh/chem-a04 recorded rested on ``WebSearch`` alone.

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

import numpy as np
from scipy.optimize import brentq

from nb3x8_gaps import NB3X8_CLUSTERS, NB3X8_HT_BULK, NB3X8_LT_BULK
from nb3x8_magnetometry import EXPERIMENT, MEV_PER_K, predicted_transition_K
from nb3x8_susceptibility import n2_spectrum, susceptibility
from nb3x8_thermo import entropy
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

# All three candidates chem-g78 left open are now tested (candidate3_verdict, candidate1_verdict,
# candidate2_verdict). Kept, empty, so the honest-boundary gate can assert nothing is silently
# claimed: a future candidate goes here with its pre-registered check until it is run.
UNTESTED_CANDIDATES: dict = {}


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




# ---------------------------------------------------------------------------------------------
# CANDIDATE 2 (chem-b0y): in-plane exchange, RPA closure vs the sourced arXiv:2508.17352 J1/J2
# ---------------------------------------------------------------------------------------------

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
    print("chem-phh/chem-cfl/chem-b0y: Nb3X8 Tc-overprediction attribution -- three candidates, all tested\n")

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

    print("\nCandidate 1 (cooperative/first-order structural transition, F_phase(T) crossing):")
    for name in ("Nb3Cl8", "Nb3Br8"):
        v = candidate1_verdict(name)
        print(f"  {name}: needs E_lattice = {v['needed_meV']:.2f} meV to match measured Tc "
              f"({100*v['needed_frac_of_bound']:.1f}% of the {v['bound_meV']:.1f} meV exfoliation-"
              f"energy bound, arXiv:2412.03360) -> {'KILLED' if v['killed'] else 'survives'}")
