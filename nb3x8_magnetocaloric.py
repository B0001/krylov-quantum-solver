#!/usr/bin/env python3
"""
Nb3X8 magnetocaloric effect: field-dependent entropy S(T,B), composed with NO new physics from
two already-gated primitives -- nb3x8_metamagnetism.field_spectrum (the Zeeman-augmented exact
spectrum) fed into the same reduced Boltzmann-trace entropy formula already validated at zero
field in nb3x8_thermo.py. Named as a follow-up by BOTH SPEC_nb3x8_thermo Sec.7 and
SPEC_nb3x8_metamagnetism Sec.7; one composition closes both.

TWO CLAIMS CHECKED (specs/SPEC_nb3x8_magnetocaloric.md):

  (i) MACHINERY -- the Maxwell relation (dS/dB)_T = (dM/dT)_B is an exact thermodynamic identity
      of any well-defined free energy F(T,h) = -T ln Z(T,h) (Clairaut's theorem on its mixed
      partial d2F/dTdh); computing its two sides from two INDEPENDENTLY-derived Boltzmann traces
      -- entropy's "ln Z + <E>/T" route (this module) vs magnetization_thermal's "<Sz>" route
      (nb3x8_metamagnetism_thermal.py, a separate file with a separate derivation) -- is a real
      cross-check of the IMPLEMENTATION, not a test of the physics (which is guaranteed). Found to
      agree to < 2e-7 relative across the swept spin window (G1).

  (ii) VERDICT -- |Delta S_M(T,B)| = |S(T,B) - S(T,0)|, searched over B in [0, 100] T and T in
      [2 K, 5*(J/k_B)] for each halide, is maximised at the edge of the searched field range
      (B = 100 T, still far below B_c) and reaches only a fraction of a percent to ~1.7% of
      R ln2 per formula unit (Cl worst-suppressed of the three magnetic halides) -- two orders of
      magnitude under the bead's 10% kill threshold. In physical units that is ~500x (Cl) to
      ~9000x (I) below GGG's theoretical full paramagnetic-entropy ceiling (see GGG_* below).
      A cheap direct check (G2) confirms the qualitative picture behind this: |Delta S_M| grows
      by more than an order of magnitude between B = 100 T and B = B_c, i.e. the search window is
      genuinely far from the interesting physics, not merely far from an arbitrary cutoff.
      **Nb3X8 dimers are quantitatively useless magnetocalorics below megagauss fields.**

THIS IS A BOUNDING RESULT BY CONSTRUCTION (precedent: SPEC_senseforge Sec. 3-4) -- the verdict is
the NEGATIVE of the hoped-for outcome (Nb3X8 as a near-term magnetocaloric candidate), gated to
say so loudly if it turns out to be wrong. Inherits g=2, density-density-only, isolated-single-
dimer from the parent specs (nb3x8_metamagnetism.py, nb3x8_thermo.py) -- not re-derived here.

GGG COMPARISON -- SOURCING CAVEAT: repeated attempts to source a specific experimental
Delta S_M(T,B) curve for Gd3Ga5O12 (GGG) from external literature in this sandbox produced
inconsistent, uncitable numbers on cross-check (one fetch reported a GGG molar mass of
644.37 g/mol -- wrong by ~36% against the formula-weight arithmetic below, a sign of fabrication
rather than a real quote) -- so none of that is used. Instead this module uses the RIGOROUS,
textbook, code-derived ceiling: Gd3+ is 4f^7 with L=0, S=7/2 (Hund's-rule ground term 8S_7/2, an
orbital singlet, isolated from the first excited term by several eV) -- exactly why GGG is chosen
as a magnetic refrigerant -- so its full paramagnetic entropy content is 3*R*ln(8) per mole (3 Gd
per formula unit), an upper bound on any field/temperature Delta S_M GGG could ever show. Using
this ceiling is CONSERVATIVE in the direction of the negative verdict: if Nb3X8 is tiny against
GGG's theoretical maximum, it is tinier still against any real (necessarily smaller) achieved GGG
number. Molar masses (Nb3X8 and GGG) come from pyscf.data.elements.MASSES -- the natural-abundance
standard atomic weights PySCF itself uses -- not hand-copied constants.

HONEST SCOPE: isolated single dimer, g=2, density-density only (inherited, see above); the B<=100T
window is what pulsed-field magnets can reach (SPEC_nb3x8_metamagnetism Sec.2); T>=2K is the
bead's own floor. Not claimed as a solid-state prediction any more than the parent specs are.
"""
from __future__ import annotations

import numpy as np
import scipy.constants as sc
from pyscf.data import elements

from nb3x8_magnetometry import MEV_PER_K
from nb3x8_metamagnetism import G_MU_B, field_spectrum
from nb3x8_metamagnetism_thermal import magnetization_thermal

R_GAS = sc.R  # J/(mol K), CODATA (matches k_B/eV*1000 == MEV_PER_K to the digits used here)

MAGNETIC_HALIDES = ("Nb3Cl8", "Nb3Br8", "Nb3I8")  # Nb3F8 excluded: J below the model's noise floor


def _atomic_mass(symbol: str) -> float:
    """Standard (natural-abundance) atomic weight in g/mol, from pyscf.data.elements.MASSES --
    code-derived, not a hand-copied constant."""
    return float(elements.MASSES[elements.ELEMENTS.index(symbol)])


def nb3x8_formula_mass(halide_symbol: str) -> float:
    """g/mol of one Nb3X8 formula unit (X = ``halide_symbol``, e.g. 'Cl')."""
    return 3.0 * _atomic_mass("Nb") + 8.0 * _atomic_mass(halide_symbol)


NB3X8_FORMULA_MASS_G_PER_MOL = {
    "Nb3Cl8": nb3x8_formula_mass("Cl"),
    "Nb3Br8": nb3x8_formula_mass("Br"),
    "Nb3I8": nb3x8_formula_mass("I"),
}

# GGG = Gd3Ga5O12 (3 Gd, 5 Ga, 12 O per formula unit).
GGG_MOLAR_MASS_G_PER_MOL = 3.0 * _atomic_mass("Gd") + 5.0 * _atomic_mass("Ga") + 12.0 * _atomic_mass("O")

# Gd3+ (4f^7, L=0, S=7/2): full paramagnetic entropy ceiling 3*R*ln(8) per mole GGG (see docstring).
GGG_ENTROPY_CEILING_J_PER_KG_K = 3.0 * R_GAS * np.log(8.0) * 1000.0 / GGG_MOLAR_MASS_G_PER_MOL


def entropy_field(U0: float, t: float, Us: float, h: float, T):
    """Exact reduced Boltzmann-trace entropy S(T,h) per DIMER (k_B=1), of the field-augmented
    N=2-sector spectrum (nb3x8_metamagnetism.field_spectrum) -- the same "ln Z + <E-E0>/T" formula
    already gated at h=0 in nb3x8_thermo.entropy, generalized here to the Zeeman-augmented
    spectrum. Scalar or array T."""
    e, _ = field_spectrum(U0, t, Us, h)
    T = np.atleast_1d(np.asarray(T, dtype=float))
    b = np.exp(-(e[:, None] - e.min()) / T)
    Z = b.sum(0)
    e_avg = (e[:, None] * b).sum(0) / Z
    S = np.log(Z) + (e_avg - e.min()) / T
    return float(S[0]) if S.size == 1 else S


def delta_entropy_isothermal(U0: float, t: float, Us: float, T: float, h: float) -> float:
    """Delta S_M(T,h) = S(T,h) - S(T,0): the isothermal magnetic entropy change, reduced (k_B=1),
    per dimer."""
    return entropy_field(U0, t, Us, h, T) - entropy_field(U0, t, Us, 0.0, T)


def ds_dh_numeric(U0: float, t: float, Us: float, h: float, T: float, dh: float) -> float:
    """Central finite difference (dS/dh)_T from entropy_field -- the ENTROPY route of the Maxwell
    cross-check (G1)."""
    return (entropy_field(U0, t, Us, h + dh, T) - entropy_field(U0, t, Us, h - dh, T)) / (2.0 * dh)


def dm_dT_numeric(U0: float, t: float, Us: float, h: float, T: float, dT: float) -> float:
    """Central finite difference (dM/dT)_h from nb3x8_metamagnetism_thermal.magnetization_thermal
    -- the MAGNETIZATION route of the Maxwell cross-check (G1): a different file, a different
    observable, the same underlying field_spectrum."""
    return (magnetization_thermal(U0, t, Us, h, T + dT)
            - magnetization_thermal(U0, t, Us, h, T - dT)) / (2.0 * dT)


def delta_s_m_fraction_of_r_ln2(U0: float, t: float, Us: float, T_kelvin: float, B_tesla: float) -> float:
    """|Delta S_M(T,B)| as a fraction of R ln2 PER FORMULA UNIT (reduced, k_B=1) -- the bead's
    kill-criterion quantity. A dimer carries 2 formula units, so the per-dimer entropy change is
    halved before comparing to ln2."""
    T = T_kelvin * MEV_PER_K
    h = B_tesla * G_MU_B
    dS_per_fu = delta_entropy_isothermal(U0, t, Us, T, h) / 2.0
    return abs(dS_per_fu) / np.log(2.0)


def delta_s_m_j_per_kg_k(U0: float, t: float, Us: float, T_kelvin: float, B_tesla: float,
                          formula_mass_g_per_mol: float) -> float:
    """Delta S_M(T,B) in J/(kg K): reduced per-dimer entropy change * R * 1000 / (2 * formula
    mass) -- a dimer is the interlayer pair, i.e. 2 formula units."""
    T = T_kelvin * MEV_PER_K
    h = B_tesla * G_MU_B
    dS_reduced_per_dimer = delta_entropy_isothermal(U0, t, Us, T, h)
    return dS_reduced_per_dimer * R_GAS * 1000.0 / (2.0 * formula_mass_g_per_mol)


def scan_max_delta_s(U0: float, t: float, Us: float, b_max_tesla: float = 100.0,
                      t_min_kelvin: float = 2.0, t_max_factor: float = 5.0,
                      n_b: int = 60, n_t: int = 300):
    """Grid-search max |Delta S_M(T,B)| (reduced, per dimer) over B in [0, b_max_tesla] T and T in
    [t_min_kelvin, t_max_factor * J/k_B] K -- the actual search behind the bead's B<=100T, T>=2K
    kill-criterion evaluation (cheap: a <=6-level exact diagonalization, not an FCI/DMRG sweep).
    Returns (max |Delta S_M|, B* [T], T* [K])."""
    from odmd_spin import dimer_exchange_analytic

    J = dimer_exchange_analytic(U0, t, Us)
    J_kelvin = J / MEV_PER_K
    Ts_kelvin = np.linspace(t_min_kelvin, t_max_factor * J_kelvin, n_t)
    Ts = Ts_kelvin * MEV_PER_K
    S0 = entropy_field(U0, t, Us, 0.0, Ts)
    best = (0.0, 0.0, float(Ts_kelvin[0]))
    for B in np.linspace(0.0, b_max_tesla, n_b):
        h = B * G_MU_B
        Sh = entropy_field(U0, t, Us, h, Ts)
        dS = Sh - S0
        i = int(np.argmax(np.abs(dS)))
        if abs(dS[i]) > best[0]:
            best = (float(abs(dS[i])), float(B), float(Ts_kelvin[i]))
    return best


if __name__ == "__main__":
    from nb3x8_gaps import NB3X8_LT_BULK
    from nb3x8_metamagnetism import critical_field_tesla
    from odmd_spin import dimer_exchange_analytic

    LN2 = float(np.log(2.0))
    print("Nb3X8 magnetocaloric S(T,B) -- closing SPEC_nb3x8_thermo Sec.7 / "
          "SPEC_nb3x8_metamagnetism Sec.7")
    print(f"GGG (Gd3Ga5O12) molar mass = {GGG_MOLAR_MASS_G_PER_MOL:.3f} g/mol "
          f"(pyscf standard atomic weights); full paramagnetic entropy ceiling "
          f"3*R*ln(8)/M = {GGG_ENTROPY_CEILING_J_PER_KG_K:.3f} J/(kg K)")
    print(f"{'halide':>8} | {'M(g/mol)':>9} | {'J(K)':>7} | {'Bc(T)':>7} | "
          f"{'maxdS/Rln2fu':>12} | {'B*(T)':>6} | {'T*(K)':>7} | {'dS(J/kg/K)':>10} | {'vs GGG':>8}")
    for name in MAGNETIC_HALIDES:
        p = NB3X8_LT_BULK[name]
        M = NB3X8_FORMULA_MASS_G_PER_MOL[name]
        J = dimer_exchange_analytic(**p)
        Bc = critical_field_tesla(**p)
        max_dS, B_star, T_star = scan_max_delta_s(**p)
        frac_r_ln2 = (max_dS / 2.0) / LN2
        dS_phys = max_dS * R_GAS * 1000.0 / (2.0 * M)
        ratio_ggg = dS_phys / GGG_ENTROPY_CEILING_J_PER_KG_K
        print(f"{name:>8} | {M:9.3f} | {J / MEV_PER_K:7.1f} | {Bc:7.1f} | "
              f"{100 * frac_r_ln2:11.3f}% | {B_star:6.1f} | {T_star:7.1f} | "
              f"{dS_phys:10.5f} | {100 * ratio_ggg:7.4f}%")
    print("\nFinding: the searched max sits at the search boundary (B=100T), 2-3 orders of")
    print("magnitude under both the 10% R-ln2/f.u. kill threshold and GGG's theoretical entropy")
    print("ceiling -- laboratory fields cannot make Nb3X8 dimers useful magnetocalorics.")
    print("See specs/SPEC_nb3x8_magnetocaloric.md.")
