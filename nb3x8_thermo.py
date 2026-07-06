#!/usr/bin/env python3
"""
Magnetic heat capacity C_m(T) and entropy S_m(T) of the Nb3X8 downfolded dimers -- completing the
thermodynamic triad (nb3x8_susceptibility gave chi(T); this adds C and S).

Exact N=2 Boltzmann trace over the same extended-Hubbard dimer (cRPA parameters of arXiv:2501.10320,
via nb3x8_gaps), reduced units k_B = 1, energies & T in meV, entropy per dimer in units of k_B:

    C_m(T) = (<E^2> - <E>^2) / T^2          S_m(T) = ln Z + <E - E_0>/T

REPRODUCTION of standard magnetothermodynamics (the two-level Schottky anomaly; the R ln(2S+1)
magnetic-entropy plateau -- Gopal, *Specific Heats at Low Temperatures*, 1966) applied to new
ab-initio parameters. The Schottky/entropy fingerprints of localized moments are textbook; the
family-wide ab-initio numbers and the charge-scale plateau boundary are the contribution. Prior
Nb3X8 heat-capacity data exists experimentally (Sheckelton 2017 -- the ~90 K Nb3Cl8 transition).

THE FINDING (specs/SPEC_nb3x8_thermo.md):
  * The magnetic Schottky peak sits at T ~ 0.352 J for the WHOLE family (0.351-0.357), pinned by
    the analytic two-level (singlet g0=1, triplet g1=3) result -- a J-scale fingerprint distinct
    from chi(T)'s peak at 0.625 J. Cross-tie: C-peak / chi-peak = 0.3515/0.625 ~ 0.56 is universal
    and material-independent (both are exact two-level features).
  * The localized-moment entropy plateau S = R ln 4 per dimer (= R ln 2 per S=1/2 cluster) is a
    clean plateau ONLY when the charge scale E_s >> J: its flatness (min |dS/dlnT| between the spin
    and charge features) and its deviation from ln 4 both grow monotonically as E_s/J shrinks --
    sharp for Nb3Cl8 (E_s/J ~ 17, dev -1.3%), gone for Nb3I8 (E_s/J ~ 3, dev +17%, S already
    climbing to ln 6). This is the SAME charge-scale boundary that bounds the Bleaney-Bowers regime
    in nb3x8_susceptibility, now read off the entropy.

HONEST SCOPE: isolated single dimer -- no inter-dimer coupling, no phonon/lattice heat capacity, no
structural transition (real Nb3Cl8 has a first-order one at ~90 K); density-density interactions
only; magnetic (spin+charge-sector) contribution only. Nb3F8 (J ~ 0.05 meV, below the model's own
neglected terms) has no resolvable spin feature and is excluded. A reference table, not a solid-state
heat-capacity prediction.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar

from nb3x8_susceptibility import ionic_singlet_energy, n2_spectrum
from odmd_spin import dimer_exchange_analytic


def heat_capacity(U0: float, t: float, Us: float, T):
    """Exact reduced magnetic heat capacity C_m(T) = Var(E)/T^2 over the N=2 sector (k_B=1)."""
    e, _ = n2_spectrum(U0, t, Us)
    T = np.atleast_1d(np.asarray(T, dtype=float))
    b = np.exp(-(e[:, None] - e.min()) / T)
    Z = b.sum(0)
    e_avg = (e[:, None] * b).sum(0) / Z
    e2_avg = ((e[:, None] ** 2) * b).sum(0) / Z
    C = (e2_avg - e_avg ** 2) / T ** 2
    return float(C[0]) if C.size == 1 else C


def entropy(U0: float, t: float, Us: float, T):
    """Exact reduced magnetic entropy S_m(T) = ln Z + <E-E_0>/T (k_B=1, per dimer). S(0)=0 (singlet
    ground), S(sector -> inf) = ln 6 (the full N=2 manifold)."""
    e, _ = n2_spectrum(U0, t, Us)
    T = np.atleast_1d(np.asarray(T, dtype=float))
    b = np.exp(-(e[:, None] - e.min()) / T)
    Z = b.sum(0)
    e_avg = (e[:, None] * b).sum(0) / Z
    S = np.log(Z) + (e_avg - e.min()) / T
    return float(S[0]) if S.size == 1 else S


def two_level_schottky_ratio(g_ratio: float = 3.0) -> float:
    """Analytic Schottky peak T_pk/Delta for a two-level system with degeneracies (g0, g1=g_ratio*g0)
    and gap Delta. For the singlet(1)/triplet(3) dimer, g_ratio=3 -> ~0.3515."""
    def neg_C(x):                                          # x = T/Delta
        f = g_ratio * np.exp(-1.0 / x)
        return -(1.0 / x ** 2) * f / (1.0 + f) ** 2
    return float(minimize_scalar(neg_C, bounds=(0.05, 2.0), method="bounded").x)


def schottky_peak_temperature(U0: float, t: float, Us: float) -> float:
    """Temperature (meV) of the magnetic (spin) Schottky maximum of C_m(T) -- searched in the spin
    window below the charge scale, so the ionic peak near E_s is excluded."""
    J = dimer_exchange_analytic(U0, t, Us)
    res = minimize_scalar(lambda T: -heat_capacity(U0, t, Us, T),
                          bounds=(0.1 * J, 1.0 * J), method="bounded")
    return float(res.x)


def entropy_plateau(U0: float, t: float, Us: float):
    """(S_plateau, flatness) of the localized-moment entropy plateau between the spin gap and the
    charge scale. ``flatness`` = min |dS/d ln T| on [2J, E_s/3] (0 = perfect R ln 4 plateau);
    ``S_plateau`` is S there. A clean plateau (flat, ~ln 4) needs E_s >> J."""
    J = dimer_exchange_analytic(U0, t, Us)
    E_s = ionic_singlet_energy(U0, t, Us)
    ln_t = np.linspace(np.log(2.0 * J), np.log(E_s / 3.0), 4000)
    S = entropy(U0, t, Us, np.exp(ln_t))
    dS = np.gradient(S, ln_t)
    i = int(np.argmin(np.abs(dS)))
    return float(S[i]), float(np.abs(dS[i]))


if __name__ == "__main__":
    from nb3x8_gaps import NB3X8_LT_BULK
    from nb3x8_magnetometry import chi_max_temperature

    LN4 = np.log(4.0)
    r_analytic = two_level_schottky_ratio(3.0)
    print("Nb3X8 magnetic heat capacity & entropy -- exact N=2 trace (isolated dimer, reproduction)")
    print(f"analytic two-level (1,3) Schottky peak: T_pk/J = {r_analytic:.4f};  "
          f"C-peak/chi-peak = {r_analytic / 0.625:.3f} (universal)")
    print(f"{'set':>7} | {'J(meV)':>8} | {'E_s/J':>6} | {'C_pk(meV)':>9} | {'T_pk/J':>7} | "
          f"{'T_pk(K)':>8} | {'S_plat/ln4':>10} | {'flatness':>8} | {'C/chi peak':>10}")
    for name, p in NB3X8_LT_BULK.items():
        J = dimer_exchange_analytic(**p)
        if J < 1.0:
            print(f"{name:>7} | {J:8.3f} |   ---  (J~0, no resolvable spin feature)")
            continue
        E_s = ionic_singlet_energy(**p)
        Tpk = schottky_peak_temperature(**p)
        Splat, flat = entropy_plateau(**p)
        ratio = Tpk / chi_max_temperature(**p)
        print(f"{name:>7} | {J:8.2f} | {E_s / J:6.1f} | {Tpk:9.2f} | {Tpk / J:7.3f} | "
              f"{Tpk / 0.0861733:8.0f} | {Splat / LN4:10.3f} | {flat:8.3f} | {ratio:10.3f}")
    print("\nFinding: Schottky peak at ~0.352 J (universal); R ln4/dimer (R ln2/cluster) entropy")
    print("plateau is clean only for E_s >> J (Cl) and gone for the iodide (E_s/J~3) -- the same")
    print("charge-scale boundary as chi(T). See specs/SPEC_nb3x8_thermo.md.")
