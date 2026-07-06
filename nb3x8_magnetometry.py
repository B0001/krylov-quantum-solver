#!/usr/bin/env python3
"""
Nb3X8 interlayer exchange vs measured magnetometry -- a parameter-free prediction against the lab.

The response trilogy (odmd_spin) extracted the interlayer singlet-triplet gap J of the Nb3X8
bilayer dimers from the paper's cRPA parameters (J = 66.2 / 119.1 meV for Cl / Br). The finite-T
study (nb3x8_susceptibility) then acknowledged -- but explicitly scoped OUT -- the real materials'
magnetic transitions (Sheckelton et al. 2017: the ~90 K Nb3Cl8 interlayer-singlet transition;
Haraguchi et al. 2017; Nb3Br8 ~382 K). This module closes that loop: it turns the acknowledged
experimental contact into a FALSIFIABLE comparison, with no fit parameters.

The physical picture is exact for this test: the low-T nonmagnetic state of Nb3Cl8/Nb3Br8 IS the
interlayer dimerization of adjacent-layer Nb3 clusters into singlets, so the very coupling J that
nb3x8_susceptibility/odmd_spin compute is the one that sets the singlet-formation temperature. A
Heisenberg (Bleaney-Bowers) dimer's susceptibility peaks at k_B T_max ~ 0.625 J; that T_max is the
ab-initio prediction of the transition SCALE.

THE FINDING (specs/SPEC_nb3x8_magnetometry.md): the isolated-bilayer J gets the chemistry right and
the number wrong, in a specific, quantified way --
  * ORDERING / SCALE (the win): T_max(Cl) < T_max(Br) matches the observed 90 K < 382 K, and each
    prediction is within an order of magnitude of experiment.
  * THE MISS (the cluster->solid gap): the isolated dimer OVERPREDICTS the measured transition by
    5.3x (Cl) and 2.3x (Br) -- an overcoupling that WEAKENS monotonically down the halide series,
    exactly the isolated-cluster -> cooperative-lattice renormalization (a single dimer over-counts
    the coupling; the real transition is a cooperative lattice effect that sets in below the
    single-dimer chi maximum).
  * TWO DISTINCT COUPLINGS (the second finding): the interlayer J does NOT set the high-T
    Curie-Weiss theta. -J/4 = -192 K for Cl, 15x the measured theta_W = -13.1 K. The measured
    theta_W is governed by the WEAK IN-PLANE exchange this isolated bilayer dimer does not contain
    -- interlayer J sets T_c, in-plane exchange sets theta_W.

HONEST SCOPE: isolated bilayer dimer, no in-plane kagome exchange, no phonons, no cooperative /
first-order structural transition, density-density interactions only. The model sets *scales*, it
does not reproduce a first-order transition. Nb3I8 has no interlayer-singlet transition (a
different, moment-retaining ground state) and is excluded from the experimental comparison. The
experimental targets are quoted with primary-source citations; this is a comparison, not a fit.
"""
from __future__ import annotations

from scipy.optimize import minimize_scalar

from nb3x8_gaps import NB3X8_LT_BULK
from nb3x8_susceptibility import curie_weiss_theta, susceptibility
from odmd_spin import dimer_exchange_analytic

MEV_PER_K = 0.0861733  # meV per kelvin (k_B); K = meV / MEV_PER_K

# Measured magnetometry of the real solids, with primary sources. Only the members with an
# established interlayer-singlet transition are listed (Nb3I8 has none -- see docstring).
#   Tc_K   : temperature of the magnetic -> nonmagnetic (interlayer-singlet) transition [K]
#   theta_K: high-T Curie-Weiss Weiss temperature [K] (None where not quoted here)
EXPERIMENT = {
    "Nb3Cl8": dict(Tc_K=90.0, theta_K=-13.1,
                   ref="Sheckelton et al., Inorg. Chem. Front. 4, 481 (2017); arXiv:1701.05528"),
    "Nb3Br8": dict(Tc_K=382.0, theta_K=None,
                   ref="Haraguchi et al., Inorg. Chem. 56, 3483 (2017)"),
}


def chi_max_temperature(U0: float, t: float, Us: float) -> float:
    """Temperature (meV) of the exact dimer susceptibility maximum -- the ab-initio prediction of
    the singlet-formation scale. Maximizes the exact N=2 chi(T); ~0.625 J for a clean spin gap."""
    J = dimer_exchange_analytic(U0, t, Us)
    lo, hi = 0.05 * max(J, 1e-6), 3.0 * max(J, 1e-6)
    res = minimize_scalar(lambda T: -susceptibility(U0, t, Us, T),
                          bounds=(lo, hi), method="bounded")
    return float(res.x)


def predicted_transition_K(U0: float, t: float, Us: float) -> float:
    """The chi-maximum temperature in kelvin -- the predicted transition scale of the solid."""
    return chi_max_temperature(U0, t, Us) / MEV_PER_K


def overprediction_factor(name: str) -> float:
    """Predicted transition temperature / measured Tc (dimensionless). > 1 => overprediction."""
    p = NB3X8_LT_BULK[name]
    return predicted_transition_K(**p) / EXPERIMENT[name]["Tc_K"]


def theta_over_measured(name: str) -> float:
    """(-J/4 in K) / measured theta_W -- how far the interlayer J is from the theta_W coupling."""
    J = dimer_exchange_analytic(**NB3X8_LT_BULK[name])
    theta_pred_K = curie_weiss_theta(J) / MEV_PER_K
    return theta_pred_K / EXPERIMENT[name]["theta_K"]


if __name__ == "__main__":
    print("Nb3X8 interlayer exchange vs measured magnetometry (isolated bilayer dimer, no fit)")
    print(f"{'material':>8} | {'J(meV)':>8} | {'Tmax(meV)':>9} | {'pred Tc(K)':>10} | "
          f"{'obs Tc(K)':>9} | {'overpred':>8} | {'-J/4(K)':>8} | {'obs theta':>9}")
    for name, p in NB3X8_LT_BULK.items():
        J = dimer_exchange_analytic(**p)
        Tmax = chi_max_temperature(**p)
        predK = predicted_transition_K(**p)
        thetaK = curie_weiss_theta(J) / MEV_PER_K
        e = EXPERIMENT.get(name)
        obs_tc = f"{e['Tc_K']:.0f}" if e else "--"
        over = f"{overprediction_factor(name):.2f}x" if e else "--"
        obs_th = f"{e['theta_K']}" if e and e["theta_K"] is not None else "--"
        print(f"{name:>8} | {J:8.2f} | {Tmax:9.2f} | {predK:10.0f} | {obs_tc:>9} | "
              f"{over:>8} | {thetaK:8.0f} | {obs_th:>9}")
    print("\nFinding: isolated-dimer J reproduces the Cl<Br ordering and the scale, but overpredicts")
    print("Tc by 5.3x (Cl) / 2.3x (Br) -- a cluster->solid renormalization that weakens down the")
    print("series; and -J/4 (=-192 K, Cl) != measured theta_W (-13.1 K): interlayer J sets Tc, weak")
    print("in-plane exchange sets theta_W. See specs/SPEC_nb3x8_magnetometry.md.")
