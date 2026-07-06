"""
Acceptance gates G1-G4 for specs/SPEC_nb3x8_susceptibility.md (finite-T susceptibility study).

Test-first: ``nb3x8_susceptibility`` does not exist yet, so this file is RED until the spec is
implemented. Claim: the exact-spectrum molar susceptibility chi(T) of the downfolded Nb3X8 dimer
reproduces the analytic Bleaney-Bowers singlet-triplet law up to a temperature set by the CHARGE
scale E_s (the first ionic singlet), not by the exchange J -- and deviates first, at the lowest
reduced temperature, for the iodides (smallest spin/charge separation E_s/J ~ 3).

This is a REPRODUCTION of standard magnetochemistry (Bleaney-Bowers 1952; Van Vleck; the exact
Hubbard-dimer thermodynamics, Carrascal et al. arXiv:1502.05038) applied to new ab-initio
parameters -- prior Nb3X8 chi(T) exists (Sheckelton 2017, Haraguchi 2017, Grytsiuk/Rosner
arXiv:2305.04854). Isolated single dimer, density-density only. PySCF/qiskit, no block2;
`make gates` runs it in its own process.
"""
import numpy as np

from nb3x8_gaps import NB3X8_CLUSTERS
from nb3x8_susceptibility import (
    bleaney_bowers,
    bb_deviation_temperature,
    curie_weiss_theta,
    ionic_singlet_energy,
    n2_spectrum,
    susceptibility,
)
from odmd_spin import dimer_exchange_analytic

# E_s/J > 10: clean Curie-Weiss window; the two iodides (E_s/J ~ 3) are the special case.
IODIDES = ("I  LT-bulk", "I  LT-bil")


def test_G1_exact_spectrum_and_pure_singlet_ground():
    """The N=2 triplet sits at J; the ground state is a pure singlet (<Sz_tot^2>=0, NOT the
    0.759 staggered moment of odmd_spin); chi*T -> 1/3 as T -> inf."""
    for name, p in NB3X8_CLUSTERS.items():
        energies, sz2 = n2_spectrum(**p)
        J = dimer_exchange_analytic(**p)
        assert abs(energies[0]) < 1e-9, (name, energies[0])           # ground at 0
        assert sz2[0] < 1e-9, (name, sz2[0])                          # pure singlet
        triplet = energies[1:4]                                       # 3-fold at J
        assert np.allclose(triplet, J, atol=1e-6), (name, triplet, J)
        # T->inf limit: chi*T -> <Sz^2> averaged over the full N=2 manifold = 2/6 = 1/3 exactly
        # (2 of the 6 states carry Sz^2=1). Gate the exact invariant, plus the numeric approach.
        assert abs(np.mean(sz2) - 1.0 / 3.0) < 1e-9, (name, np.mean(sz2))
        chiT_hi = susceptibility(**p, T=1e8) * 1e8
        assert abs(chiT_hi - 1.0 / 3.0) < 1e-4, (name, chiT_hi)


def test_G2_reproduces_bleaney_bowers_in_spin_regime():
    """At T = E_s/20 (below the charge scale) the exact chi matches Bleaney-Bowers < 1e-3. The
    chi*T Curie constant approaches the coupled-pair value 0.5 only in the J/T -> 0 idealization
    (at finite T/J, BB itself gives 2/(3+e^{J/T}) < 0.5), so it is gated to 0.5 only for the
    members with a clean deep-spin window (J/T <= 0.01 and T <= E_s/20 both reachable) -- a
    revision made during implementation, since the strict '=0.5' over-specified the finite-T
    value (see spec G2)."""
    for name, p in NB3X8_CLUSTERS.items():
        J = dimer_exchange_analytic(**p)
        E_s = ionic_singlet_energy(**p)
        T = E_s / 20.0
        chi = susceptibility(**p, T=T)
        chi_bb = bleaney_bowers(J, T)
        assert abs(chi - chi_bb) / chi_bb < 1e-3, (name, chi, chi_bb)
        T_plat = 100.0 * J                                            # J/T = 0.01
        if T_plat <= E_s / 20.0:                                      # deep-spin window exists
            assert abs(susceptibility(**p, T=T_plat) * T_plat - 0.5) < 2e-3, name


def test_G3_curie_weiss_theta_is_minus_J_over_4():
    """theta = -J/4 analytically, and a C/(T-theta) fit of the exact chi over the spin window
    recovers it (<5%) for the well-separated members. The iodides have no clean window (G4)."""
    for name, p in NB3X8_CLUSTERS.items():
        J = dimer_exchange_analytic(**p)
        assert abs(curie_weiss_theta(J) - (-J / 4.0)) < 1e-12
        E_s = ionic_singlet_energy(**p)
        if E_s / max(J, 1e-9) < 10.0:                                 # iodides: skip (finding)
            continue
        # fit 1/chi = (T - theta)/C over the spin window [5J, E_s/10]
        lo = max(5.0 * J, 1.0)
        Ts = np.linspace(lo, E_s / 10.0, 40)
        inv = 1.0 / susceptibility(**p, T=Ts)
        slope, intercept = np.polyfit(Ts, inv, 1)                     # 1/chi = slope*T + intercept
        theta = -intercept / slope
        assert abs(theta - (-J / 4.0)) < 0.05 * max(J, 1e-3), (name, theta, -J / 4.0)


def test_G4_charge_scale_boundary_and_iodide_ordering():
    """DEFINITION OF DONE: the 5% Bleaney-Bowers deviation temperature tracks the CHARGE scale
    (T5/E_s in [0.38,0.46] for all 10 sets), and the iodides break BB at the lowest reduced
    temperature (T5/J in [1.2,1.5] for Nb3I8 LT-bulk), strictly ordered I < Br < Cl < F by E_s/J."""
    ratios_es, t5_over_j = {}, {}
    for name, p in NB3X8_CLUSTERS.items():
        J = dimer_exchange_analytic(**p)
        E_s = ionic_singlet_energy(**p)
        T5 = bb_deviation_temperature(**p, tol=0.05)
        assert 0.38 < T5 / E_s < 0.46, (name, T5 / E_s)
        ratios_es[name] = E_s / J
        t5_over_j[name] = T5 / J
    assert 1.2 < t5_over_j["I  LT-bulk"] < 1.5, t5_over_j["I  LT-bulk"]
    # strict ordering of the LT-bulk halide series by reduced onset (small E_s/J breaks first)
    series = ["I  LT-bulk", "Br LT-bulk", "Cl LT-bulk", "F  LT-bulk"]
    onsets = [t5_over_j[s] for s in series]
    assert onsets == sorted(onsets), onsets
    assert t5_over_j["Br LT-bulk"] > 3.0 and t5_over_j["F  LT-bulk"] > 1e3
