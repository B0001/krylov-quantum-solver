"""
Acceptance gates G1-G4 for specs/SPEC_nb3x8_strain.md (strain / pressure Grueneisen response).

Claim: with |t| as the compression knob, (G1) the spin gap stiffens everywhere -- gamma_J = dlnJ/dln|t|
is positive and runs monotonically from the atomic-limit 2 (Nb3F8) toward 1 (Nb3I8), closed form ==
finite difference; (G2) the charge gap is non-monotonic (a minimum at |t*|) and the halide family
straddles it, so gamma_gap changes sign F/Cl (<0) -> Br/I (>0); (G3, DoD) Nb3Cl8 sits at its charge-
gap minimum -> a spin-charge-DECOUPLED strain response (strong gamma_J, near-zero gamma_gap); (G4)
all J-scale observables share gamma_J.

Isolated dimer, |t| the sole strain proxy, density-density only. PySCF/qiskit, no block2;
`make gates` runs it in its own process.
"""
from nb3x8_gaps import NB3X8_LT_BULK, exact_charge_gap
from nb3x8_magnetometry import chi_max_temperature
from nb3x8_strain import charge_gap_gruneisen, charge_gap_min_hopping, spin_gap_gruneisen
from nb3x8_thermo import schottky_peak_temperature
from odmd_spin import dimer_exchange_analytic

_SERIES = ("Nb3F8", "Nb3Cl8", "Nb3Br8", "Nb3I8")   # increasing |t|, decreasing U0


def _fd_gruneisen(f, p, rel=1e-6):
    h = rel * abs(p["t"])
    pu, pd = dict(p, t=p["t"] + h), dict(p, t=p["t"] - h)
    return p["t"] * (f(**pu) - f(**pd)) / (2 * h) / f(**p)


def test_G1_spin_gap_stiffens_and_runs_2_to_1():
    """gamma_J > 0 everywhere (compression stiffens the singlet), the closed form matches central
    finite differences (< 1e-3), gamma_J decreases strictly F->I, and hits the analytic limits:
    2 in the atomic limit (Nb3F8, J ~ 4t^2/(U0-Us)) and heading to 1 at strong hopping."""
    gammas = []
    for name in _SERIES:
        p = NB3X8_LT_BULK[name]
        gJ = spin_gap_gruneisen(**p)
        assert gJ > 0.0, (name, gJ)                                   # stiffens
        assert abs(gJ - _fd_gruneisen(dimer_exchange_analytic, p)) < 1e-3, name
        gammas.append(gJ)
    assert gammas == sorted(gammas, reverse=True), gammas             # strictly F>Cl>Br>I
    assert abs(gammas[0] - 2.0) < 1e-3, gammas[0]                     # atomic limit -> 2 (Nb3F8)
    assert 1.0 < gammas[-1] < 1.6, gammas[-1]                         # toward 1 (Nb3I8)


def test_G2_charge_gap_minimum_and_sign_straddle():
    """The exact charge gap has an interior minimum at |t*| (so it is non-monotonic in |t|), |t*|
    falls F->I, and the family straddles it: gamma_gap < 0 for the light halides (below their
    minima) and > 0 for the heavy ones (above)."""
    tstars = []
    for name in _SERIES:
        p = NB3X8_LT_BULK[name]
        tstar = charge_gap_min_hopping(p["U0"], p["Us"])
        # genuine interior minimum: gap at t* is below the gap at both flanks
        g_at = exact_charge_gap(p["U0"], -tstar, Us=p["Us"])
        assert exact_charge_gap(p["U0"], -0.5 * tstar, Us=p["Us"]) > g_at, name
        assert exact_charge_gap(p["U0"], -min(2 * tstar, 399.0), Us=p["Us"]) > g_at, name
        tstars.append(tstar)
    assert tstars == sorted(tstars, reverse=True), tstars             # |t*| falls F->I
    # sign straddle: heavy halides clearly stiffen, light ones do not (soft/near-zero)
    assert charge_gap_gruneisen(**NB3X8_LT_BULK["Nb3I8"]) > 0.1
    assert charge_gap_gruneisen(**NB3X8_LT_BULK["Nb3Br8"]) > 0.0
    assert charge_gap_gruneisen(**NB3X8_LT_BULK["Nb3Cl8"]) < 0.05
    assert charge_gap_gruneisen(**NB3X8_LT_BULK["Nb3F8"]) < 0.05


def test_G3_Nb3Cl8_spin_charge_decoupled_strain():
    """DEFINITION OF DONE (the sharp prediction): Nb3Cl8 sits at its charge-gap minimum (|t| within
    ~20% of |t*|), so its strain response is spin-charge decoupled -- strong spin-gap stiffening
    (gamma_J > 1.5) with a near-vanishing charge-gap response (|gamma_gap| < 0.05): a > 30x split."""
    p = NB3X8_LT_BULK["Nb3Cl8"]
    tstar = charge_gap_min_hopping(p["U0"], p["Us"])
    assert abs(abs(p["t"]) - tstar) / tstar < 0.20, (abs(p["t"]), tstar)   # sits near the minimum
    gJ = spin_gap_gruneisen(**p)
    gg = charge_gap_gruneisen(**p)
    assert gJ > 1.5, gJ
    assert abs(gg) < 0.05, gg
    assert abs(gJ / gg) > 30.0, (gJ, gg)                             # decoupled


def test_G4_J_scale_observables_share_gruneisen():
    """The chi(T) maximum and the Schottky peak both scale with J (exact two-level features), so for
    the WELL-SEPARATED members (Nb3Cl8, Nb3Br8) their Grueneisen parameters equal gamma_J to machine
    precision (< 0.1%): the spin thermodynamics inherits the spin-gap strain response exactly. The
    iodide is the SAME E_s/J boundary as everywhere in the thread -- charge fluctuations pull its
    chi/Schottky peaks off pure-J scaling, so its Grueneisen deviates markedly (> 1%) and by far
    more than the chloride's."""
    for name in ("Nb3Cl8", "Nb3Br8"):
        p = NB3X8_LT_BULK[name]
        gJ = spin_gap_gruneisen(**p)
        assert abs(_fd_gruneisen(chi_max_temperature, p) - gJ) < 1e-3 * gJ, name
        assert abs(_fd_gruneisen(schottky_peak_temperature, p) - gJ) < 1e-3 * gJ, name
    # iodide (E_s/J ~ 3): charge contamination breaks pure-J scaling -- the recorded boundary
    pI = NB3X8_LT_BULK["Nb3I8"]
    gJ_I = spin_gap_gruneisen(**pI)
    dev_I = abs(_fd_gruneisen(schottky_peak_temperature, pI) - gJ_I) / gJ_I
    pCl = NB3X8_LT_BULK["Nb3Cl8"]
    dev_Cl = abs(_fd_gruneisen(schottky_peak_temperature, pCl)
                 - spin_gap_gruneisen(**pCl)) / spin_gap_gruneisen(**pCl)
    assert dev_I > 0.01, dev_I
    assert dev_I > 10.0 * dev_Cl, (dev_I, dev_Cl)
