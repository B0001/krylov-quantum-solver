"""
Acceptance gates G1-G4 for specs/SPEC_visibility_law.md (the predictive shot-cost law).

Test-first: ``visibility_law`` does not exist yet, so this file is RED until the spec is
implemented. Claim: the thrice-recorded qualitative visibility rule is a quantitative,
TRANSFERABLE law -- for the unnormalized correlator C_k = <psi0|O e^{-iktauH} O|psi0>, a line of
weight w costs shots* = 2(2-1/dim)/sigma*^2 with sigma* = w sqrt(dm)/(c(sqrt d + sqrt m)), i.e.
shots* ~ 1/(w^2 K). Gated three ways: the -2 log-log slope over four orders of magnitude in w,
single-calibration transfer to every other line (including components of a multi-line signal),
and the 1/K depth scaling. Boundary (found in probing): line attribution needs a tolerance below
the line spacing, else the strong neighbor masquerades as the weak line (45x too early).

All RNG seeded -> deterministic. PySCF/qiskit, no block2; `make gates` isolates it.
"""
import numpy as np

from hybrid_quantum_solver.model_hamiltonians import ModelIntegrals
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from nb3x8_gaps import NB3X8_LT_BULK, dimer_cluster_integrals
from odmd_optical import dimer_polarization
from odmd_spectral import ladder_operator, reference_signal
from visibility_law import crossover_shots, detect_line, predicted_shots

ORDER = ("Nb3F8", "Nb3Cl8", "Nb3Br8", "Nb3I8")
_CACHE = {}


def _dimer_correlator(name, K=16):
    if (name, K) not in _CACHE:
        base = dimer_cluster_integrals(**NB3X8_LT_BULK[name])
        mh = ModelIntegrals(base.h1, base.eri, 0.0, (1, 1), 2).to_hamiltonian()
        w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
        psi_hf = np.asarray(mh.hf_state().data, dtype=complex)
        pops = np.abs(V.conj().T @ psi_hf) ** 2
        psi0 = V[:, int(np.flatnonzero(pops > 1e-8)[0])]
        s, tau, mu, nrm2 = reference_signal(mh, dimer_polarization() @ psi0, K)
        _CACHE[(name, K)] = (nrm2 * s, tau, nrm2)          # line at 0 in the centered frame
    return _CACHE[(name, K)]


def _h4_correlator():
    if "h4" not in _CACHE:
        mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7")
        psi_hf = np.asarray(mh.hf_state().data, dtype=complex)
        raw = ladder_operator("-", 0, 8) @ psi_hf
        s, tau, mu, nrm2 = reference_signal(mh, raw, 32)
        w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
        pops = np.abs(V.conj().T @ (raw / np.sqrt(nrm2))) ** 2 * nrm2
        agg = {}
        for e, wt in zip(w, pops):
            key = round(float(e), 8)
            agg[key] = agg.get(key, 0.0) + float(wt)
        lines = sorted(agg.items(), key=lambda kv: -kv[1])[:3]
        _CACHE["h4"] = (nrm2 * s, tau, mu, lines)
    return _CACHE["h4"]


def _dimer_crossovers():
    if "stars" not in _CACHE:
        _CACHE["stars"] = {}
        for name in ORDER:
            C, tau, w = _dimer_correlator(name)
            _CACHE["stars"][name] = (w, crossover_shots(C, tau, 0.0, 16))
    return _CACHE["stars"]


def test_G1_minus_two_power_law():
    """Crossovers span >= 6 orders of magnitude in shots; log-log slope in [-2.05, -1.95]."""
    stars = _dimer_crossovers()
    ws = np.log10([stars[n][0] for n in ORDER])
    ss = np.log10([stars[n][1] for n in ORDER])
    assert ss.max() - ss.min() >= 6.0, (ss.min(), ss.max())
    slope = float(np.polyfit(ws, ss, 1)[0])
    assert -2.05 < slope < -1.95, slope


def test_G2_one_calibration_predicts_everything():
    """DEFINITION OF DONE: calibrate the prefactor on Nb3Br8 alone; every other line -- the
    three remaining dimers AND all three components of the multi-line H4 signal -- is predicted
    within a factor 1.5."""
    stars = _dimer_crossovers()
    w_br, star_br = stars["Nb3Br8"]
    cal = star_br / predicted_shots(w_br, 16, 16)
    for name in ("Nb3F8", "Nb3Cl8", "Nb3I8"):
        w, star = stars[name]
        ratio = star / (cal * predicted_shots(w, 16, 16))
        assert 1 / 1.5 < ratio < 1.5, (name, ratio)
    C, tau, mu, lines = _h4_correlator()
    for e_line, wt in lines:
        star = crossover_shots(C, tau, e_line - mu, 256)
        ratio = star / (cal * predicted_shots(wt, 32, 256))
        assert 1 / 1.5 < ratio < 1.5, (e_line, wt, ratio)


def test_G3_depth_buys_shots_linearly():
    """The EDGE component scales as 1/K: shots*(K=8)/shots*(K=32) in [3, 5] with attribution
    wide open (measured 3.83; law: 4). The full protocol's tight attribution adds a
    pole-accuracy cost concentrated at shallow depth (revised during implementation -- reality
    showed two separable effects): > 30% overhead at K=8, < 20% at K=32."""
    C8, tau8, _ = _dimer_correlator("Nb3I8", K=8)
    C32, tau32, _ = _dimer_correlator("Nb3I8", K=32)
    edge8 = crossover_shots(C8, tau8, 0.0, 16, tol_frac=0.5)
    edge32 = crossover_shots(C32, tau32, 0.0, 16, tol_frac=0.5)
    assert 3.0 < edge8 / edge32 < 5.0, edge8 / edge32
    full8 = crossover_shots(C8, tau8, 0.0, 16)
    full32 = crossover_shots(C32, tau32, 0.0, 16)
    assert full8 / edge8 > 1.3, full8 / edge8              # attribution costs at shallow K...
    assert full32 / edge32 < 1.2, full32 / edge32          # ...and is nearly free at depth


def test_G4_false_positives_and_the_attribution_boundary():
    """Protocol soundness: false-positive rate <= 2% at w=0. Boundary: with tol_frac=0.1 (above
    the 0.215 Ha line spacing) the H4 middle line's apparent crossover falls > 10x below the law
    -- the strong neighbor masquerades unless the tolerance sits below the line spacing."""
    sigma = 1e-3
    fp = np.mean([detect_line(np.zeros(16, complex), 1.0, sigma, 0.0, sd)
                  for sd in range(200)])
    assert fp <= 0.02, fp
    C, tau, mu, lines = _h4_correlator()
    e_mid, w_mid = lines[1]
    star_sloppy = crossover_shots(C, tau, e_mid - mu, 256, tol_frac=0.1)
    assert star_sloppy < predicted_shots(w_mid, 32, 256) / 10.0, star_sloppy
    star_tight = crossover_shots(C, tau, e_mid - mu, 256)      # default tol_frac=0.03
    assert star_tight > star_sloppy * 10.0, (star_tight, star_sloppy)
