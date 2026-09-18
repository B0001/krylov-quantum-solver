"""
Acceptance gates for specs/SPEC_nb3x8_alloy_phase.md (Nb3(Br1-xIx)8 alloy magnetic phase diagram).

The spec's headline prediction -- the Heisenberg-validity boundary sits at x_c in (0.45, 0.65) --
was FALSIFIED during implementation, and G3 below encodes the kill rather than the claim
(specs/README.md step 5). The measured boundary is x_c = 0.296953: 41% iodine, not 45-65%.

  G1  endpoint consistency (unchanged, passes): spin-ODMD on the alloy at x = 0 and x = 1
      reproduces odmd_spin's own gated Nb3Br8/Nb3I8 exchange constants 119.11 / 245.92 meV, and
      the published 4-figure values 119.1 / 245.9, far inside the specced 0.1%.
  G2  smooth, strictly monotonic, non-singular interpolation (unchanged, passes) -- with the
      nuance the spec did not state: U0, t and Us are all strictly DECREASING in x, so |t| rises
      while U0 - Us falls. The two move the deviation the same way, which is why Delta climbs so
      fast.
  G3  THE KILL. x_c is pinned at the measured 0.296953 (+/- 1e-6) against an exact closed form,
      and asserted to lie OUTSIDE the spec's (0.45, 0.65) window. The tolerance was NOT widened.
  G3b why the window was wrong: x_c is extremely soft in the threshold convention (0.052 at 15%,
      0.297 at 20%, 0.639 at 30%). The spec's window is what a 25-30% threshold would give, not
      the 20% the spec itself defined.
  G4  singlet ground state everywhere (unchanged, passes), S^2 < 1e-6 across the sweep -- and the
      S^2 operator is proved non-vacuous (the kicked triplet carries S^2 = 2 at exactly J).

Deterministic: no RNG anywhere (exact statevector + closed forms).
PySCF/qiskit, no block2; `make gates` runs it in its own process.
"""
import numpy as np

from nb3x8_alloy import (
    HEISENBERG_THRESHOLD,
    alloy_exchange_constant,
    alloy_ground_state,
    alloy_heisenberg_deviation,
    alloy_heisenberg_exchange,
    alloy_parameters,
    critical_alloy_fraction,
    critical_alloy_fraction_closed_form,
    map_magnetic_phase_diagram,
    total_spin_operator,
)
from nb3x8_gaps import NB3X8_LT_BULK
from odmd_spin import dimer_exchange_analytic, dimer_staggered_moment

# The measured boundary that replaces the spec's (0.45, 0.65) window.
X_C_MEASURED = 0.2969530815
_CACHE = {}


def _diagram():
    """The 21-point sweep, built once (each point is a 16x16 eigh -- cheap, but not free x4)."""
    if "rows" not in _CACHE:
        _CACHE["rows"] = map_magnetic_phase_diagram(steps=21)
    return _CACHE["rows"]


def test_G1_endpoint_consistency():
    """x = 0 and x = 1 reproduce Nb3Br8 / Nb3I8 exactly: 119.108 / 245.920 meV.

    Measured relative error vs odmd_spin's gated analytic exchange: < 1e-14 (the ODMD line IS the
    exact singlet-triplet splitting). Against the published 4-figure 119.1 / 245.9 meV: 6.8e-5 and
    8.0e-5 -- the specced 0.1% holds with three orders of magnitude to spare, so the rounding in
    those figures never mattered.
    """
    for x, name, published in ((0.0, "Nb3Br8", 119.1), (1.0, "Nb3I8", 245.9)):
        assert alloy_parameters(x) == NB3X8_LT_BULK[name], x
        j, _ = alloy_exchange_constant(x)
        exact = dimer_exchange_analytic(**NB3X8_LT_BULK[name])
        assert abs(j - exact) < 1e-12 * exact, (name, j, exact)
        assert abs(j - published) < 1e-3 * published, (name, j, published)


def test_G2_smooth_strictly_monotonic_interpolation():
    """t(x), U0(x), Us(x) are strictly monotonic and J(x) is smooth with no singular points.

    All three parameters strictly DECREASE with x (t: -169.40 -> -218.20, U0: 1186.60 -> 787.00,
    Us: 342.00 -> 258.50). So |t| grows while the charge-transfer scale U0 - Us shrinks: the
    deviation's single group |2t/(U0-Us)| is strictly increasing (0.401 -> 0.826), which is the
    mechanism behind G3. J(x) itself is strictly increasing, 119.11 -> 245.92 meV, with exactly
    one spectral line at every x (no degeneracy collapse, no singular matrix).
    """
    rows = _diagram()
    for key in ("t", "U0", "Us"):
        seq = [r[key] for r in rows]
        assert all(b < a for a, b in zip(seq, seq[1:])), (key, seq)
    ratio = [abs(2 * r["t"] / (r["U0"] - r["Us"])) for r in rows]
    assert all(b > a for a, b in zip(ratio, ratio[1:])), ratio
    assert abs(ratio[0] - 0.4011) < 1e-3 and abs(ratio[-1] - 0.8257) < 1e-3, (ratio[0], ratio[-1])

    js = [r["J_meV"] for r in rows]
    assert all(np.isfinite(js)) and all(b > a for a, b in zip(js, js[1:])), js
    # second differences of a smooth curve on a uniform grid stay tiny next to the curve's range
    assert np.max(np.abs(np.diff(js, 2))) < 0.01 * (js[-1] - js[0]), np.diff(js, 2)
    # local-moment fraction falls monotonically 0.890 -> 0.759: moment eaten by charge fluctuations
    moments = [r["local_moment_fraction"] for r in rows]
    assert all(b < a for a, b in zip(moments, moments[1:])), moments
    assert abs(moments[0] - 0.8900) < 1e-3 and abs(moments[-1] - 0.7590) < 1e-3, moments


def test_G3_critical_fraction_falsifies_the_specced_window():
    """THE KILL -- replaces SPEC v1.0 G3 ("x_c in (0.45, 0.65)").

    Measured x_c = 0.296953 at the spec's own 20% threshold, pinned here to +/- 1e-6 against an
    exact closed form (bisection and closed form agree to 2.4e-15). That is 0.153 BELOW the
    specced window, which this gate asserts explicitly so the falsification cannot silently decay
    back into a pass. Delta(x) is strictly increasing, so the crossing is unique -- there is no
    second root hiding inside (0.45, 0.65).
    """
    xc = critical_alloy_fraction()
    assert abs(xc - critical_alloy_fraction_closed_form()) < 1e-12, xc
    assert abs(xc - X_C_MEASURED) < 1e-6, xc
    assert not (0.45 < xc < 0.65), f"spec v1.0 window (0.45, 0.65) would have held: x_c={xc}"

    grid = np.linspace(0.0, 1.0, 2001)
    dev = [alloy_heisenberg_deviation(float(x)) for x in grid]
    assert all(b > a for a, b in zip(dev, dev[1:])), "Delta(x) not strictly increasing"
    assert abs(dev[0] - 0.1410) < 1e-3 and abs(dev[-1] - 0.4653) < 1e-3, (dev[0], dev[-1])
    # the single crossing, located on the grid, agrees with the closed form
    crossings = np.flatnonzero(np.diff(np.array(dev) >= HEISENBERG_THRESHOLD))
    assert len(crossings) == 1 and abs(grid[crossings[0]] - xc) <= grid[1] - grid[0], crossings
    # and J there -- the number to quote alongside x_c
    assert abs(alloy_exchange_constant(xc)[0] - 150.1466) < 1e-3


def test_G3b_the_boundary_is_soft_in_the_threshold():
    """WHY the specced window was wrong: x_c is a property of the threshold convention.

    Heisenberg is already 14.1% off at the bromide endpoint, so a 20% threshold has ~6 points of
    headroom and is crossed early. Measured x_c vs threshold: 15% -> 0.0524, 18% -> 0.2074,
    20% -> 0.2970, 25% -> 0.4861, 30% -> 0.6395. The spec's (0.45, 0.65) is exactly the 25-30%
    band -- the window was consistent with a stricter-sounding but weaker criterion than the 20%
    the spec defined. This gate pins that diagnosis, so "just move the threshold" is visible as
    the choice it is rather than a tuning knob.
    """
    measured = {0.15: 0.0524, 0.18: 0.2074, 0.20: 0.2970, 0.25: 0.4861, 0.30: 0.6395}
    got = {d: critical_alloy_fraction(d) for d in measured}
    for d, expect in measured.items():
        assert abs(got[d] - expect) < 1e-4, (d, got[d], expect)
    # ~0.039 in x per point of threshold: a soft boundary, not a sharp transition
    slope = (got[0.30] - got[0.15]) / (0.30 - 0.15)
    assert 3.5 < slope < 4.5, slope
    # the diagnosis: 25% and 30% thresholds DO land in the spec's window; the specced 20% does not
    assert all(0.45 < got[d] < 0.65 for d in (0.25, 0.30)), got
    assert not (0.45 < got[0.20] < 0.65), got[0.20]


def test_G4_ground_state_is_a_singlet_everywhere():
    """S^2 < 1e-6 at every x, so the staggered-magnetization kick is structurally valid.

    Measured max S^2 over the 21-point sweep: ~1e-30 (machine zero). The operator is NOT
    vacuously zero: the state the kick actually populates is the m = 0 triplet with S^2 = 2.000,
    sitting exactly J above the singlet -- which is both the non-vacuity check and an independent
    confirmation of J from raw eigenvalues rather than from ODMD.
    """
    s2 = total_spin_operator()
    sz = dimer_staggered_moment()
    worst = 0.0
    for r in _diagram():
        mh, psi0 = alloy_ground_state(r["x"])
        worst = max(worst, abs(float(np.real(psi0.conj() @ (s2 @ psi0)))))

        # non-vacuity + independent J: the kicked state is a pure triplet at energy E0 + J
        kicked = sz @ psi0
        nrm = float(np.linalg.norm(kicked))
        assert nrm > 0.5, (r["x"], nrm)
        kicked = kicked / nrm
        assert abs(float(np.real(kicked.conj() @ (s2 @ kicked))) - 2.0) < 1e-9, r["x"]
        h = mh.qubit_hamiltonian.to_matrix()
        gap = float(np.real(kicked.conj() @ (h @ kicked) - psi0.conj() @ (h @ psi0)))
        assert abs(gap - r["J_meV"]) < 1e-9 * r["J_meV"], (r["x"], gap, r["J_meV"])
    assert worst < 1e-6, worst


def test_G5_heisenberg_always_overestimates():
    """Rider: J_Heis > J at every x -- the perturbative form is a one-sided error, never a wobble.

    4t^2/(U0-Us) is the leading term of sqrt(W^2/4 + 4t^2) - W/2, whose remainder is strictly
    negative, so the deviation is a systematic overestimate. Gated because it is what lets G3
    treat Delta as a monotone validity measure rather than an absolute error that could cross
    zero.
    """
    for r in _diagram():
        assert alloy_heisenberg_exchange(r["x"]) > r["J_meV"], r["x"]
        assert r["J_heisenberg_meV"] > r["J_meV"], r["x"]
