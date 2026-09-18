#!/usr/bin/env python3
"""
Nb3(Br(1-x)I(x))8 alloy magnetic phase diagram -- where Heisenberg superexchange dies.

Spin-ODMD (``odmd_spin``) on a virtual-crystal alloy of the two Nb3X8 dimer endpoints: linearly
interpolate the cRPA parameters (U0, t, Us) of Nb3Br8 and Nb3I8, rebuild the generalized Hubbard
dimer at each alloy fraction x, and kick the half-filled ground state with the staggered
magnetization S1z - S2z. The single line that comes back is the interlayer exchange J(x); comparing
it to the textbook Heisenberg superexchange 4t^2/(U0-Us) locates the concentration where localized
moment physics stops being a valid description.

THE FINDING -- and it killed the spec's own prediction. SPEC v1.0 asserted the 20%-deviation
boundary sits at x_c in (0.45, 0.65). It does not:

    x_c = 0.296953   (Nb3Br0.59 I0.41 -- 41% iodine, NOT 45-65%)

Heisenberg is already 14.1% wrong at the bromide endpoint, so it only has ~6 points of headroom to
spend; it burns through them by x ~ 0.30, less than a third of the way across the alloy series.
The window in the spec corresponds to a 25-30% threshold, not the 20% the spec itself defined.

x_c has an exact closed form, which is what makes the kill airtight rather than a numerics
argument. The deviation depends on the parameters through the SINGLE dimensionless group
a = (2t / (U0 - Us))^2:

    Delta(a) = 2a / (sqrt(1 + 4a) - 1) - 1        so   Delta = D  <=>  |2t / (U0 - Us)| = sqrt(D + D^2)

which is linear in x on both sides and solves in closed form. Bisection and the closed form agree
to 1.7e-15. The same algebra exposes how SOFT the boundary is: x_c runs 0.05 -> 0.30 -> 0.64 as the
threshold moves 15% -> 20% -> 30%. x_c is a property of the threshold convention at least as much
as of the material, and no gate here pretends otherwise.

Secondary numbers (all in data/nb3x8_alloy_diagram.csv): J rises monotonically 119.11 -> 245.92 meV
across the series (150.15 meV at x_c), while the local-moment fraction ||Sz|psi0>||^2 falls
0.890 -> 0.759 -- charge fluctuations eating the moment is the same physics that breaks Heisenberg.

HONEST SCOPE (specs/SPEC_nb3x8_alloy_phase.md) -- read this before quoting x_c anywhere:
  * This is the VIRTUAL CRYSTAL APPROXIMATION, not an ab-initio alloy calculation. Linearly
    interpolating cRPA integrals models a fictitious averaged halide. There is NO site disorder, no
    local halide configuration around a given Nb3 trimer, no lattice relaxation or bowing, no
    change in the breathing-mode distortion, and no cRPA screening recomputed for the mixed
    compound. Real alloys bow; VCA cannot produce bowing by construction.
  * x_c is therefore a MODEL PREDICTION UNDER A STATED APPROXIMATION, not a measured or ab-initio
    property of Nb3(Br,I)8. It is falsifiable by a real cRPA calculation on a supercell alloy, and
    that is the point -- but it is not a substitute for one.
  * Inherits every caveat of odmd_spin: the ISOLATED interlayer dimer only (no in-plane kagome
    exchange, no band broadening -- see the coordination correction in nb3x8_gaps), density-density
    interactions only, exact statevector (no shots, no device noise).
  * The 20% deviation threshold is a convention, not a phase transition. Nothing is non-analytic at
    x_c; Delta(x) is smooth and strictly increasing everywhere. "Phase diagram" here means a
    validity boundary for an approximation, not a thermodynamic phase boundary.
  * Two endpoints, one interpolation path, one family. Nb3Cl8/Nb3F8 alloys are not covered.
"""
from __future__ import annotations

import math

import numpy as np
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.properties import AngularMomentum
from scipy.optimize import brentq

from hybrid_quantum_solver.model_hamiltonians import ModelIntegrals
from nb3x8_gaps import NB3X8_LT_BULK, dimer_cluster_integrals
from odmd_spin import (
    dimer_exchange_analytic,
    dimer_exchange_heisenberg,
    spin_excitation_lines,
)

# Alloy endpoints: x = 0 is Nb3Br8, x = 1 is Nb3I8 (LT bulk cRPA, arXiv:2501.10320 Table I).
ENDPOINTS = (NB3X8_LT_BULK["Nb3Br8"], NB3X8_LT_BULK["Nb3I8"])
_KEYS = ("U0", "t", "Us")

# Deviation threshold defining the Heisenberg-validity boundary. A CONVENTION (see HONEST SCOPE).
HEISENBERG_THRESHOLD = 0.20


def alloy_parameters(x: float) -> dict[str, float]:
    """Virtual-crystal (linear) interpolation of the cRPA parameters at alloy fraction ``x``."""
    br, i = ENDPOINTS
    return {k: (1.0 - x) * br[k] + x * i[k] for k in _KEYS}


def alloy_integrals(x: float) -> ModelIntegrals:
    """The half-filled Hubbard-dimer integrals ``(h1, eri, e_core, nelec, norb)`` of the alloy.

    Returns the repo's :class:`ModelIntegrals` (a dataclass carrying exactly those five fields),
    not the bare 4-tuple SPEC v1.0 wrote -- reusing ``dimer_cluster_integrals`` beats re-deriving
    the cluster, and every consumer in this repo takes a ``ModelIntegrals``.
    """
    return dimer_cluster_integrals(**alloy_parameters(x))


def alloy_ground_state(x: float):
    """``(mh, psi0)``: the alloy qubit Hamiltonian and its HALF-FILLED singlet ground state.

    Full-Fock-space diagonalization returns the one-electron bonding state for this model (see
    ``model_hamiltonians``), so the half-filled sector is selected by overlap with |HF>, exactly as
    ``odmd_spin``/``odmd_optical`` do.
    """
    base = alloy_integrals(x)
    mh = ModelIntegrals(base.h1, base.eri, 0.0, (1, 1), 2).to_hamiltonian()
    _, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    pops = np.abs(V.conj().T @ np.asarray(mh.hf_state().data, dtype=complex)) ** 2
    return mh, V[:, int(np.flatnonzero(pops > 1e-8)[0])]


def alloy_exchange_constant(x: float, n: int = 24) -> tuple[float, float]:
    """Spin-ODMD interlayer exchange ``(J, local_moment_fraction)`` of the alloy, in meV.

    ``n`` is the ODMD signal depth (SPEC v1.0 called it ``krylov_dim``). The v1.0 ``shots``
    argument is dropped: ``odmd_spin`` is an exact-statevector path with no shot model, and an
    ignored ``shots=`` parameter would be a lie in the signature.
    """
    mh, psi0 = alloy_ground_state(x)
    om, wt = spin_excitation_lines(mh, reference=psi0, n=n)
    return float(om[0]), float(wt[0])


def alloy_heisenberg_exchange(x: float) -> float:
    """Perturbative Heisenberg superexchange ``4 t(x)^2 / (U0(x) - Us(x))`` in meV."""
    return dimer_exchange_heisenberg(**alloy_parameters(x))


def alloy_heisenberg_deviation(x: float) -> float:
    """Fractional error of Heisenberg against the exact exchange, ``|J_Heis - J| / J``."""
    p = alloy_parameters(x)
    j = dimer_exchange_analytic(**p)
    return abs(dimer_exchange_heisenberg(**p) - j) / j


def critical_alloy_fraction(threshold: float = HEISENBERG_THRESHOLD) -> float:
    """Bisect for the ``x`` where the Heisenberg deviation first reaches ``threshold``."""
    return float(brentq(lambda x: alloy_heisenberg_deviation(x) - threshold, 0.0, 1.0, xtol=1e-14))


def critical_alloy_fraction_closed_form(threshold: float = HEISENBERG_THRESHOLD) -> float:
    """Closed-form ``x_c`` -- the ground-truth reference the bisection is checked against.

    Delta depends on the parameters only through ``a = (2t / (U0 - Us))^2``, via
    ``Delta = 2a/(sqrt(1+4a) - 1) - 1``. Inverting gives ``|2t| = sqrt(D + D^2) * (U0 - Us)``, and
    both sides are linear in ``x``, so the crossing is the root of a linear equation.
    """
    br, i = ENDPOINTS
    r = math.sqrt(threshold + threshold * threshold)
    t0, dt = 2.0 * br["t"], 2.0 * (i["t"] - br["t"])
    w0, dw = br["U0"] - br["Us"], (i["U0"] - i["Us"]) - (br["U0"] - br["Us"])
    return (r * w0 + t0) / (-dt - r * dw)  # |2t(x)| = -2t(x); t < 0 throughout


def map_magnetic_phase_diagram(steps: int = 21) -> list[dict[str, float]]:
    """Sweep ``x`` over ``[0, 1]``: parameters, J(x), J_Heis(x), deviation, local-moment fraction."""
    rows = []
    for x in np.linspace(0.0, 1.0, steps):
        p = alloy_parameters(float(x))
        j, moment = alloy_exchange_constant(float(x))
        jh = dimer_exchange_heisenberg(**p)
        rows.append({"x": float(x), **{k: p[k] for k in _KEYS}, "J_meV": j,
                     "J_heisenberg_meV": jh, "heisenberg_deviation": abs(jh - j) / j,
                     "local_moment_fraction": moment,
                     "heisenberg_valid": float(abs(jh - j) / j < HEISENBERG_THRESHOLD)})
    return rows


def total_spin_operator():
    """S^2 of the two-orbital dimer (qiskit-nature's AngularMomentum, JW-mapped)."""
    op = AngularMomentum(2).second_q_ops()["AngularMomentum"]
    return JordanWignerMapper().map(op).to_matrix(sparse=True).tocsc()


if __name__ == "__main__":
    import pathlib

    import polars as pl

    rows = map_magnetic_phase_diagram(steps=21)
    out = pathlib.Path("data/nb3x8_alloy_diagram.csv")
    out.parent.mkdir(exist_ok=True)
    pl.DataFrame(rows).write_csv(out)

    xc, xc_cf = critical_alloy_fraction(), critical_alloy_fraction_closed_form()
    print("Nb3(Br(1-x)I(x))8 virtual-crystal alloy -- interlayer exchange and Heisenberg validity")
    print(f"{'x':>5} {'t':>9} {'U0':>9} {'Us':>8} {'J [meV]':>9} {'J_Heis':>9} {'dev':>7} "
          f"{'moment':>7}")
    for r in rows:
        flag = "" if r["heisenberg_valid"] else "  <- Heisenberg invalid"
        print(f"{r['x']:5.2f} {r['t']:9.3f} {r['U0']:9.3f} {r['Us']:8.3f} {r['J_meV']:9.3f} "
              f"{r['J_heisenberg_meV']:9.3f} {r['heisenberg_deviation']*100:6.2f}% "
              f"{r['local_moment_fraction']:7.3f}{flag}")
    print(f"\nx_c (Delta = {HEISENBERG_THRESHOLD:.0%}) = {xc:.6f}  "
          f"[closed form {xc_cf:.6f}, agree to {abs(xc - xc_cf):.1e}]")
    print(f"  J(x_c) = {alloy_exchange_constant(xc)[0]:.2f} meV")
    print(f"  SPEC v1.0 predicted x_c in (0.45, 0.65) -- FALSIFIED by {0.45 - xc:.3f} in x.")
    print("\nx_c is soft: it is a property of the threshold convention as much as of the material.")
    for d in (0.15, 0.18, 0.20, 0.25, 0.30):
        print(f"  threshold {d:.0%} -> x_c = {critical_alloy_fraction(d):.4f}")
    print(f"\nwrote {out}  (VCA -- no disorder, no bowing, isolated dimer; see HONEST SCOPE)")
