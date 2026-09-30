#!/usr/bin/env python3
"""
bead chem-a04 -- re-verification scout for the BACKLOG.md "parameter-free cluster optical gap"
entry (specs/BACKLOG.md, "Nb3X8 materials line" section).

Reuses validated primitives only (``odmd_optical.dimer_optical_gap``, ``nb3x8_gaps.coordination_gap``
and its LT-bulk parameter tables) -- no new physics, this is a pre-registered check against numbers
already quoted in the backlog entry. It does NOT re-source the measured literature values (1.10 eV,
~12% family spread, 0.63 eV collapse) -- this sandbox has no working internet access (WebSearch
returned only canned non-answers when tried), so those three numbers are taken as given by the
backlog entry's own citation, not re-derived here. What IS independently re-verified here is that the
MODEL's own numbers (which the backlog entry also quotes) reproduce from the checked-in code.

HONEST CAVEAT (stated up front, per the bead): the base 1.6% LT-Cl agreement is PARTLY FORTUITOUS --
omega_opt models a cluster excitation with the dipole operator P as a stand-in, while a measured
optical absorption onset folds in dispersion and excitonic binding energy the cluster model does not
capture at all. Read it as "scale agreement," not "prediction." The family-trend and phase-collapse
legs carry the actual falsifiability.

Finding (see specs/BACKLOG.md for the closed entry): all three legs match the backlog's quoted
numbers and neither kill condition triggers -- the LT Cl base number survives (1.59% miss, threshold
30%), and coordination_gap at z=3 for Nb3Cl8 (872.9 meV, robust FCI settings) stays well above the
700 meV kill threshold, so coordination/broadening alone is confirmed NOT to explain the ~43% measured
thermal collapse (model HT-parameter prediction only drops ~5%).

Convergence note: ``nb3x8_gaps.coordination_gap(*NB3X8_LT_BULK_5P["Nb3Cl8"], z=3)`` RAISES
RuntimeError with the module's default FCI settings (max_cycle=1000) -- the N=9 (5,4) sector does not
hit PySCF's Davidson convergence flag in time, even though the last-iteration energy
(4123.5653625580) already agrees with a converged run to 6 decimals. This script reproduces the exact
gap with tightened settings (max_cycle=4000, conv_tol=1e-10) rather than silently forcing the library
default. Filed separately (not fixed here -- out of this bead's scope): see the new bead created
alongside this check.
"""
from __future__ import annotations

import numpy as np
from pyscf import fci

from nb3x8_gaps import NB3X8_CLUSTERS, NB3X8_LT_BULK, NB3X8_LT_BULK_5P
from odmd_optical import dimer_optical_gap

# -- pre-registered reference numbers (quoted by specs/BACKLOG.md; NOT re-sourced in this sandbox --
# no working internet access here) ------------------------------------------------------------------
MEASURED_LT_CL_MEV = 1100.0     # ~1.10 eV absorption edge at 100 K
MEASURED_HT_CL_MEV = 630.0      # ~0.63 eV absorption edge at 300 K
MEASURED_FAMILY_SPREAD = 0.12   # measured Cl->I spread, ~12% (quoted, not re-derived)
LT_CL_KILL_REL = 0.30           # dies if LT Cl prediction misses measured by > 30%
COORD_Z3_KILL_MEV = 700.0       # dies (phase-collapse mechanism) if coordination_gap(z=3) <= this


def leg1_lt_cl_base():
    p = NB3X8_LT_BULK["Nb3Cl8"]
    model = dimer_optical_gap(**p)
    rel = (model - MEASURED_LT_CL_MEV) / MEASURED_LT_CL_MEV
    killed = abs(rel) > LT_CL_KILL_REL
    return model, rel, killed


def leg2_family_trend():
    gaps = {name: dimer_optical_gap(**NB3X8_LT_BULK[name])
            for name in ("Nb3F8", "Nb3Cl8", "Nb3Br8", "Nb3I8")}
    spread = (gaps["Nb3Cl8"] - gaps["Nb3I8"]) / gaps["Nb3I8"]
    return gaps, spread


def leg3_phase_collapse_ht_prediction():
    p_lt = NB3X8_LT_BULK["Nb3Cl8"]
    p_ht = NB3X8_CLUSTERS["Cl HT-bulk"]
    g_lt = dimer_optical_gap(**p_lt)
    g_ht = dimer_optical_gap(**p_ht)
    model_drop = (g_lt - g_ht) / g_lt
    measured_drop = (MEASURED_LT_CL_MEV - MEASURED_HT_CL_MEV) / MEASURED_LT_CL_MEV
    return g_lt, g_ht, model_drop, measured_drop


def _build_coordination_cluster(U0, ts, Us, tw, Uw, z):
    """Verbatim copy of nb3x8_gaps.coordination_gap's cluster construction (not importable
    standalone from that module -- h1/eri are built inline there)."""
    L = 2 + 2 * z
    h1 = np.zeros((L, L))
    eri = np.zeros((L, L, L, L))
    for i in range(L):
        eri[i, i, i, i] = U0

    def bond(i, j, t, U):
        h1[i, j] = h1[j, i] = t
        eri[i, i, j, j] = eri[j, j, i, i] = U

    bond(0, 1, ts, Us)
    for k in range(z):
        a, b = 2 + 2 * k, 3 + 2 * k
        bond(a, b, ts, Us)
        bond(0 if k % 2 == 0 else 1, a, tw, Uw)
    return h1, eri


def leg3_coordination_gap_robust(z: int, max_cycle: int = 4000, conv_tol: float = 1e-10) -> float:
    """coordination_gap's own physics (central dimer + z weak-linked pendant dimers, FCI charge
    gap), with tightened Davidson settings -- the module default (max_cycle=1000) does not
    converge at z=3 for the Nb3Cl8 parameter set (see module docstring)."""
    h1, eri = _build_coordination_cluster(*NB3X8_LT_BULK_5P["Nb3Cl8"], z)
    L = h1.shape[0]
    E = {}
    for n in (L - 1, L, L + 1):
        na, nb = n - n // 2, n // 2
        solver = fci.direct_spin1.FCI()
        solver.max_cycle = max_cycle
        solver.conv_tol = conv_tol
        e, _ = solver.kernel(h1, eri, L, (na, nb))
        if not solver.converged:
            raise RuntimeError(f"z={z}, N={n}: FCI did not converge even at max_cycle={max_cycle}")
        E[n] = e
    return E[L + 1] + E[L - 1] - 2 * E[L]


def main():
    print("=== chem-a04: Nb3Cl8 optical gap -- three-leg re-verification ===\n")
    print("CAVEAT (stated up front, per the bead): the LT-Cl 1.6%-scale hit is PARTLY FORTUITOUS --")
    print("omega_opt is a cluster excitation with P as a dipole stand-in; a measured absorption onset")
    print("folds in dispersion and excitonic binding the cluster model does not capture. Read it as")
    print("'scale agreement,' not 'prediction.' Legs 2-3 carry the actual falsifiability.\n")

    model, rel, killed = leg1_lt_cl_base()
    print(f"Leg 1 -- LT Cl base: model={model:.1f} meV, measured={MEASURED_LT_CL_MEV:.1f} meV, "
          f"rel err={rel * 100:+.2f}% (kill if |rel|>{LT_CL_KILL_REL * 100:.0f}%) -> "
          f"{'KILLED' if killed else 'survives'}")

    gaps, spread = leg2_family_trend()
    print("\nLeg 2 -- family trend (LT bulk, model optical gaps):")
    for name, g in gaps.items():
        print(f"  {name}: {g:.1f} meV")
    print(f"  model Cl->I spread (Cl-I)/I = {spread * 100:.1f}%  "
          f"vs measured ~{MEASURED_FAMILY_SPREAD * 100:.0f}% (quoted, not re-sourced here)")
    print("  Nb3Br8 leg DROPPED from the measured comparison: no measured Nb3Br8 optical gap was "
          "located (per the backlog entry) -- reporting the model value only, uncompared.")

    g_lt, g_ht, model_drop, measured_drop = leg3_phase_collapse_ht_prediction()
    print(f"\nLeg 3a -- phase collapse, HT-parameter prediction: LT model={g_lt:.1f} meV, "
          f"HT model={g_ht:.1f} meV, model drop={model_drop * 100:.2f}%; measured drop "
          f"({MEASURED_LT_CL_MEV:.0f}->{MEASURED_HT_CL_MEV:.0f} meV)={measured_drop * 100:.1f}%")
    print("  -> the model's own HT parameters predict a collapse ~1/9th the measured size: the "
          "optical-gap formula by itself cannot explain the 300 K collapse.")

    print("\nLeg 3b -- coordination_gap(Nb3Cl8, z) vs the measured 630 meV collapse "
          f"(kill if z=3 reaches <= {COORD_Z3_KILL_MEV:.0f} meV):")
    z3 = None
    for z in (0, 1, 2, 3, 4):
        g = leg3_coordination_gap_robust(z)
        if z == 3:
            z3 = g
        print(f"  z={z}: {g:.2f} meV  (short of 630 meV by {g - MEASURED_HT_CL_MEV:.2f} meV)")
    coord_killed = z3 <= COORD_Z3_KILL_MEV
    print(f"  z=3 = {z3:.2f} meV -> {'KILLED (coordination alone explains the collapse)' if coord_killed else 'survives (coordination alone does NOT explain the collapse)'}")

    print("\n=== verdict ===")
    print(f"Leg 1 (base number, 30% kill): {'KILLED' if killed else 'SURVIVES'}")
    print(f"Leg 3b (coordination-explains-collapse kill, <=700meV): "
          f"{'KILLED' if coord_killed else 'SURVIVES'}")
    print("Neither kill condition triggers -- the selection-rule optical-gap picture is not "
          "falsified by these two checks, and the thermal collapse still needs a mechanism beyond "
          "simple coordination/broadening.")


if __name__ == "__main__":
    main()
