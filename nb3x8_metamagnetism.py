#!/usr/bin/env python3
"""
Metamagnetism of the Nb3X8 dimer -- a field-driven singlet-to-triplet ground-state crossing.

Adds a uniform Zeeman term -h*Sz,tot to the same exactly-diagonalizable extended-Hubbard dimer
used by nb3x8_susceptibility.py/odmd_spin.py (cRPA parameters of arXiv:2501.10320). Because
[H0, Sz] = 0, the fully-polarized Sz=+1 triplet member is the ONLY two-same-spin Slater
determinant in a 2-orbital space -- so it is a field-independent-apart-from-the-rigid-shift exact
eigenstate, while the Sz=0 singlet ground state is untouched by the field. That predicts an exact
closed-form crossing at h_c = J (the interlayer exchange already gated in SPEC_odmd_spin),
reproduced here by DIRECT diagonalization of the full field-augmented matrix -- not assumed from
block symmetry.

THE FINDING (specs/SPEC_nb3x8_metamagnetism.md G4): expressed as a magnetic field via
B = h/(g*mu_B) (g=2), Nb3Cl8 (572 T) and Nb3Br8 (1029 T) sit above the best documented
non-destructive pulsed-field record (100.75 T, Los Alamos 2012) but within the best documented
destructive electromagnetic flux-compression record (1200 T, U. Tokyo 2018); Nb3I8 (2124 T)
exceeds even that. A real observation would need either the lattice's cooperative coupling to be
much weaker than this isolated-cluster J -- consistent with the 2.3-5.3x Tc overcoupling already
found in SPEC_nb3x8_magnetometry -- or destructive megagauss techniques.

HONEST SCOPE: isolated single dimer, T=0 ground-state crossing only, density-density interactions
only, g=2 assumed (same convention as nb3x8_susceptibility.py). Nb3F8's h_c (~0.44 T) is not a
physical prediction -- its J (0.051 meV) is already flagged in SPEC_odmd_spin as below the
model's own neglected non-density-density terms.
"""
from __future__ import annotations

import numpy as np
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp

from nb3x8_gaps import dimer_cluster_integrals

_MAPPER = JordanWignerMapper()

# JW block order for the 2-orbital dimer: q0,q1 = up orbitals; q2,q3 = down orbitals (matches
# nb3x8_susceptibility.py's convention).
_SZ = _MAPPER.map(FermionicOp({"+_0 -_0": 0.5, "+_1 -_1": 0.5, "+_2 -_2": -0.5, "+_3 -_3": -0.5},
                              num_spin_orbitals=4)).to_matrix()
_N = _MAPPER.map(FermionicOp({f"+_{i} -_{i}": 1.0 for i in range(4)},
                             num_spin_orbitals=4)).to_matrix()

# g=2 Bohr magneton in meV/T (same g=2 spin-only convention as nb3x8_susceptibility.py).
G_MU_B = 0.115767636


def zeeman_ground_state(U0: float, t: float, Us: float, h: float) -> tuple[float, float]:
    """(energy, <Sz>) of the N=2-sector ground state of H0 - h*Sz,tot, by DIRECT diagonalization
    of the full field-augmented matrix (no assumed block structure)."""
    H0 = dimer_cluster_integrals(U0, t, Us).to_hamiltonian().qubit_hamiltonian.to_matrix()
    Hh = H0 - h * _SZ
    w, V = np.linalg.eigh(Hh)
    n = np.real(np.einsum("ji,jk,ki->i", V.conj(), _N, V))
    idx = np.flatnonzero(np.abs(n - 2.0) < 1e-6)
    gi = idx[np.argmin(w[idx])]
    sz = float(np.real(V[:, gi].conj() @ _SZ @ V[:, gi]))
    return float(w[gi]), sz


def magnetization(U0: float, t: float, Us: float, h: float) -> float:
    """<Sz> of the ground state at field h (meV) -- the metamagnetic order parameter."""
    return zeeman_ground_state(U0, t, Us, h)[1]


def critical_field_numeric(U0: float, t: float, Us: float, tol: float = 1e-9) -> float:
    """Bisect on <Sz> to find the field h_c (meV) where the ground state crosses from singlet
    (Sz=0) to the polarized triplet member (Sz=1). No closed form assumed."""
    from odmd_spin import dimer_exchange_analytic

    J = dimer_exchange_analytic(U0, t, Us)
    lo, hi = 0.0, 2.0 * J
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if magnetization(U0, t, Us, mid) < 0.5:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def critical_field_tesla(U0: float, t: float, Us: float) -> float:
    """B_c = J / (g*mu_B) in Tesla (g=2), using the closed-form J directly."""
    from odmd_spin import dimer_exchange_analytic

    return dimer_exchange_analytic(U0, t, Us) / G_MU_B


if __name__ == "__main__":
    from nb3x8_gaps import NB3X8_LT_BULK
    from odmd_spin import dimer_exchange_analytic

    NONDESTRUCTIVE_RECORD_T = 100.75   # Los Alamos, March 2012
    DESTRUCTIVE_RECORD_T = 1200.0      # U. Tokyo, 2018

    print("Nb3X8 metamagnetism -- field-driven singlet->triplet crossing:")
    print(f"{'material':>8} | {'J (meV)':>10} | {'h_c numeric':>12} | {'B_c (T)':>10} | feasibility")
    for name, p in NB3X8_LT_BULK.items():
        J = dimer_exchange_analytic(**p)
        h_c = critical_field_numeric(**p)
        B_c = critical_field_tesla(**p)
        if B_c < NONDESTRUCTIVE_RECORD_T:
            note = "below even non-destructive record (below model's noise floor)"
        elif B_c < DESTRUCTIVE_RECORD_T:
            note = "needs destructive flux-compression; beyond non-destructive record"
        else:
            note = "beyond even the destructive world record"
        print(f"{name:>8} | {J:10.4f} | {h_c:12.6f} | {B_c:10.2f} | {note}")
