#!/usr/bin/env python3
"""
Finite-temperature magnetic susceptibility chi(T) of the Nb3X8 downfolded dimer clusters.

A thermodynamic study (not a method rung -- exact trace at 4 qubits is trivial; no quantum
advantage, and claiming one would be dishonest). From the same exactly-diagonalizable
extended-Hubbard dimer as nb3x8_gaps.py (cRPA parameters of arXiv:2501.10320), the molar magnetic
susceptibility follows from the Van Vleck zero-field trace over the half-filling (N=2) spectrum:

    chi(T) = <S_z,tot^2>_thermal / T          (reduced units: k_B = 1, energies & T in meV,
                                               chi in meV^-1, g=2 folded into the emu conversion)

REPRODUCTION of standard magnetochemistry applied to new ab-initio parameters: the
Bleaney-Bowers singlet-triplet law (Bleaney & Bowers, Proc. R. Soc. A 214, 451, 1952; Kahn,
Molecular Magnetism, 1993) and its deviation as T approaches the charge scale are the exact
two-site Hubbard-dimer thermodynamics (Anderson, Phys. Rev. 115, 2, 1959; Carrascal et al.,
arXiv:1502.05038) -- known, expected, NOT a discovery. Prior Nb3X8 chi(T) exists experimentally
(Sheckelton et al., Inorg. Chem. Front. 4, 481, 2017 -- the ~90 K Nb3Cl8 singlet transition;
Haraguchi et al., Inorg. Chem. 56, 3483, 2017) and theoretically (Grytsiuk/Rosner,
arXiv:2305.04854). This module contributes the family-wide exact-spectrum table and the clean
charge-scale validity boundary, cross-checking that known physics.

THE FINDING (specs/SPEC_nb3x8_susceptibility.md G4): the pure-spin (Bleaney-Bowers) picture holds
up to a temperature set by the CHARGE scale E_s (the first ionic singlet), not the exchange J:
the 5% deviation temperature is T ~ 0.40 * E_s across the whole family, so the iodides
(E_s/J ~ 3, spin and charge least separated) break Bleaney-Bowers at the lowest reduced
temperature (T5/J ~ 1.3, vs > 3 for Br and > 10^3 for F).

HONEST SCOPE: isolated single dimer -- no inter-dimer coupling, no phonons, no structural
transition (real Nb3Cl8 has a first-order one at ~90 K); density-density interactions only; the
emu conversion assumes g=2 (spin-only). A reference table, not a solid-state prediction.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp

from nb3x8_gaps import dimer_cluster_integrals

_MAPPER = JordanWignerMapper()

# JW block order for the 2-orbital dimer: q0,q1 = up orbitals; q2,q3 = down orbitals.
_SZ = _MAPPER.map(FermionicOp({"+_0 -_0": 0.5, "+_1 -_1": 0.5, "+_2 -_2": -0.5, "+_3 -_3": -0.5},
                              num_spin_orbitals=4)).to_matrix().real
_N = _MAPPER.map(FermionicOp({f"+_{i} -_{i}": 1.0 for i in range(4)},
                             num_spin_orbitals=4)).to_matrix().real

# chi[emu/mol] = EMU_PER_REDUCED * chi_reduced[meV^-1].  Derivation: chi_reduced = <Sz^2>/T[meV];
# chi[emu/mol] = (N_A g^2 mu_B^2 / k_B) <Sz^2> / T[K] = 1.5006 <Sz^2>/T[K] (g=2); with
# T[K] = T[meV]/0.0861733, this is 1.5006 * 0.0861733 * chi_reduced = 0.12931 * chi_reduced.
EMU_PER_REDUCED = 1.5006 * 0.0861733     # ~= 0.12931 emu*meV/mol


@lru_cache(maxsize=None)
def n2_spectrum(U0: float, t: float, Us: float):
    """(energies, sz2) of the half-filling (N=2) sector, energies in meV relative to its ground.

    ``sz2[i]`` = <i|S_z,tot^2|i> (diagonal, since [H, S_z] = 0). The structure is: singlet ground
    (sz2=0), 3-fold triplet at J (sz2 = 1,0,1), then the ionic/charge-transfer singlets at ~E_s.
    Cached: the Hamiltonian build + diagonalization is reused across the many chi(T)/onset calls.
    """
    H = dimer_cluster_integrals(U0, t, Us).to_hamiltonian().qubit_hamiltonian.to_matrix()
    w, V = np.linalg.eigh(H)
    n = np.real(np.einsum("ji,jk,ki->i", V.conj(), _N, V))
    keep = np.abs(n - 2.0) < 1e-6
    e = w[keep].real
    sz2 = np.real(np.einsum("ji,jk,ki->i", V[:, keep].conj(), _SZ @ _SZ, V[:, keep]))
    order = np.argsort(e)
    e, sz2 = e[order], sz2[order]
    return e - e[0], sz2


def ionic_singlet_energy(U0: float, t: float, Us: float) -> float:
    """E_s: the first N=2 level above the J-triplet (the charge/ionic scale that bounds the
    pure-spin regime). The triplet is the 3-fold degenerate manifold at J."""
    e, _ = n2_spectrum(U0, t, Us)
    triplet = e[3]                                    # top of the 3-fold triplet (e[1:4] ~ J)
    return float(e[e > triplet + 1e-6][0])


def susceptibility(U0: float, t: float, Us: float, T):
    """Exact reduced molar susceptibility chi(T) = <S_z,tot^2>_thermal / T (meV^-1); T in meV,
    scalar or array."""
    e, sz2 = n2_spectrum(U0, t, Us)
    T = np.asarray(T, dtype=float)
    scalar = T.ndim == 0
    Tv = np.atleast_1d(T)
    beta_e = e[:, None] / Tv                            # (levels, nT)
    boltz = np.exp(-(beta_e - beta_e.min(axis=0)))     # shift for numerical stability
    sz2_avg = (sz2[:, None] * boltz).sum(axis=0) / boltz.sum(axis=0)
    chi = sz2_avg / Tv
    return float(chi[0]) if scalar else chi


def bleaney_bowers(J: float, T):
    """Analytic reduced susceptibility of an isotropic S=1/2 dimer: chi = (2/T)/(3 + e^{J/T})
    (Bleaney-Bowers; J = E_triplet - E_singlet, antiferromagnetic for J > 0)."""
    T = np.asarray(T, dtype=float)
    chi = (2.0 / T) / (3.0 + np.exp(J / T))
    return float(chi) if chi.ndim == 0 else chi


def curie_weiss_theta(J: float) -> float:
    """Weiss temperature of the S=1/2 AF dimer, theta = -J/4 (reduced, meV; high-T expansion of
    Bleaney-Bowers, H = +J S1.S2, z=1)."""
    return -J / 4.0


def bb_deviation_temperature(U0: float, t: float, Us: float, tol: float = 0.05,
                             t_max_factor: float = 1.5) -> float:
    """Lowest T (meV) at which |chi_exact - chi_BB| / chi_BB first reaches ``tol``.

    Bisection between the spin plateau and ~t_max_factor*E_s; chi_BB uses the exact J. Returns
    the deviation onset -- gated (G4) to be ~0.40*E_s across the family (charge-scale, not J)."""
    from odmd_spin import dimer_exchange_analytic

    J = dimer_exchange_analytic(U0, t, Us)
    E_s = ionic_singlet_energy(U0, t, Us)

    def dev(T):
        return abs(susceptibility(U0, t, Us, T) - bleaney_bowers(J, T)) / bleaney_bowers(J, T)

    # Find the FIRST upward crossing of ``tol`` by a coarse (vectorized) scan, then bisect inside
    # that bracket (the physical onset temperature). Start well below the onset (~0.4*E_s): at low
    # T both chi and chi_BB vanish with the same singlet-triplet exponential, so the deviation ->
    # 0 there. (Starting at 2J would miss the onset for the iodides, whose onset sits below 2J.)
    grid = np.linspace(E_s / 50.0, t_max_factor * E_s, 2000)
    d = np.abs(susceptibility(U0, t, Us, grid) - bleaney_bowers(J, grid)) / bleaney_bowers(J, grid)
    above = np.flatnonzero(d >= tol)
    if above.size == 0:
        return float(grid[-1])
    i = above[0]
    if i == 0:
        return float(grid[0])
    lo, hi = grid[i - 1], grid[i]                      # dev(lo) < tol <= dev(hi)
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if dev(mid) >= tol:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    from nb3x8_gaps import NB3X8_CLUSTERS
    from odmd_spin import dimer_exchange_analytic

    T_room = 25.85                                     # meV ~ 300 K
    print("Nb3X8 dimer magnetic susceptibility -- exact N=2 Boltzmann trace")
    print("(isolated cluster; reproduction of Bleaney-Bowers with ab-initio (t,U) parameters)")
    print(f"{'set':>11} | {'J(meV)':>8} | {'E_s(meV)':>8} | {'E_s/J':>7} | "
          f"{'chi(300K)':>10} | {'mu_eff':>6} | {'theta_CW':>8} | {'T5%/E_s':>7} | {'T5%/J':>8}")
    for name, p in NB3X8_CLUSTERS.items():
        J = dimer_exchange_analytic(**p)
        E_s = ionic_singlet_energy(**p)
        chi_r = susceptibility(**p, T=T_room)
        mu_eff = np.sqrt(3.0 * T_room * chi_r)
        T5 = bb_deviation_temperature(**p, tol=0.05)
        print(f"{name:>11} | {J:8.3f} | {E_s:8.1f} | {E_s / J:7.1f} | "
              f"{chi_r:10.3e} | {mu_eff:6.3f} | {curie_weiss_theta(J):8.3f} | "
              f"{T5 / E_s:7.3f} | {T5 / J:8.2f}")
    print(f"\nemu conversion: chi[emu/mol] = {EMU_PER_REDUCED:.5f} * chi_reduced[meV^-1]  (g=2)")
