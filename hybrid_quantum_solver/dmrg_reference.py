#!/usr/bin/env python3
"""
dmrg_reference.py -- classical reference energies for active-space benchmarks.

For active spaces small enough that FCI is exact (≲ 16-18 orbitals), exact FCI is the gold
standard and is what `benchmark_n2.py` uses. For LARGER active spaces (transition-metal scale,
beyond FCI's reach), DMRG is the reference. This module provides both behind one interface so the
benchmark ladder can be extended in an environment that has `block2` installed.

  fci_energy(...)        exact active-space FCI via PySCF        (always available; tested here)
  dmrg_energy(...)       DMRG via block2 / pyblock2              (requires: pip install block2)
  reference_energy(...)  DMRG if available (method="auto"), else exact FCI

Integral convention (identical to ``molecular_hamiltonian.build_hamiltonian_from_integrals``):
  h1     : (norb, norb) one-electron MO integrals       = cas.get_h1eff()[0]
  eri    : (norb,)*4 chemist (pq|rs), 4-index restored  = ao2mo.restore(1, cas.get_h2eff(), norb)
  e_core : constant (nuclear + frozen-core energy)      = cas.get_h1eff()[1]
  n_elec : (n_alpha, n_beta) active electrons           = cas.nelecas

VALIDATED: with ``block2`` installed, ``dmrg_energy`` reproduces exact FCI to ~1e-10 Ha on small
active spaces (LiH CAS(4,5); N2 CAS up to (12,12) in ``benchmark_dmrg.py``), and carries the
reference alone past FCI's determinant reach (N2 CAS(14,14), ~1.2e7 determinants). The
``fci_energy`` fallback is exercised by the test suite, which also pins the shared integral
convention DMRG relies on.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np


def _as_pair(n_elec) -> Tuple[int, int]:
    if isinstance(n_elec, (int, np.integer)):
        n_beta = int(n_elec) // 2
        return (int(n_elec) - n_beta, n_beta)
    return (int(n_elec[0]), int(n_elec[1]))


def fci_energy(h1: np.ndarray, eri: np.ndarray, n_elec, e_core: float = 0.0) -> float:
    """Exact active-space FCI total energy via PySCF's direct solver on raw integrals."""
    from pyscf import fci

    norb = h1.shape[0]
    na, nb = _as_pair(n_elec)
    energy, _ = fci.direct_spin1.kernel(np.asarray(h1), np.asarray(eri), norb, (na, nb),
                                        ecore=float(e_core))
    return float(energy)


def dmrg_available() -> bool:
    try:
        import pyblock2  # noqa: F401
        return True
    except Exception:
        return False


def dmrg_energy(
    h1: np.ndarray,
    eri: np.ndarray,
    n_elec,
    e_core: float = 0.0,
    *,
    bond_dims=(100, 200, 400, 600),
    n_sweeps: int = 20,
    noises=(1e-4, 1e-5, 1e-6, 0.0),
    tol: float = 1e-9,
    scratch: str = "./.dmrg_tmp",
    n_threads: int = 4,
) -> float:
    """DMRG total energy via block2 / pyblock2. Requires ``pip install block2``.

    Intended for active spaces too large for exact FCI. On small spaces it reproduces
    ``fci_energy`` (the natural validation in an environment that has block2).
    """
    try:
        from pyblock2.driver.core import DMRGDriver, SymmetryTypes
    except Exception as exc:  # pragma: no cover - exercised only without block2
        raise ImportError(
            "dmrg_energy requires the 'block2' package (pip install block2). "
            "For active spaces where exact FCI is tractable, use fci_energy(...) instead."
        ) from exc

    norb = h1.shape[0]
    na, nb = _as_pair(n_elec)
    driver = DMRGDriver(scratch=scratch, symm_type=SymmetryTypes.SU2, n_threads=n_threads)
    driver.initialize_system(n_sites=norb, n_elec=na + nb, spin=na - nb, orb_sym=None)
    mpo = driver.get_qc_mpo(h1e=np.asarray(h1), g2e=np.asarray(eri), ecore=float(e_core), iprint=0)
    ket = driver.get_random_mps(tag="KET", bond_dim=int(bond_dims[0]), nroots=1)
    energy = driver.dmrg(
        mpo, ket,
        n_sweeps=n_sweeps,
        bond_dims=[int(b) for b in bond_dims],
        noises=list(noises),
        thrds=[tol] * len(bond_dims),
        iprint=0,
    )
    return float(energy)


@dataclass
class ExtrapResult:
    """Bond-dimension extrapolation of a DMRG energy (see specs/SPEC_hchain_tdl.md)."""
    energy: float                                   # E(D -> infinity)
    stderr: float                                   # standard error of the extrapolated energy
    per_D: List[Tuple[int, float, float]] = field(default_factory=list)  # (D, discarded_weight, E)
    method: str = "dweight"                         # "dweight" (E vs truncation error) or "invD"


def dmrg_energy_extrapolated(
    h1: np.ndarray,
    eri: np.ndarray,
    n_elec,
    e_core: float = 0.0,
    *,
    bond_dims=(200, 400, 800, 1600),
    n_sweeps_per: int = 8,
    protocol: str = "perD",
    sweeps_per_stage: int = 4,
    scratch: str = "./.dmrg_tmp",
    n_threads: int = 4,
    seed=None,
) -> ExtrapResult:
    """DMRG energy extrapolated to infinite bond dimension via the discarded-weight rule.

    The converged energy ``E(D)`` is linear in the discarded weight ``delta(D)`` near convergence,
    so a least-squares fit of ``E`` vs ``delta`` extrapolated to ``delta = 0`` gives ``E(D -> inf)``
    (White/Chan). Falls back to an ``E`` vs ``1/D`` fit (``method='invD'``) if the discarded weights
    are unusable (all ~0 or non-monotone).

    ``protocol``:
      * ``"perD"`` (default) -- a SEPARATE converged DMRG per bond dimension (warm-started, final
        sweeps at zero noise). Cleanest truncation points; ~``len(bond_dims)``x the sweeps.
      * ``"ramp"`` -- ONE DMRG run whose schedule holds each ``D`` for ``sweeps_per_stage`` sweeps;
        the per-stage points come from block2's ``get_dmrg_results()``. Much cheaper (see
        specs/SPEC_singleramp.md); intermediate-``D`` points are slightly less converged.
    """
    try:
        from pyblock2.driver.core import DMRGDriver, SymmetryTypes
    except Exception as exc:  # pragma: no cover - only without block2
        raise ImportError("dmrg_energy_extrapolated requires block2 (pip install block2).") from exc
    if protocol not in ("perD", "ramp"):
        raise ValueError(f"protocol must be 'perD' or 'ramp', got {protocol!r}")

    norb = h1.shape[0]
    na, nb = _as_pair(n_elec)
    drv = DMRGDriver(scratch=scratch, symm_type=SymmetryTypes.SU2, n_threads=n_threads)
    if seed is not None:
        try:
            import block2
            block2.Random.rand_seed(int(seed))   # seeds the global RNG used by get_random_mps
        except Exception:
            pass  # converged DMRG energy is independent of the random initial MPS anyway
    drv.initialize_system(n_sites=norb, n_elec=na + nb, spin=na - nb, orb_sym=None)
    mpo = drv.get_qc_mpo(h1e=np.asarray(h1), g2e=np.asarray(eri), ecore=float(e_core), iprint=0)
    ket = drv.get_random_mps(tag="KET", bond_dim=int(bond_dims[0]), nroots=1)

    per_D: List[Tuple[int, float, float]] = []
    if protocol == "ramp":
        # one run; schedule holds each target D for sweeps_per_stage sweeps
        schedule = [int(D) for D in bond_dims for _ in range(sweeps_per_stage)]
        nsw = len(schedule)
        noises = ([1e-4, 1e-5, 1e-6] + [0.0] * nsw)[:nsw]
        drv.dmrg(mpo, ket, n_sweeps=nsw, bond_dims=schedule, noises=noises,
                 thrds=[1e-12] * nsw, iprint=0)
        r_dims, r_dws, r_energies = drv.get_dmrg_results()   # per distinct bond dim
        by_D = {int(d): (float(w), float(np.asarray(e).ravel()[0]))
                for d, w, e in zip(r_dims, r_dws, r_energies)}
        per_D = [(int(D), by_D[int(D)][0], by_D[int(D)][1]) for D in bond_dims if int(D) in by_D]
    else:
        # anneal noise to zero so the final sweeps report a clean truncation error
        noises = ([1e-4, 1e-5, 1e-6] + [0.0] * max(0, n_sweeps_per - 3))[:n_sweeps_per]
        for D in bond_dims:
            e = drv.dmrg(
                mpo, ket, n_sweeps=n_sweeps_per, bond_dims=[int(D)] * n_sweeps_per,
                noises=noises, thrds=[1e-12] * n_sweeps_per, iprint=0,
            )
            dw = float(drv._dmrg.discarded_weights[-1])
            per_D.append((int(D), dw, float(e)))

    Ds = np.array([p[0] for p in per_D], dtype=float)
    dws = np.array([p[1] for p in per_D], dtype=float)
    Es = np.array([p[2] for p in per_D], dtype=float)

    method = "dweight"
    x = dws
    usable = (len(per_D) >= 2 and not np.allclose(dws, 0.0)
              and np.all(np.diff(dws) <= 1e-12))    # monotone non-increasing
    if not usable:
        x, method = 1.0 / Ds, "invD"

    if len(per_D) >= 3:
        coef, cov = np.polyfit(x, Es, 1, cov=True)
        energy = float(coef[1])
        stderr = float(np.sqrt(max(cov[1, 1], 0.0)))
    elif len(per_D) == 2:
        coef = np.polyfit(x, Es, 1)          # exact line through 2 points; no covariance
        energy, stderr = float(coef[1]), 0.0
    else:
        energy, stderr = float(Es[-1]), 0.0

    return ExtrapResult(energy=energy, stderr=stderr, per_D=per_D, method=method)


def thermodynamic_limit_fit(ns, e_per_atom) -> Tuple[float, float]:
    """Extrapolate per-atom energy to n -> infinity: E(n)/n ~ e_inf + a/n (open-chain surface term).

    Returns ``(e_inf, stderr)`` from a linear fit of ``e_per_atom`` vs ``1/n``.
    """
    x = 1.0 / np.asarray(ns, dtype=float)
    y = np.asarray(e_per_atom, dtype=float)
    if len(x) < 2:
        return float(y[-1]), 0.0
    if len(x) == 2:
        return float(np.polyfit(x, y, 1)[1]), 0.0
    coef, cov = np.polyfit(x, y, 1, cov=True)
    return float(coef[1]), float(np.sqrt(max(cov[1, 1], 0.0)))


def reference_energy(h1, eri, n_elec, e_core: float = 0.0, method: str = "auto", **kwargs):
    """Return ``(energy, method_used)``. ``method`` in {"auto", "fci", "dmrg"}.

    "auto" uses DMRG when block2 is importable, otherwise exact FCI.
    """
    if method == "fci":
        return fci_energy(h1, eri, n_elec, e_core), "fci"
    if method == "dmrg":
        return dmrg_energy(h1, eri, n_elec, e_core, **kwargs), "dmrg"
    if method == "auto":
        if dmrg_available():
            return dmrg_energy(h1, eri, n_elec, e_core, **kwargs), "dmrg"
        return fci_energy(h1, eri, n_elec, e_core), "fci"
    raise ValueError(f"unknown method {method!r}; expected 'auto', 'fci', or 'dmrg'")
