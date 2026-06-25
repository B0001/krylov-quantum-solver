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

NOTE: the ``dmrg_energy`` path is written against the standard pyblock2 ``DMRGDriver`` API but was
not executed in the environment where this module was written (``block2`` could not be installed
there). The ``fci_energy`` path it falls back to is exercised by the test suite, which also pins
the shared integral convention DMRG relies on.
"""
from __future__ import annotations

from typing import Tuple

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
