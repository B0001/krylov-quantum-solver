#!/usr/bin/env python3
"""
model_hamiltonians.py -- lattice/cRPA model parameters -> universal active-space integrals.

This is the loader the session handoff (CLAUDE_CODE_HANDOFF.md) asks for: a bridge that maps a
tight-binding hopping matrix plus a Hubbard/cRPA interaction into the project's universal
``(h1, eri, e_core, nelec, norb)`` tuple, so the *same* validated solver stack
(``build_hamiltonian_from_integrals`` -> ``QuantumKrylovSolver`` / ``SampleKrylovSolver``) runs on
model Hamiltonians without going through heavy-atom PySCF/CASCI.

Why this matters for Nb3X8
--------------------------
The Nb3X8 Model Database (Aretz, Grytsiuk, Strand, van Loon, Rosner 2025; ref [86] of
``arXiv:2501.10320``, PRX DOI 10.1103/wr7w-nfhg) publishes, for every compound x structure, the
single-particle hopping matrix and the full rank-4 cRPA Coulomb tensor U_ijkl in a molecular-orbital
basis. Those map *directly* onto our interface: ``h1 <- hopping``, ``eri <- U_ijkl`` (chemist
notation). That lets the stack be validated against expert DMFT/Hubbard-I results with no open-shell
Nb SCF at all.

The bulk low-temperature electronic structure of Nb3X8 collapses (under breathing-mode
trimerization) to **one flat band per monolayer**, half-filled; bulk = two dimerized trimer MOs --
a generalized **Hubbard dimer** (bonding/antibonding inter-layer pair), singlet ground state. So the
smallest faithful model is the two-site, two-electron Hubbard dimer, which has an exact analytic
ground-state energy -- the falsifiable reference this module is validated against.

Honest scope
------------
* The cRPA/DB path (``hubbard_integrals`` with a full rank-4 tensor, or ``load_from_nb3x8_database``
  with parsed ``hopping``/``coulomb`` arrays) is exercised against PySCF FCI in the matching
  particle-number sector -- a real ground-truth check of the *mapping*.
* The named-compound shortcut builds the published bulk-LT Hubbard *dimer* from
  :data:`NB3X8_BULK_DIMER_PARAMS`. It is validated only as the analytic Hubbard dimer. The
  database's own DMFT/Hubbard-I **gap** is not bundled here (the DB file format is not in the repo;
  contact m.roesner@science.ru.nl) -- so a compound's reported gap is a *model prediction*, not yet a
  validated-against-DMFT number. Do not present it as the latter.

A subtlety that is itself a finding: a Hubbard model lives at a *fixed* filling.
``MolecularHamiltonian.ground_state_energy()`` diagonalizes the full Fock space and, for U/t above
~3, returns the *one-electron* bonding state (E = -t) instead of the half-filled singlet. Use the
particle-number-conserving solvers (``QuantumKrylovSolver`` / ``SampleKrylovSolver``) or
:func:`fixed_filling_energy` -- never the bare full-space diagonalization -- for a fixed-filling
model. This is the concrete content of SKQD checklist item 5 (U(1)/electron-number conservation).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Tuple, Union

import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import (
    MolecularHamiltonian,
    build_hamiltonian_from_integrals,
)

# 2018 CODATA Hartree <-> eV, matching PySCF's internal constant.
HARTREE_PER_EV = 1.0 / 27.211386245988
_UNIT_TO_HA = {"Ha": 1.0, "hartree": 1.0, "eV": HARTREE_PER_EV, "meV": HARTREE_PER_EV / 1000.0}

Interaction = Union[float, np.ndarray]


@dataclass
class ModelIntegrals:
    """The universal active-space tuple, plus convenience constructors.

    ``(h1, eri, e_core, nelec, norb)`` is the interface shared across the pipeline; ``h1``/``eri``/
    ``e_core`` are in Hartree (converted on load), ``nelec`` is ``(n_alpha, n_beta)``.
    """

    h1: np.ndarray
    eri: np.ndarray
    e_core: float
    nelec: Tuple[int, int]
    norb: int

    def as_tuple(self) -> Tuple[np.ndarray, np.ndarray, float, Tuple[int, int], int]:
        return (self.h1, self.eri, self.e_core, self.nelec, self.norb)

    def to_hamiltonian(self) -> MolecularHamiltonian:
        """Build the qubit Hamiltonian via the vetted Jordan-Wigner path."""
        return build_hamiltonian_from_integrals(self.h1, self.eri, self.nelec, self.e_core)


def _split_nelec(nelec: Union[int, Tuple[int, int]]) -> Tuple[int, int]:
    """Normalize an electron count to ``(n_alpha, n_beta)`` (extra electron -> alpha)."""
    if isinstance(nelec, (tuple, list)):
        if len(nelec) != 2:
            raise ValueError("nelec tuple must be (n_alpha, n_beta)")
        return int(nelec[0]), int(nelec[1])
    n = int(nelec)
    if n < 0:
        raise ValueError("nelec must be non-negative")
    return (n - n // 2, n // 2)


def hubbard_integrals(
    hopping: np.ndarray,
    interaction: Interaction,
    nelec: Union[int, Tuple[int, int]],
    *,
    e_core: float = 0.0,
    units: str = "Ha",
) -> ModelIntegrals:
    """Map a tight-binding hopping matrix + a Hubbard/cRPA interaction into the universal tuple.

    Args:
        hopping: ``(norb, norb)`` one-body matrix, used as ``h1`` verbatim -- put on-site energies on
            the diagonal and hopping amplitudes (e.g. ``-t``) off-diagonal, with whatever sign
            convention you intend. Must be symmetric (Hermitian, real).
        interaction: the two-body term, in one of three forms:
            * scalar ``U`` -- on-site Hubbard ``U`` applied to every site;
            * length-``norb`` vector -- per-site ``U_i``;
            * rank-4 ``(norb,)*4`` tensor -- a full Coulomb tensor ``U_ijkl`` in **chemist notation**
              ``(ij|kl)`` (the cRPA / Nb3X8-database case), used as ``eri`` verbatim.
        nelec: total electron count, or an explicit ``(n_alpha, n_beta)`` pair, fixing the filling.
        e_core: constant energy offset (Hartree-frame; converted with ``units``).
        units: ``"Ha"`` (default), ``"eV"``, or ``"meV"`` -- all inputs are scaled to Hartree so the
            solver, which assumes atomic units, stays consistent.

    Returns:
        A :class:`ModelIntegrals` ready for ``.to_hamiltonian()`` or ``.as_tuple()``.

    Note:
        A Hubbard model is solved at *fixed filling*. Diagonalize with the number-conserving solvers
        or :func:`fixed_filling_energy`, not the full-Fock-space ``ground_state_energy()`` -- see the
        module docstring.
    """
    try:
        scale = _UNIT_TO_HA[units]
    except KeyError:
        raise ValueError(f"units must be one of {sorted(_UNIT_TO_HA)}, got {units!r}")

    h1 = np.asarray(hopping, dtype=float)
    if h1.ndim != 2 or h1.shape[0] != h1.shape[1]:
        raise ValueError(f"hopping must be square (norb, norb), got {h1.shape}")
    norb = h1.shape[0]
    if not np.allclose(h1, h1.T, atol=1e-12):
        raise ValueError("hopping must be symmetric")
    h1 = h1 * scale

    interaction = np.asarray(interaction, dtype=float)
    if interaction.ndim == 0:                                   # scalar U -> on-site, every site
        eri = np.zeros((norb,) * 4)
        idx = np.arange(norb)
        eri[idx, idx, idx, idx] = float(interaction)
    elif interaction.ndim == 1:                                 # per-site U_i
        if interaction.shape[0] != norb:
            raise ValueError(f"per-site interaction must have length {norb}, got {interaction.shape}")
        eri = np.zeros((norb,) * 4)
        idx = np.arange(norb)
        eri[idx, idx, idx, idx] = interaction
    elif interaction.ndim == 4:                                 # full rank-4 Coulomb tensor (cRPA)
        if interaction.shape != (norb,) * 4:
            raise ValueError(f"rank-4 interaction must have shape {(norb,) * 4}, got {interaction.shape}")
        eri = interaction.copy()
    else:
        raise ValueError("interaction must be scalar, length-norb vector, or rank-4 tensor")
    eri = eri * scale

    return ModelIntegrals(
        h1=h1,
        eri=eri,
        e_core=float(e_core) * scale,
        nelec=_split_nelec(nelec),
        norb=norb,
    )


def fixed_filling_energy(model: ModelIntegrals) -> float:
    """Exact ground-state energy of ``model`` **in its fixed particle-number sector** (PySCF FCI).

    This is the correct reference for a lattice model: unlike the full-Fock-space
    ``MolecularHamiltonian.ground_state_energy()``, it stays at the requested electron number.
    """
    from hybrid_quantum_solver.dmrg_reference import fci_energy

    return fci_energy(model.h1, model.eri, model.nelec, model.e_core)


# -- analytic references --------------------------------------------------------------------------

def hubbard_dimer_energy(t: float, U: float) -> float:
    """Exact ground-state energy of the half-filled (2-electron) two-site Hubbard model.

    ``E0 = (U - sqrt(U**2 + 16 t**2)) / 2`` (singlet). Limits: ``U=0 -> -2|t|`` (both electrons in the
    bonding orbital); ``U -> inf -> -4 t**2 / U`` (Heisenberg superexchange ``-J``). The triplet sits
    at 0, so this energy is also minus the singlet-triplet gap. Returned in the same units as ``t, U``.
    """
    return 0.5 * (U - math.sqrt(U * U + 16.0 * t * t))


def hubbard_dimer_gap(t: float, U: float) -> float:
    """Singlet-triplet gap of the half-filled two-site Hubbard model (triplet at 0)."""
    return -hubbard_dimer_energy(t, U)


def hubbard_chain_integrals(
    n_sites: int,
    U: float,
    t: float = 1.0,
    *,
    closed_shell: bool = True,
    e_core: float = 0.0,
    units: str = "Ha",
) -> ModelIntegrals:
    """Half-filled 1D Hubbard ring: nearest-neighbour hopping ``t``, on-site ``U``, as integrals.

    ``h1`` is the ``n_sites x n_sites`` ring hopping matrix (``-t`` on the nearest-neighbour
    off-diagonals and the wrap). ``closed_shell`` picks the boundary phase that closes the
    non-interacting shell -- **periodic** when ``n_sites/2`` is odd, **antiperiodic** (``+t`` on the
    wrap) when even. That removes the open-shell even/odd-L zigzag in the finite-size energy, so the
    per-site energy converges smoothly (``~1/L^2``) to the thermodynamic limit and can be compared
    to :func:`lieb_wu_energy`. Filling is fixed to half (``n_sites/2`` electrons per spin); ``n_sites``
    must be even.
    """
    if n_sites % 2:
        raise ValueError("n_sites must be even (half-filling)")
    h = np.zeros((n_sites, n_sites))
    for i in range(n_sites - 1):
        h[i, i + 1] = h[i + 1, i] = -t
    if n_sites > 2:
        wrap = -t if (not closed_shell or (n_sites // 2) % 2 == 1) else +t
        h[0, n_sites - 1] = h[n_sites - 1, 0] = wrap
    return hubbard_integrals(h, U, nelec=(n_sites // 2, n_sites // 2), e_core=e_core, units=units)


def lieb_wu_energy(U: float, t: float = 1.0) -> float:
    """Exact half-filled 1D Hubbard ground-state energy **per site** in the thermodynamic limit.

    The Bethe-ansatz result (Lieb & Wu, Phys. Rev. Lett. 20, 1445, 1968):

        e0 / t = -4 * integral_0^inf  J0(w) J1(w) / ( w (1 + exp(w U / (2 t))) )  dw

    with ``J0, J1`` Bessel functions. Limits: ``U=0 -> -4 t / pi`` (free fermions); ``U -> inf -> 0``
    (Mott insulator, ``~ -4 ln2 t^2 / U`` superexchange). This is the analytic reference the finite-
    size :func:`hubbard_chain_integrals` energies are extrapolated against (see specs/SPEC_hubbard_bethe.md).
    """
    from scipy.integrate import quad
    from scipy.special import j0, j1

    def integrand(w):
        return j0(w) * j1(w) / (w * (1.0 + np.exp(np.clip(w * U / (2.0 * t), 0.0, 700.0))))

    value, _ = quad(integrand, 0.0, np.inf, limit=400)
    return -4.0 * t * value


# -- Nb3X8 -------------------------------------------------------------------------------------

# Bulk low-temperature generalized Hubbard *dimer* parameters (meV). The bulk LT manifold is two
# dimerized trimer MOs -> a 2-site, 2-electron Hubbard dimer whose strong bond is the inter-layer
# hopping t_perp_strong. Source: arXiv:2501.10320 (PRX 10.1103/wr7w-nfhg), Table of bulk parameters.
NB3X8_BULK_DIMER_PARAMS = {
    "Nb3I8": {
        "U": 787.0,                # on-site (per-trimer-MO) Hubbard U  [meV]
        "t": 218.2,                # intra-dimer (inter-layer strong) hopping |t_perp^s|  [meV]
        "t_intralayer": 8.3,       # in-plane hopping t_parallel  [meV] (not in the 2-site dimer)
        "t_perp_weak": 24.6,       # weak inter-layer hopping  [meV]
        "correlation": "weak (U/t ~ 3.6; obstructed-atomic band insulator)",
        "source": "arXiv:2501.10320 / PRX 10.1103/wr7w-nfhg",
    },
}


def load_from_nb3x8_database(
    compound: str = "Nb3I8",
    *,
    hopping: Optional[np.ndarray] = None,
    coulomb: Optional[Interaction] = None,
    nelec: Optional[Union[int, Tuple[int, int]]] = None,
    e_core: float = 0.0,
    units: str = "meV",
) -> ModelIntegrals:
    """Map Nb3X8 model parameters into the universal active-space tuple.

    Two modes:

    * **Database mode** -- pass ``hopping`` (the single-particle matrix) and ``coulomb`` (a scalar,
      per-orbital vector, or the full rank-4 cRPA tensor ``U_ijkl``), as parsed from the Nb3X8 Model
      Database (ref [86] of ``arXiv:2501.10320``). This is the intended drop-in: ``h1 <- hopping``,
      ``eri <- coulomb``. ``nelec`` is required (the database fixes the filling, e.g. one electron
      per trimer). The mapping is validated against PySCF FCI; see the spec.

    * **Named-compound mode** -- pass only ``compound`` to build the published bulk low-T generalized
      Hubbard **dimer** from :data:`NB3X8_BULK_DIMER_PARAMS`: a 2-MO, 2-electron model with
      ``h1 = [[0, -t], [-t, 0]]`` and on-site ``U``. Half-filled (singlet ground state). This path is
      validated only as the analytic Hubbard dimer -- the DB's DMFT/Hubbard-I gap is not bundled.

    Args:
        compound: key into :data:`NB3X8_BULK_DIMER_PARAMS` (named-compound mode).
        hopping, coulomb, nelec: database-mode inputs (override the named-compound shortcut).
        e_core: constant offset in ``units``.
        units: units of all inputs (default ``"meV"`` to match the database / the published tables).

    Returns:
        A :class:`ModelIntegrals` in Hartree.
    """
    if hopping is not None or coulomb is not None:
        if hopping is None or coulomb is None or nelec is None:
            raise ValueError("database mode needs hopping, coulomb, and nelec together")
        return hubbard_integrals(hopping, coulomb, nelec, e_core=e_core, units=units)

    if compound not in NB3X8_BULK_DIMER_PARAMS:
        raise KeyError(
            f"unknown compound {compound!r}; known: {sorted(NB3X8_BULK_DIMER_PARAMS)}. "
            "For other compounds/structures, pass hopping + coulomb parsed from the Nb3X8 database."
        )
    p = NB3X8_BULK_DIMER_PARAMS[compound]
    t, U = p["t"], p["U"]
    hop = np.array([[0.0, -t], [-t, 0.0]])
    return hubbard_integrals(hop, U, nelec=(1, 1), e_core=e_core, units=units)


if __name__ == "__main__":
    from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

    print("Half-filled Hubbard dimer  (Krylov vs analytic, t=1):")
    for U in (0.0, 2.0, 4.0, 8.0, 20.0):
        model = hubbard_integrals(np.array([[0.0, -1.0], [-1.0, 0.0]]), U, nelec=2)
        e_kry = QuantumKrylovSolver(model.to_hamiltonian()).solve(8).energy
        print(f"  U={U:5.1f}  E_krylov={e_kry:+.6f}  analytic={hubbard_dimer_energy(1.0, U):+.6f}")

    print("\nNb3I8 bulk-LT Hubbard dimer (meV), from the published parameters:")
    p = NB3X8_BULK_DIMER_PARAMS["Nb3I8"]
    print(f"  U={p['U']} meV, t={p['t']} meV  ->  singlet-triplet gap = "
          f"{hubbard_dimer_gap(p['t'], p['U']):.1f} meV   ({p['correlation']})")
