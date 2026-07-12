"""Validated-envelope caps for CertChem-M1.

``check_caps`` is the single gate that refuses systems outside the envelope the solver core
is validated on, raising :class:`CapExceededError` that *names* the violated cap. Zero solver
dependencies — this is pure input validation.
"""

from __future__ import annotations

from typing import Any

from .contract import CapExceededError

#: Hard cap on spin-orbitals (= 2 x active spatial orbitals). Beyond this the ODMD/Krylov
#: path is not validated against exact references, so we refuse rather than guess.
MAX_SPIN_ORBITALS = 16

#: Bases the solver core is validated on (see benchmark_krylov.py / benchmark_n2.py).
#: Compared case-insensitively.
ALLOWED_BASES = frozenset({"sto-3g", "6-31g", "cc-pvdz"})


def check_caps(molecule: Any, basis: str, cas: tuple[int, int]) -> None:
    """Validate ``(molecule, basis, cas)`` against the solver's envelope.

    ``cas`` is ``(n_electrons, n_active_orbitals)``. Raises :class:`CapExceededError` naming
    the first violated cap; returns ``None`` when every cap passes.
    """
    if molecule is None:
        raise CapExceededError("molecule is None", cap="molecule")

    if basis is None or basis.lower() not in ALLOWED_BASES:
        raise CapExceededError(
            f"basis {basis!r} not in validated set {sorted(ALLOWED_BASES)}",
            cap="basis",
        )

    n_electrons, n_orbitals = cas
    if n_orbitals <= 0:
        raise CapExceededError(
            f"active orbitals must be positive, got {n_orbitals}", cap="active_orbitals"
        )

    spin_orbitals = 2 * n_orbitals
    if spin_orbitals > MAX_SPIN_ORBITALS:
        raise CapExceededError(
            f"{spin_orbitals} spin-orbitals exceeds cap of {MAX_SPIN_ORBITALS} "
            f"(active orbitals {n_orbitals} > {MAX_SPIN_ORBITALS // 2})",
            cap="spin_orbitals",
        )

    if n_electrons < 0 or n_electrons > spin_orbitals:
        raise CapExceededError(
            f"{n_electrons} electrons cannot occupy {spin_orbitals} spin-orbitals",
            cap="electron_count",
        )
