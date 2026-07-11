"""Term-stream IR (task 1): a Hamiltonian as an ordered stream of Pauli terms.

Every compiler rung is a pure function ``TermStreamIR -> TermStreamIR``. The IR carries a
content hash (:func:`hamiltonian_hash`) that is invariant under reordering of the *input* terms,
so the audit trail identifies *which Hamiltonian* a circuit implements regardless of how the
terms happened to be listed. The hash is byte-compatible with ``chemcheck.canonical_hamiltonian_
sha256`` — the same operator gets the same ``hamiltonian_hash`` in both packages.

Zero solver deps at import; qiskit is imported lazily only in the SparsePauliOp bridges.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from typing import Any

#: Coefficients are quantized to this many decimals before hashing/serialization (determinism).
_COEFF_DIGITS = 12


def _canonical_terms(terms: list[tuple[str, complex]]) -> list[tuple[str, float, float]]:
    """Sorted, coefficient-quantized ``(label, Re, Im)`` stream — the canonical form."""
    return sorted(
        (label, round(float(c.real), _COEFF_DIGITS), round(float(c.imag), _COEFF_DIGITS))
        for label, c in terms
    )


def hamiltonian_hash(terms: list[tuple[str, complex]]) -> str:
    """SHA-256 of the canonical term stream. Invariant under input reordering."""
    canon = _canonical_terms(terms)
    return hashlib.sha256(json.dumps(canon, separators=(",", ":")).encode()).hexdigest()


@dataclass(frozen=True)
class TermStreamIR:
    """An immutable ordered stream of Pauli terms plus qubit count.

    ``terms`` preserves the order a transform produced (ordering *is* a compiler lever); the hash
    is computed on the canonical (sorted) form so it does not depend on that order.
    """

    num_qubits: int
    terms: tuple[tuple[str, complex], ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        for label, _ in self.terms:
            if len(label) != self.num_qubits:
                raise ValueError(
                    f"Pauli label {label!r} has length {len(label)} != num_qubits "
                    f"{self.num_qubits}"
                )

    @property
    def hamiltonian_hash(self) -> str:
        """Content hash of the Hamiltonian this IR represents (order-independent)."""
        return hamiltonian_hash(list(self.terms))

    def reordered(self, order: tuple[int, ...]) -> TermStreamIR:
        """Return a copy with terms permuted by ``order`` (same Hamiltonian, same hash)."""
        return TermStreamIR(self.num_qubits, tuple(self.terms[i] for i in order))

    # --- SparsePauliOp bridges (lazy qiskit import) --------------------------------------

    @classmethod
    def from_sparse_pauli_op(cls, op: Any) -> TermStreamIR:
        terms = tuple(
            (p.to_label(), complex(c)) for p, c in zip(op.paulis, op.coeffs)
        )
        return cls(num_qubits=op.num_qubits, terms=terms)

    def to_sparse_pauli_op(self) -> Any:
        from qiskit.quantum_info import SparsePauliOp

        labels = [label for label, _ in self.terms]
        coeffs = [c for _, c in self.terms]
        return SparsePauliOp(labels, coeffs)

    # --- round-trip serialization --------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "num_qubits": self.num_qubits,
            "terms": [[label, c.real, c.imag] for label, c in self.terms],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> TermStreamIR:
        terms = tuple((label, complex(re, im)) for label, re, im in data["terms"])
        return cls(num_qubits=int(data["num_qubits"]), terms=terms)
