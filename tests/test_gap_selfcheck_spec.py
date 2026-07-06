"""
Acceptance gates G1-G4 for specs/SPEC_gap_selfcheck.md.

Claim: the certified gap bracket can self-verify WITHOUT an oracle -- a bracket is corroborated iff
it overlaps the deep-anchor (intersection of the deepest brackets), and this ADAPTIVELY catches the
premise-failure regime (rejects M=4 on H4/N2 where certified_gaps escapes, accepts M=4 on LiH where
the premise already holds), yielding an oracle-free gap interval validated to contain the exact
reachable gap and to exclude the shallow outliers.

Exact statevector, sector-restricted. PySCF/qiskit, no block2; `make gates` runs it in its own
process.
"""
import numpy as np
import pytest

from certified_gaps import gap_bracket_ladder, reachable_gap
from gap_selfcheck import (
    anchor_interval,
    corroborated_flags,
    self_checked_gap,
)
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

CASES = {
    "H4": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
    "LiH": dict(atom="Li 0 0 0; H 0 0 1.6", active_electrons=2, active_orbitals=5),
    "N2": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
}
_DIMS = (4, 6, 8, 12, 16, 20, 24)


@pytest.fixture(scope="module")
def ladders():
    out = {}
    for name, spec in CASES.items():
        mh = build_molecular_hamiltonian(**spec)
        solver = QuantumKrylovSolver(mh)
        out[name] = (reachable_gap(mh), gap_bracket_ladder(mh, _DIMS, solver=solver))
    return out


def test_G1_adaptive_rejection_of_premise_failures(ladders):
    """The certificate is ADAPTIVE: it rejects the M=4 bracket for the multireference cases (H4, N2,
    where the lower certificate escapes) but ACCEPTS M=4 for LiH (premise already holds), and
    accepts EVERY M >= 6 for all systems. Oracle-free -- uses only cross-depth consistency."""
    for name, (_, brackets) in ladders.items():
        flags = dict(zip(_DIMS, corroborated_flags(brackets)))
        for m in (6, 8, 12, 16, 20, 24):
            assert flags[m], (name, m)                       # every deep bracket corroborated
    h4 = dict(zip(_DIMS, corroborated_flags(ladders["H4"][1])))
    n2 = dict(zip(_DIMS, corroborated_flags(ladders["N2"][1])))
    lih = dict(zip(_DIMS, corroborated_flags(ladders["LiH"][1])))
    assert h4[4] is False and n2[4] is False                 # premise failures caught
    assert lih[4] is True                                    # shallow-but-valid NOT falsely rejected


def test_G2_oracle_free_interval_covers_the_truth(ladders):
    """DEFINITION OF DONE: the self-checked interval (intersection of corroborated brackets) is
    non-empty, finite, and CONTAINS the exact reachable gap for every system -- an oracle-free error
    bar validated to cover the truth."""
    for name, (gap, brackets) in ladders.items():
        lo, hi = self_checked_gap(brackets)
        assert np.isfinite(lo) and np.isfinite(hi) and lo <= hi, (name, lo, hi)
        assert lo - 1e-9 <= gap <= hi + 1e-9, (name, lo, gap, hi)


def test_G3_self_check_repairs_the_naive_estimate(ladders):
    """The self-check strictly helps: naively intersecting ALL brackets (including the M=4 outlier)
    gives an EMPTY interval for the premise-failure cases, so it is useless; dropping the
    uncorroborated M=4 restores a non-empty interval that contains the gap."""
    for name in ("H4", "N2"):
        gap, brackets = ladders[name]
        naive_lo = max(b.gap_lower for b in brackets)        # intersection of all (incl. M=4)
        naive_hi = min(b.gap_upper for b in brackets)
        assert naive_lo > naive_hi, (name, naive_lo, naive_hi)   # naive all-intersection is EMPTY
        lo, hi = self_checked_gap(brackets)                  # corroborated-only
        assert lo <= hi and lo - 1e-9 <= gap <= hi + 1e-9, (name, lo, gap, hi)


def test_G4_anchor_robustness_and_minimal_premise(ladders):
    """Robustness + the honest boundary. The corroboration is stable to the anchor depth: using the
    deepest 2 vs deepest 3 brackets gives identical flags. And the minimal premise the certificate
    rests on -- that the deep anchor is itself self-consistent (non-empty intersection) -- holds;
    that anchor is taken on trust (necessary, not sufficient; see spec)."""
    for name, (_, brackets) in ladders.items():
        assert corroborated_flags(brackets, k=2) == corroborated_flags(brackets, k=3), name
        lo, hi = anchor_interval(brackets, k=2)
        assert lo <= hi, (name, lo, hi)                      # deep anchor self-consistent
        lo3, hi3 = anchor_interval(brackets, k=3)
        assert lo3 <= hi3, (name, lo3, hi3)
