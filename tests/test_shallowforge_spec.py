"""Spec gate for ShallowForge M1 (specs/SPEC_shallowforge.md, tasks 1-4)."""

import json
import subprocess
import sys
from pathlib import Path

import pytest
from jsonschema import ValidationError

from shallowforge import (
    TermStreamIR,
    TransformEntry,
    baseline_cx_per_step,
    build_manifest,
    hamiltonian_hash,
    step_fidelity,
    validate_manifest,
)

_BASELINES = json.loads((Path(__file__).resolve().parents[1]
                         / "shallowforge" / "baselines.json").read_text())
_GOLDEN = {
    "T0": dict(atom="H 0 0 0; H 0 0 0.74"),
    "T1": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
    "T2": dict(atom="Li 0 0 0; H 0 0 1.6", active_electrons=2, active_orbitals=5),
    "T3": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
}


def _ir(name):
    from hybrid_quantum_solver import build_molecular_hamiltonian
    mh = build_molecular_hamiltonian(basis="sto3g", **_GOLDEN[name])
    return mh.qubit_hamiltonian, TermStreamIR.from_sparse_pauli_op(mh.qubit_hamiltonian)


# --- Gate 1: IR hash + round-trip -------------------------------------------------------


def test_ir_imports_without_solver_deps():
    code = (
        "import shallowforge.ir, sys;"
        "bad=[m for m in ('pyscf','qiskit','numpy','scipy') if m in sys.modules];"
        "assert not bad, bad; print('clean')"
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert "clean" in proc.stdout


def test_hash_invariant_under_input_reordering():
    terms = [("XX", 0.5 + 0j), ("ZI", -1.2 + 0j), ("IZ", 0.3 + 0j)]
    h1 = hamiltonian_hash(terms)
    h2 = hamiltonian_hash(list(reversed(terms)))
    assert h1 == h2
    ir = TermStreamIR(2, tuple(terms))
    assert ir.hamiltonian_hash == h1
    assert ir.reordered((2, 0, 1)).hamiltonian_hash == h1  # same Hamiltonian, same hash


def test_ir_roundtrips_dict_and_sparse_pauli_op():
    ir = TermStreamIR(2, (("XX", 0.5 + 0j), ("ZI", -1.2 + 0j)))
    assert TermStreamIR.from_dict(ir.to_dict()) == ir
    back = TermStreamIR.from_sparse_pauli_op(ir.to_sparse_pauli_op())
    assert back.hamiltonian_hash == ir.hamiltonian_hash


def test_hash_matches_chemcheck_for_same_operator():
    from chemcheck.tiers import canonical_hamiltonian_sha256

    op, ir = _ir("T0")

    class _MH:  # canonical_hamiltonian_sha256 only touches .qubit_hamiltonian
        qubit_hamiltonian = op

    assert ir.hamiltonian_hash == canonical_hamiltonian_sha256(_MH())


def test_ir_rejects_mismatched_label_length():
    with pytest.raises(ValueError):
        TermStreamIR(2, (("XXX", 1.0 + 0j),))


# --- Gate 2: manifest emitter -----------------------------------------------------------


def test_manifest_validates_against_schema():
    m = build_manifest(
        hamiltonian_hash="deadbeef",
        stack=[
            TransformEntry("term_ordering", {"key": "magnitude"}, 0.0, lossless=True),
            TransformEntry("trotter_order", {"order": 2}, 0.4, lossless=False),
        ],
        cx_per_step=3268, depth=900, ancillas=0, solver_version="0.1.0",
    )
    validate_manifest(m)  # raises on invalid
    assert m["totals"]["predicted_epsilon_total_mha"] == pytest.approx(0.4)
    assert m["totals"]["cx_at_epsilon_claim"]["epsilon_mha"] == 1.6
    assert m["totals"]["cx_at_epsilon_claim"]["cx_per_step"] == 3268


def test_lossless_transform_with_nonzero_epsilon_rejected():
    with pytest.raises(ValueError):
        TransformEntry("pauli_grouping", {}, 0.5, lossless=True)


def test_unknown_transform_rejected():
    with pytest.raises(ValueError):
        TransformEntry("magic_wand", {}, 0.0, lossless=True)


def test_malformed_manifest_fails_schema():
    with pytest.raises(ValidationError):
        validate_manifest({"manifest_version": "x"})  # missing required fields


# --- Gate 4: frozen baselines -----------------------------------------------------------


@pytest.mark.parametrize("name", ["T0", "T1", "T2", "T3"])
def test_baseline_cx_reproduced(name):
    _, ir = _ir(name)
    cx = baseline_cx_per_step(ir, _BASELINES["reference_dt"])
    assert cx == _BASELINES["systems"][name]["baseline_cx_per_step"]


def test_n2_second_order_headline_recorded():
    # The honest finding: the published "~6,500" is second-order, not the first-order baseline.
    t3 = _BASELINES["systems"]["T3"]
    assert t3["baseline_cx_per_step"] == 3268
    assert t3["second_order_reps1_cx_per_step"] == 6534


# --- Gate 3: correctness verifier (fidelity) --------------------------------------------


@pytest.mark.parametrize("name", ["T0", "T1", "T2"])
def test_step_fidelity_resolves_trotter_error(name):
    _, ir = _ir(name)
    dt = _BASELINES["reference_dt"]
    # First-order at the reference step is genuinely APPROXIMATE (infidelity > 1e-8) — this
    # verifier is not exact-evolution-in-disguise (the failure mode build_trotter_step warns of).
    f_first = step_fidelity(ir, dt, order=1, reps=1)
    assert f_first < 1.0 - 1e-8, f"{name}: first-order fidelity {f_first} looks exact"
    # The reference protocol (second-order Suzuki, 2 reps) clears the correctness gate.
    f_ref = step_fidelity(ir, dt, order=2, reps=2)
    assert f_ref >= 1.0 - 1e-8, f"{name}: reference-protocol fidelity {f_ref} below gate"
    assert f_ref > f_first
