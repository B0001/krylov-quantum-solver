"""Gates for specs/SPEC_scf_conv_tol.md (chem-ayr): does the SCF stopping tolerance move any number
the certified arc consumes?

Every threshold below is transcribed from the spec's section 5, which was committed BEFORE this file
existed. Both tolerances are built in this one process (the default through the builder's own default,
the tight one through `conv_tol=TIGHT_SCF_CONV_TOL`) so run-to-run noise cannot confound the diff.

Regenerate the table recorded in the spec's results section (`t.residue_sweep()` prints the a=1.10 sweep):
    uv run --no-sync python -c "import sys; sys.path.insert(0, 'tests'); \
import test_scf_conv_tol_spec as t; t.table()"
"""
import inspect
from functools import lru_cache

import numpy as np
import pytest

from hf_overlap_certificate import exact_reachable_overlap
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from reachability import (
    REACHABLE_TOL_CERTIFIED,
    TIGHT_SCF_CONV_TOL,
    _dense_hf_projection,
    reachable_eigenpairs,
)

LEAK_TOL = 1e-12        # == test_chained_overlap_spec.LEAK_TOL (test modules are not importable)
LEAK_DEPTH = 6


def _square(a):
    return f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0"


def _chain(d):
    return f"H 0 0 0; H 0 0 {d}; H 0 0 {2 * d}; H 0 0 {3 * d}"


# Clean set: 8 qubits and below, artifact absent at the default. Not gated: square a=1.10 (its
# default-tolerance residue is platform-dependent), linear H6 (12 qubits), LiH/N2 CAS tiers.
CLEAN = {
    "H2 0.74": "H 0 0 0; H 0 0 0.74",
    "H2 2.0": "H 0 0 0; H 0 0 2.0",
    "H4 chain 0.9": _chain(0.9),
    "H4 chain 1.0": _chain(1.0),
    "H4 chain 2.0": _chain(2.0),
    "H4 square 1.0": _square(1.0),
    "H4 square 1.05": _square(1.05),
    "H4 square 1.2": _square(1.2),
    "H4 square 1.4": _square(1.4),
}
# Symmetric-SCF square H4 where the default-tolerance residue is large enough to be platform-robust.
WITNESS = {"H4 square 1.35": _square(1.35), "H4 square 1.19": _square(1.19)}


@lru_cache(maxsize=None)
def _build(atom, conv_tol):
    """conv_tol=None means: do not pass it, i.e. the builder's own default."""
    if conv_tol is None:
        return build_molecular_hamiltonian(atom=atom)
    return build_molecular_hamiltonian(atom=atom, conv_tol=conv_tol)


@lru_cache(maxsize=None)
def _measure(atom, conv_tol):
    mh = _build(atom, conv_tol)
    w, _, pops = _dense_hf_projection(mh)
    return {
        "E0": mh.ground_state_energy(),
        "EHF": mh.hf_energy,
        "ov": exact_reachable_overlap(mh),
        "E0reach": float(reachable_eigenpairs(mh)[0][0]),
        "p0": float(pops[0]),
        "i_reach": int(np.where(pops > REACHABLE_TOL_CERTIFIED)[0][0]),
        "n_reach": int(np.sum(pops > REACHABLE_TOL_CERTIFIED)),
    }


@lru_cache(maxsize=None)
def _leak(atom, conv_tol):
    """||P_unreach v|| for the depth-LEAK_DEPTH ground Ritz vector (chained-overlap R2b quantity)."""
    mh = _build(atom, conv_tol)
    _, vecs, pops = _dense_hf_projection(mh)
    v = np.asarray(QuantumKrylovSolver(mh).eigenstates(LEAK_DEPTH, n_states=1)[1][0], dtype=complex)
    v = v / np.linalg.norm(v)
    return float(np.linalg.norm(vecs[:, pops <= REACHABLE_TOL_CERTIFIED].conj().T @ v))


def _pair(atom):
    return _measure(atom, None), _measure(atom, TIGHT_SCF_CONV_TOL)


# --- G1: the clean set does not move (claim a) ----------------------------------------------------

@pytest.mark.parametrize("name", list(CLEAN))
def test_G1_clean_geometry_does_not_move(name):
    d, t = _pair(CLEAN[name])
    assert abs(d["E0"] - t["E0"]) < 1e-9, (name, d["E0"], t["E0"])
    assert abs(d["EHF"] - t["EHF"]) < 1e-8, (name, d["EHF"], t["EHF"])
    assert abs(d["ov"] - t["ov"]) < 1e-4, (name, d["ov"], t["ov"])
    assert d["i_reach"] == t["i_reach"], (name, d["i_reach"], t["i_reach"])


# --- G2: the witnesses move, the eigenvalue does not (claim b; DEFINITION OF DONE) -----------------

@pytest.mark.parametrize("name", list(WITNESS))
def test_G2_witness_moves_but_the_eigenvalue_does_not(name):
    d, t = _pair(WITNESS[name])
    assert abs(d["E0"] - t["E0"]) < 1e-9, (name, d["E0"], t["E0"])           # same eigenvalue
    assert d["p0"] > REACHABLE_TOL_CERTIFIED, (name, d["p0"])                # admitted as reachable
    assert d["ov"] < 1e-3, (name, d["ov"])                                   # ...at residue level
    assert abs(d["ov"] - t["ov"]) > 0.1, (name, d["ov"], t["ov"])
    assert abs(d["E0reach"] - t["E0reach"]) > 0.010, (name, d["E0reach"], t["E0reach"])


# --- G3: the residue collapses at TIGHT_SCF_CONV_TOL (claim c) -------------------------------------

@pytest.mark.parametrize("name", list(WITNESS))
def test_G3_residue_collapses_at_the_tight_tolerance(name):
    atom = WITNESS[name]
    assert _measure(atom, TIGHT_SCF_CONV_TOL)["p0"] < 1e-20, name
    assert _leak(atom, TIGHT_SCF_CONV_TOL) < LEAK_TOL, name
    assert _leak(atom, None) > LEAK_TOL, (name, "exclusion no longer warranted at the default")


# --- G4: the default is pinned ---------------------------------------------------------------------

def test_G4_builder_default_conv_tol_is_pinned():
    default = inspect.signature(build_molecular_hamiltonian).parameters["conv_tol"].default
    assert default == 1e-9, default


def table():
    """Print the measurements the spec's results section records (not a gate)."""
    print(f"{'geometry':15s} {'|dE0|':>8s} {'|dE_HF|':>8s} {'ov default':>10s} {'ov tight':>9s} "
          f"{'|d ov|':>8s} {'p0 default':>10s} {'p0 tight':>9s} {'|dE0reach|':>10s} {'idx':>5s} "
          f"{'nreach':>6s} {'leak default':>12s} {'leak tight':>10s}")
    for name, atom in {**CLEAN, **WITNESS}.items():
        d, t = _pair(atom)
        print(f"{name:15s} {abs(d['E0'] - t['E0']):8.1e} {abs(d['EHF'] - t['EHF']):8.1e} "
              f"{d['ov']:10.4g} {t['ov']:9.4g} {abs(d['ov'] - t['ov']):8.1e} "
              f"{d['p0']:10.2e} {t['p0']:9.2e} "
              f"{abs(d['E0reach'] - t['E0reach']):10.3e} {d['i_reach']:>2d}/{t['i_reach']:<2d} "
              f"{d['n_reach']:>3d}/{t['n_reach']:<2d} "
              f"{_leak(atom, None):12.2e} {_leak(atom, TIGHT_SCF_CONV_TOL):10.2e}", flush=True)


def residue_sweep(a=1.10):
    """Print the lowest level's HF population vs conv_tol at square-H4 side `a` (not a gate). This is
    the platform-dependent quantity behind SPEC_reachability_tolerance section 10."""
    print(f"a={a}: lowest-level HF population p0 vs conv_tol; window = (1e-10, 1e-8)")
    for ct in (1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13):
        p0 = _measure(_square(a), ct)["p0"]
        print(f"  conv_tol={ct:7.0e}  p0={p0:10.3e}  in window: {1e-10 < p0 < 1e-8}", flush=True)
