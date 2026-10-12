"""Gates for specs/SPEC_symmetric_scf.md (chem-n3w): fix the SCF-convergence residue at its
source, via an opt-in ``symmetric_scf`` flag on ``build_molecular_hamiltonian``, rather than
vetoing it downstream (the fix-forward ``SPEC_eigenstate_reachability`` §0 deferred as a
follow-up).

G6/G7 touch SCF trajectories directly and must be thread-pinned the way
``test_hamiltonian_determinism.py`` (chem-8tl) pins them: multithreaded BLAS reductions are a
separate, already-documented source of run-to-run bit noise in this exact SCF. Set
OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1 in the environment
this file runs in (``make gates`` / CI already does not parallelize within a pytest process; for a
manual single-file run use:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
        uv run pytest tests/test_symmetric_scf_spec.py -v
"""
import inspect
import sys
import warnings
from functools import lru_cache

import numpy as np
import pytest

sys.path.insert(0, "tests")
import test_eigenstate_reachability_spec as ter  # noqa: E402  (path insert above)

from hf_overlap_certificate import exact_reachable_overlap  # noqa: E402
from hybrid_quantum_solver.molecular_hamiltonian import (  # noqa: E402
    _point_group,
    build_molecular_hamiltonian,
)
from reachability import (  # noqa: E402
    REACHABLE_TOL_CERTIFIED,
    TIGHT_SCF_CONV_TOL,
    _dense_hf_projection,
    reachable_eigenpairs,
)

ASYM_H4 = "H 0 0 0; H 1.1 0.2 0.3; H 1.9 1.3 0.1; H 0.3 1.7 0.9"  # true C1: no mirror plane
WITNESS_A = 1.2  # square H4 side; Linux-native residue witness (SPEC_symmetric_scf R2)


def _square(a):
    return f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0"


def _chain(d):
    return f"H 0 0 0; H 0 0 {d}; H 0 0 {2 * d}; H 0 0 {3 * d}"


# Reused verbatim from SPEC_scf_conv_tol.CLEAN (test modules are not importable as a package here).
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
SIDES = (1.05, 1.10, 1.19, 1.2, 1.35, 1.4)


@lru_cache(maxsize=None)
def _build(atom, conv_tol, symmetric_scf):
    return build_molecular_hamiltonian(atom=atom, conv_tol=conv_tol, symmetric_scf=symmetric_scf)


# --- G1: defaults pinned --------------------------------------------------------------------------

def test_G1_defaults_are_pinned():
    params = inspect.signature(build_molecular_hamiltonian).parameters
    assert params["conv_tol"].default == 1e-9
    assert params["symmetric_scf"].default is False


# --- G2: point-group detection matches PySCF, no SCF ----------------------------------------------

@pytest.mark.parametrize("atom, expect", [
    ("H 0 0 0; H 0 0 0.74", "Dooh"),
    (_chain(1.0), "Dooh"),
    (_square(1.2), "D2h"),
    (ASYM_H4, "C1"),
])
def test_G2_point_group_matches_pyscf(atom, expect):
    assert _point_group(atom, "sto3g", 0, 0) == expect


# --- G3: no-op on an asymmetric (C1) geometry ------------------------------------------------------

def test_G3_noop_on_c1_geometry():
    off = _build(ASYM_H4, 1e-9, False)
    on = _build(ASYM_H4, 1e-9, True)
    assert on.build_args["conv_tol"] == off.build_args["conv_tol"] == 1e-9
    assert on.ground_state_energy() == off.ground_state_energy()
    assert on.hf_energy == off.hf_energy


# --- G4: harmless on the already-validated clean set -----------------------------------------------

@pytest.mark.parametrize("name", list(CLEAN))
def test_G4_harmless_on_clean_geometry(name):
    atom = CLEAN[name]
    off = _build(atom, 1e-9, False)
    on = _build(atom, 1e-9, True)
    assert abs(on.ground_state_energy() - off.ground_state_energy()) < 1e-9, name
    assert abs(on.hf_energy - off.hf_energy) < 1e-8, name


# --- G5: equivalent to the manual TIGHT_SCF_CONV_TOL workaround ------------------------------------

@pytest.mark.parametrize("a", SIDES)
def test_G5_equivalent_to_manual_tight_tolerance(a):
    atom = _square(a)
    opt = _build(atom, 1e-6, True)                      # loose input, flag floors it
    tight = _build(atom, TIGHT_SCF_CONV_TOL, False)      # the existing manual workaround
    assert opt.build_args["conv_tol"] == tight.build_args["conv_tol"] == TIGHT_SCF_CONV_TOL
    assert opt.ground_state_energy() == tight.ground_state_energy()
    assert opt.hf_energy == tight.hf_energy
    assert exact_reachable_overlap(opt) == exact_reachable_overlap(tight)


# --- G6: residue collapse at the Linux-native witness, thread-pinned (definition of done, 1/2) ----

def test_G6_residue_collapses_at_the_witness():
    atom = _square(WITNESS_A)
    for _ in range(5):
        off = build_molecular_hamiltonian(atom=atom, conv_tol=1e-9)
        _, _, pops_off = _dense_hf_projection(off)
        assert pops_off[0] > REACHABLE_TOL_CERTIFIED, pops_off[0]
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter("always")
            reachable_eigenpairs(off, tol=REACHABLE_TOL_CERTIFIED)
        assert any(issubclass(r.category, RuntimeWarning) for r in rec)

        on = build_molecular_hamiltonian(atom=atom, conv_tol=1e-9, symmetric_scf=True)
        _, _, pops_on = _dense_hf_projection(on)
        assert pops_on[0] < 1e-20, pops_on[0]
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter("always")
            reachable_eigenpairs(on, tol=REACHABLE_TOL_CERTIFIED)
        assert not any(issubclass(r.category, RuntimeWarning) for r in rec)


# --- G7: QKSD lands on the physical level, not the artifact (definition of done, 2/2) -------------

def test_G7_qksd_lands_on_the_physical_level_not_the_artifact():
    atom = _square(WITNESS_A)
    e_ag = ter._fci(atom, "Ag")[0]
    e_b1g = ter._fci(atom, "B1g")[0]
    assert e_b1g < e_ag - 0.1, (e_ag, e_b1g)   # sanity: B1g is the lower, HF-unreachable level here

    off = build_molecular_hamiltonian(atom=atom, conv_tol=1e-9)
    e_off_tail = ter._qksd_electronic(off)[[m - 1 for m in ter.G9_M]]   # same window as G9
    assert (e_off_tail - e_b1g).min() >= -1e-9      # variational floor around the artifact level
    assert (e_off_tail < e_ag - 0.1).all()          # ...and settled nowhere near the physical one
    assert reachable_eigenpairs(off)[0][0] == pytest.approx(e_b1g, abs=1e-8)

    on = build_molecular_hamiltonian(atom=atom, conv_tol=1e-9, symmetric_scf=True)
    e_on_tail = ter._qksd_electronic(on)[[m - 1 for m in ter.G9_M]]
    assert np.abs(e_on_tail - e_ag).max() < 1e-8, e_on_tail - e_ag
    assert reachable_eigenpairs(on)[0][0] == pytest.approx(e_ag, abs=1e-8)


def _child_pinned_check():
    """Not a gate: run manually to confirm the environment is actually thread-pinned for G6/G7."""
    import os

    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS"):
        print(v, os.environ.get(v))


if __name__ == "__main__":
    _child_pinned_check()
