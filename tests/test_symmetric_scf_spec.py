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
from pyscf import lib

sys.path.insert(0, "tests")
import test_eigenstate_reachability_spec as ter  # noqa: E402  (path insert above)

import hybrid_quantum_solver.molecular_hamiltonian as mham  # noqa: E402
from hf_overlap_certificate import exact_reachable_overlap  # noqa: E402
from hybrid_quantum_solver.molecular_hamiltonian import (  # noqa: E402
    _point_group,
    build_dipole_operators,
    build_molecular_hamiltonian,
)
from reachability import (  # noqa: E402
    REACHABLE_TOL_CERTIFIED,
    TIGHT_SCF_CONV_TOL,
    _dense_hf_projection,
    reachable_eigenpairs,
)

ASYM_H4 = "H 0 0 0; H 1.1 0.2 0.3; H 1.9 1.3 0.1; H 0.3 1.7 0.9"  # true C1: no mirror plane
# Square H4 side 1.2: the Linux-native residue witness (spec R2), whose free SCF stopped on the D2h
# saddle on Linux but not on macOS. Since chem-8ls the builder follows that saddle's instability to
# the minimum on both (specs/SPEC_scf_stability_rule.md), so G6/G7 run everywhere (spec section 10).
WITNESS_A = 1.2


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
    dipole = inspect.signature(build_dipole_operators).parameters  # same defaults -> same SCF
    assert dipole["conv_tol"].default == 1e-9
    assert dipole["symmetric_scf"].default is False


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


# --- G6: no lowest-level residue left for the flag to collapse (chem-8ls; was definition of done) --

def test_G6_no_residue_with_or_without_the_flag():
    """Was (Linux only): without the flag the lowest level carried a 5.08e-9 residue and was warned
    about, and with it the residue was < 1e-20, repeatably over 5 builds each. That residue was the
    free SCF stopped on the D2h saddle. Since chem-8ls both builds are the broken-symmetry minimum,
    whose lowest level is physical, so with or without the flag p0 is the same, > 0.4, and nothing
    is warned about. Still 5 builds each: the original's repeatability check."""
    atom = _square(WITNESS_A)
    for _ in range(5):
        p0 = []
        for flag in (False, True):
            mh = build_molecular_hamiltonian(atom=atom, conv_tol=1e-9, symmetric_scf=flag)
            p0.append(float(_dense_hf_projection(mh)[2][0]))
            with warnings.catch_warnings(record=True) as rec:
                warnings.simplefilter("always")
                reachable_eigenpairs(mh, tol=REACHABLE_TOL_CERTIFIED)
            assert not any(issubclass(r.category, RuntimeWarning) for r in rec), flag
        assert p0[0] == pytest.approx(p0[1], abs=1e-4) and min(p0) > 0.4, p0


# --- G7: QKSD lands on the B1g ground state either way (chem-8ls; was definition of done) ----------

def test_G7_qksd_lands_on_the_ground_state_with_or_without_the_flag():
    """Was (Linux only): QKSD from the default build sat on B1g (the residue level), and with the
    flag on Ag. With the stability rule HF genuinely overlaps the B1g ground state, so QKSD sits on
    it, and the reference starts there, with or without the flag."""
    atom = _square(WITNESS_A)
    e_ag = ter._fci(atom, "Ag")[0]
    e_b1g = ter._fci(atom, "B1g")[0]
    assert e_b1g < e_ag - 0.1, (e_ag, e_b1g)   # sanity: B1g is the lower level here

    for flag in (False, True):
        mh = build_molecular_hamiltonian(atom=atom, conv_tol=1e-9, symmetric_scf=flag)
        tail = ter._qksd_electronic(mh)[[m - 1 for m in ter.G9_M]]   # same window as G9
        assert (tail - e_b1g).min() >= -1e-9, (flag, tail - e_b1g)   # the variational floor
        assert (tail - e_b1g).max() < 1e-8, (flag, tail - e_b1g)     # ...and converged on B1g
        assert reachable_eigenpairs(mh)[0][0] == pytest.approx(e_b1g, abs=1e-8), flag


# --- G8: the dipole runs H's SCF -- same conv_tol floor, same thread pin (chem-9ud) ---------------

@pytest.mark.parametrize("a", SIDES)
def test_G8_dipole_shares_the_hamiltonians_scf_orbitals(a, monkeypatch):
    """H and dipole built with the same symmetric_scf=True arguments, at ambient PySCF thread
    counts 2 and 8: every SCF runs at the floored conv_tol, on one thread, to bit-identical MO
    coefficients. Before chem-9ud the dipole ignored the floor (raw 1e-9) and the pin; on macOS its
    orbitals then differed from H's at 5 of these 6 sides (whole MO columns sign-flipped at
    1.05/1.2/1.35, putting mu_x/mu_y 4.4-5.4 a.u. off in spectral norm). The thread leg can only
    fail where PySCF has OpenMP (Linux): the macOS wheel has none, so lib.num_threads(n) is a
    no-op there (spec R4)."""
    runs = []  # (conv_tol, threads inside the SCF, mo_coeff), one per PySCFDriver run

    class Recorder(mham.PySCFDriver):
        # _driver_problem runs run_pyscf, the stability rule (chem-8ls), then to_problem, all inside
        # the pin. to_problem therefore sees the final orbitals that H and the dipole are built on.
        def to_problem(self, *args, **kwargs):
            problem = super().to_problem(*args, **kwargs)
            runs.append((self.conv_tol, lib.num_threads(), self._calc.mo_coeff.copy()))
            return problem

    monkeypatch.setattr(mham, "PySCFDriver", Recorder)
    ambient = lib.num_threads()
    try:
        for n in (2, 8):
            lib.num_threads(n)
            build_molecular_hamiltonian(atom=_square(a), symmetric_scf=True)
            build_dipole_operators(atom=_square(a), symmetric_scf=True)
    finally:
        lib.num_threads(ambient)
    assert len(runs) == 4                       # H, dipole, H, dipole
    ref = runs[0][2]
    assert all(mo.tobytes() == ref.tobytes() for *_, mo in runs), \
        [float(np.abs(mo - ref).max()) for *_, mo in runs]
    assert [ct for ct, _, _ in runs] == [TIGHT_SCF_CONV_TOL] * 4
    assert [t for _, t, _ in runs] == [1] * 4


def _child_pinned_check():
    """Not a gate: run manually to confirm the environment is actually thread-pinned for G6/G7."""
    import os

    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS"):
        print(v, os.environ.get(v))


if __name__ == "__main__":
    _child_pinned_check()
