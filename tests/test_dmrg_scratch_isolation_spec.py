"""
Acceptance gate for chem-c3p: concurrent DMRG gates must not share ``./.dmrg_tmp`` scratch.

``scripts/run_gates.sh`` runs spec gate files in separate PROCESSES (required for block2/pyscf
isolation -- CLAUDE.md), but in parallel. Before this fix, ``dmrg_energy`` and
``dmrg_energy_extrapolated`` both defaulted ``scratch="./.dmrg_tmp"`` -- the literal same path for
every caller that didn't pass one explicitly. Two such gates running concurrently (e.g.
``GATE_JOBS=2`` over ``test_hchain_largen2_spec.py`` + ``test_hchain_tdl_spec.py``, both default-
scratch callers) wrote block2 scratch files into the same directory and could hang
(``futex_do_wait``, near-zero CPU, 30+ min -- chem-mjz's sandbox session; see the bead).

G1-G2 are pure (no block2): they pin the ``_scratch_dir`` contract directly.
G3 needs block2 and reproduces the real failure mode: two DMRG runs launched as separate OS
processes (real concurrency, not threads) at the same time, both taking the default scratch.
"""
import os
import shutil
import subprocess
import sys

import pytest

from hybrid_quantum_solver.dmrg_reference import (
    DEFAULT_SCRATCH_ROOT,
    _scratch_dir,
    dmrg_available,
)

requires_dmrg = pytest.mark.skipif(not dmrg_available(), reason="block2 not installed")


def test_G1_default_scratch_is_unique_per_call_and_cleaned_up():
    """Two calls with scratch=None (the default) must get DIFFERENT directories, and each must be
    removed again once its ``with`` block exits -- otherwise concurrent callers can still collide,
    or ``.dmrg_tmp`` grows without bound across a long gate run."""
    with _scratch_dir(None) as d1:
        assert os.path.isdir(d1)
        assert os.path.commonpath([os.path.abspath(d1), os.path.abspath(DEFAULT_SCRATCH_ROOT)]) == \
            os.path.abspath(DEFAULT_SCRATCH_ROOT)
        with _scratch_dir(None) as d2:
            assert os.path.isdir(d2)
            assert d1 != d2                     # distinct even for two calls in the SAME process
            assert os.path.isdir(d1)             # d1 not torn down by entering d2
    assert not os.path.exists(d1)
    assert not os.path.exists(d2)


def test_G2_explicit_scratch_is_passed_through_and_never_cleaned_up():
    """An explicit ``scratch=`` path is exactly what callers that already namespace their own
    scratch rely on (nbn_dmrg_reference.py, benchmark_hubbard_lieb_wu.py): passed through unchanged,
    and left on disk after the call -- this fix must not start deleting it out from under them."""
    explicit = os.path.join(DEFAULT_SCRATCH_ROOT, "explicit_probe")
    os.makedirs(explicit, exist_ok=True)
    marker = os.path.join(explicit, "marker.txt")
    with open(marker, "w") as fh:
        fh.write("kept")
    try:
        with _scratch_dir(explicit) as d:
            assert d == explicit
        assert os.path.exists(marker), "explicit scratch must survive the call"
    finally:
        shutil.rmtree(explicit, ignore_errors=True)


def test_G2b_explicit_nested_scratch_is_created_with_missing_parents(tmp_path):
    """chem-qqu: block2 only creates the leaf scratch directory, so an explicit nested path whose
    parent does not exist yet (``./.dmrg_tmp/g8_free`` in a fresh checkout) made the first MPS save
    fail. ``_scratch_dir`` must create the whole path before handing it to the driver."""
    nested = str(tmp_path / "missing_parent" / "leaf")
    with _scratch_dir(nested) as d:
        assert d == nested
        assert os.path.isdir(d)
    assert os.path.isdir(nested), "explicit scratch must survive the call"


_CHILD_TEMPLATE = """
import contextlib
import sys
sys.path.insert(0, {repo!r})
import hybrid_quantum_solver.dmrg_reference as dr

_orig_scratch_dir = dr._scratch_dir

@contextlib.contextmanager
def _spy(scratch):
    with _orig_scratch_dir(scratch) as d:
        print("SCRATCH=" + d, flush=True)
        yield d

dr._scratch_dir = _spy

import numpy as np
t, U = 1.0, 4.0
h1 = np.array([[0.0, -t], [-t, 0.0]])
eri = np.zeros((2, 2, 2, 2))
eri[0, 0, 0, 0] = U
eri[1, 1, 1, 1] = U

e = dr.dmrg_energy(h1, eri, (1, 1), 0.0, bond_dims=(4, 8), n_sweeps=6, n_threads=1)
print(f"ENERGY={{e:.12f}}", flush=True)
"""


@requires_dmrg
def test_G3_two_concurrent_dmrg_processes_use_distinct_scratch_and_leave_none_behind(tmp_path):
    """Reproduces the real failure mode: two independent OS processes (like two `run_gates.sh`
    workers), both calling ``dmrg_energy`` with the default scratch, running AT THE SAME TIME.

    Exact reference: the half-filled two-site Hubbard dimer, ``E0 = (U - sqrt(U**2+16t**2))/2``
    (t=1, U=4 -> E0 = (4 - sqrt(80))/2 = -2.472135955...).
    """
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script = tmp_path / "child.py"
    script.write_text(_CHILD_TEMPLATE.format(repo=repo_root))

    os.makedirs(DEFAULT_SCRATCH_ROOT, exist_ok=True)
    before = set(os.listdir(DEFAULT_SCRATCH_ROOT))

    procs = [
        subprocess.Popen([sys.executable, str(script)], cwd=repo_root,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for _ in range(2)
    ]
    # Block on each process's first line: both have created (and printed) their scratch dir by
    # the time this returns, and both are still alive and mid-sweep -- a real concurrency window.
    scratch_lines = [p.stdout.readline().strip() for p in procs]
    scratch_dirs = []
    for line in scratch_lines:
        assert line.startswith("SCRATCH="), line
        d = line.split("=", 1)[1]
        scratch_dirs.append(d)
        assert os.path.isdir(d), f"{d} should exist while its DMRG run is in flight"

    assert scratch_dirs[0] != scratch_dirs[1], "concurrent default-scratch callers collided"

    outs = [p.stdout.read() for p in procs]
    rcs = [p.wait(timeout=120) for p in procs]
    for rc, out, sline in zip(rcs, outs, scratch_lines):
        assert rc == 0, f"child failed (rc={rc}):\n{sline}\n{out}"

    energies = []
    for out in outs:
        energy_lines = [ln for ln in out.splitlines() if ln.startswith("ENERGY=")]
        assert energy_lines, out
        energies.append(float(energy_lines[0].split("=", 1)[1]))

    e_exact = 0.5 * (4.0 - (4.0 ** 2 + 16.0 * 1.0 ** 2) ** 0.5)
    for e in energies:
        assert abs(e - e_exact) < 1e-8, (e, e_exact)

    for d in scratch_dirs:
        assert not os.path.exists(d), f"stale scratch left behind: {d}"
    after = set(os.listdir(DEFAULT_SCRATCH_ROOT))
    assert after == before, f"stale entries under {DEFAULT_SCRATCH_ROOT}: {after - before}"
