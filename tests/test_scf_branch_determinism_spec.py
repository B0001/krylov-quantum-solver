"""
Acceptance gates for specs/SPEC_scf_branch_determinism.md.

Claim: build_molecular_hamiltonian (via pyscf.lib.num_threads(1) pinned around its PySCFDriver
call) is now deterministic across ambient pyscf.lib.num_threads()/OMP_NUM_THREADS settings. Before
the chem-7cx fix, this was false on square H4 (sto-3g, near-degenerate): 9/10 default-thread trials
at a=1.1 landed on hf_energy=-1.7825510081..., 1/10 on a different locally-stable RHF solution,
-1.7825511822... -- floating-point non-associativity in PySCF's multi-threaded integral evaluation
flipping which branch the free (symmetry=False) SCF converges to. Same mechanism and fix pattern as
SPEC_nbn_scf_determinism.md (chem-y6w); see that spec and sandbox-handoffs/chem-owr.md for the
original discovery of this class of bug.

pyscf + qiskit only (imports hybrid_quantum_solver.molecular_hamiltonian, which imports
qiskit/qiskit_nature) -- NOT block2. Must not share a process with any DMRG gate (CLAUDE.md OpenMP
isolation note); own file, per the test_*_spec.py-runs-in-its-own-process convention. G1/G2/G4 use
subprocesses regardless, so this holds even if invoked alongside other files by accident.
"""
import hashlib
import os
import subprocess
import sys

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SQUARE_H4 = "H 0 0 0; H 1.1 0 0; H 1.1 1.1 0; H 0 1.1 0"  # chem-7cx's own witness, a = 1.1
CONV_TOL = 1e-6


def _digest(mh):
    """Bit-level fingerprint: hf_energy plus the qubit Hamiltonian's labels and coefficient bytes."""
    q = mh.qubit_hamiltonian
    h = hashlib.sha1(repr([p.to_label() for p in q.paulis]).encode())
    h.update(np.ascontiguousarray(q.coeffs).tobytes())
    return f"{mh.hf_energy!r} {h.hexdigest()}"


def _child():
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    mh = build_molecular_hamiltonian(atom=SQUARE_H4, conv_tol=CONV_TOL)
    print(_digest(mh))


def _build_in_subprocess(omp_threads=None):
    """One fresh build_molecular_hamiltonian call in its own interpreter."""
    env = dict(os.environ)
    if omp_threads is not None:
        env["OMP_NUM_THREADS"] = str(omp_threads)
    out = subprocess.run([sys.executable, __file__], env=env, check=True, capture_output=True,
                         text=True, cwd=REPO_ROOT)
    return out.stdout.strip().splitlines()[-1]


def test_G1_deterministic_across_ambient_thread_counts():
    """DEFINITION OF DONE: fresh builds at ambient OMP_NUM_THREADS 2 and 8 give a bit-identical
    digest (hf_energy + qubit Hamiltonian labels/coefficients). Before the fix this varied run to
    run at a FIXED ambient count (see G2) -- varying the ambient setting explicitly here makes the
    check itself reproducible rather than relying on getting unlucky in CI."""
    d_2 = _build_in_subprocess(omp_threads=2)
    d_8 = _build_in_subprocess(omp_threads=8)
    assert d_2 == d_8, (d_2, d_8)


def test_G2_deterministic_across_repeats_separate_processes():
    """Regression floor: ten fresh builds, each its own process at the default ambient thread
    count, are pairwise bit-identical. Directly reproduces, in full, chem-7cx's own 10-trial
    experiment (which originally found 9/10 vs 1/10) -- now made to pass every time."""
    digests = [_build_in_subprocess() for _ in range(10)]
    assert len(set(digests)) == 1, digests


def test_G3_deterministic_across_repeats_same_process_and_threads_restored():
    """The _system()-cache-masked case: ten fresh, uncached builds inside ONE process (what a
    single pytest process doing ten builds looks like, vs G2's ten separate processes) are
    pairwise bit-identical, and the finally-restore leaves pyscf's ambient thread count exactly as
    this process found it after every call (a caller elsewhere in the same process must not find
    its own thread count silently pinned to 1)."""
    from pyscf import lib

    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    ambient_before = lib.num_threads()
    digests = []
    for _ in range(10):
        mh = build_molecular_hamiltonian(atom=SQUARE_H4, conv_tol=CONV_TOL)
        digests.append(_digest(mh))
        assert lib.num_threads() == ambient_before, lib.num_threads()
    assert len(set(digests)) == 1, digests


_FLIPPING_WITNESSES = (
    "tests/test_eigenstate_reachability_spec.py::"
    "test_G1_residue_level_identified_warned_and_kept[1.1-1e-06-1e-08]",
    "tests/test_eigenstate_reachability_spec.py::test_G9_control_tight_scf_stays_on_the_ag_level[1.1]",
)


def _outcomes():
    """Run both chem-7cx witness tests in one pytest subprocess; return {nodeid: passed_bool}."""
    out = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "--no-header", "-rA", *_FLIPPING_WITNESSES],
        cwd=REPO_ROOT, capture_output=True, text=True)
    lines = out.stdout.splitlines()
    outcomes = {}
    for w in _FLIPPING_WITNESSES:
        passed = any(line.startswith("PASSED ") and w in line for line in lines)
        failed = any(line.startswith("FAILED ") and w in line for line in lines)
        assert passed ^ failed, (w, out.stdout[-2000:])
        outcomes[w] = passed
    return outcomes


def test_G4_known_flipping_witnesses_stabilize_pass_or_fail():
    """The two tests chem-7cx found flipping pass/fail across repeated full-suite runs now give
    the SAME outcome every time -- five separate pytest subprocess invocations. This does NOT
    assert which outcome: chem-g14 (a separate, pre-existing cross-platform/version divergence)
    governs whether these should be green on this host. It only asserts they stop being a coin
    flip."""
    runs = [_outcomes() for _ in range(5)]
    for w in _FLIPPING_WITNESSES:
        results = [r[w] for r in runs]
        assert len(set(results)) == 1, (w, results)


if __name__ == "__main__":
    _child()
