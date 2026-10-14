"""
Acceptance gates for specs/SPEC_reachability_scf_determinism.md.

Claim: reachability.scf_symmetry_status and reachability.hf_orbital_irreps (via
pyscf.lib.num_threads(1) pinned around their RHF .kernel() calls) are now deterministic across
ambient pyscf.lib.num_threads()/OMP_NUM_THREADS settings. Before the chem-s6c fix, this was false
on square H4 (sto-3g, a=1.10, conv_tol=1e-6): OMP_NUM_THREADS in {1, 4, 8} gave is_broken=True,
OMP_NUM_THREADS=2 gave is_broken=False -- the CLASSIFICATION itself flipped, not just low-digit
noise in dE. Same mechanism and fix pattern as SPEC_scf_branch_determinism.md (chem-7cx), applied
to reachability.py's own raw scf.RHF calls instead of build_molecular_hamiltonian's.

pyscf only (imports reachability, a top-level script, not hybrid_quantum_solver) -- NOT block2.
Must not share a process with any DMRG gate (CLAUDE.md OpenMP isolation note); own file, per the
test_*_spec.py-runs-in-its-own-process convention. G1/G2/G4 use subprocesses regardless, so this
holds even if invoked alongside other files by accident.
"""
import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SQUARE_H4_110 = "H 0 0 0; H 1.10 0 0; H 1.10 1.10 0; H 0 1.10 0"  # chem-s6c's own witness
WITNESS_CONV_TOL = 1e-6


def _child():
    import reachability

    print(repr(reachability.scf_symmetry_status(SQUARE_H4_110, conv_tol=WITNESS_CONV_TOL)))


def _call_in_subprocess(omp_threads=None):
    """One fresh scf_symmetry_status call on the witness, in its own interpreter."""
    env = dict(os.environ)
    if omp_threads is not None:
        env["OMP_NUM_THREADS"] = str(omp_threads)
    else:
        env.pop("OMP_NUM_THREADS", None)
    out = subprocess.run([sys.executable, __file__], env=env, check=True, capture_output=True,
                         text=True, cwd=REPO_ROOT)
    return out.stdout.strip().splitlines()[-1]


# --- G1: DEFINITION OF DONE -- cross-ambient-thread determinism -----------------------------------

def test_G1_deterministic_across_ambient_thread_counts():
    """Before the fix this witness flipped is_broken itself (not just dE's low digits) between
    OMP_NUM_THREADS in {1, 4, 8} (True) and 2 (False) -- see spec section 1's captured table."""
    results = {n: _call_in_subprocess(omp_threads=n) for n in (1, 2, 4, 8)}
    assert len(set(results.values())) == 1, results


# --- G2: repeat-call determinism at one ambient setting (regression floor) ------------------------

def test_G2_deterministic_across_repeats_separate_processes():
    """Ten fresh calls, each its own process at the default ambient thread count, are pairwise
    bit-identical."""
    results = [_call_in_subprocess() for _ in range(10)]
    assert len(set(results)) == 1, results


# --- G3: same-process repeat-call determinism, both functions, thread count restored --------------

def test_G3_scf_symmetry_status_deterministic_same_process_and_threads_restored():
    from pyscf import lib

    import reachability

    ambient_before = lib.num_threads()
    results = []
    for _ in range(10):
        results.append(reachability.scf_symmetry_status(SQUARE_H4_110, conv_tol=WITNESS_CONV_TOL))
        assert lib.num_threads() == ambient_before, lib.num_threads()
    assert len(set(results)) == 1, results


def test_G3_hf_orbital_irreps_deterministic_same_process_and_threads_restored():
    """LiH: not symmetry-broken, so this exercises the labelling branch's own extra RHF solve
    (mf_sym.kernel() inside hf_orbital_irreps), returning a real irrep list rather than None."""
    from pyscf import lib

    import reachability

    ambient_before = lib.num_threads()
    results = []
    for _ in range(10):
        irreps = reachability.hf_orbital_irreps("Li 0 0 0; H 0 0 1.6")
        assert irreps is not None
        results.append(tuple(irreps))
        assert lib.num_threads() == ambient_before, lib.num_threads()
    assert len(set(results)) == 1, results


# --- G4: the affected spec's gates stop flipping pass/fail -----------------------------------------

_AFFECTED_WITNESSES = (
    "tests/test_symmetry_reachability_spec.py::"
    "test_G1_filter_is_available_exactly_when_the_artifact_is_present",
    "tests/test_symmetry_reachability_spec.py::"
    "test_G3_unavailable_means_RHF_found_a_lower_solution[1.2]",
    "tests/test_symmetry_reachability_spec.py::"
    "test_G3_available_means_the_symmetric_solution_is_the_RHF_solution[1.05]",
    "tests/test_symmetry_reachability_spec.py::"
    "test_G3_available_means_the_symmetric_solution_is_the_RHF_solution[1.1]",
    "tests/test_symmetry_reachability_spec.py::"
    "test_G3_available_means_the_symmetric_solution_is_the_RHF_solution[1.35]",
)


def _outcomes():
    """Run the whole symmetry-reachability file (its own fixtures build the mh objects) in one
    pytest subprocess; return {nodeid: passed_bool} for the witnesses this spec cares about."""
    out = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "--no-header", "-rA",
         "tests/test_symmetry_reachability_spec.py"],
        cwd=REPO_ROOT, capture_output=True, text=True)
    lines = out.stdout.splitlines()
    outcomes = {}
    for w in _AFFECTED_WITNESSES:
        nodeid = w.split("::", 1)[1]
        passed = any(line.startswith("PASSED ") and nodeid in line for line in lines)
        failed = any(line.startswith("FAILED ") and nodeid in line for line in lines)
        assert passed ^ failed, (w, out.stdout[-2000:])
        outcomes[w] = passed
    return outcomes


def test_G4_affected_witnesses_stabilize_pass_or_fail():
    """The five test_symmetry_reachability_spec.py gates that call scf_symmetry_status/
    hf_orbital_irreps directly on the square-H4 family now give the SAME outcome every time --
    five separate pytest subprocess invocations of the whole file. Does NOT assert which outcome:
    chem-g14 (a separate, pre-existing cross-platform/version divergence) governs whether these
    should be green on this host (spec section 2). Only asserts they stop being a coin flip."""
    runs = [_outcomes() for _ in range(5)]
    for w in _AFFECTED_WITNESSES:
        results = [r[w] for r in runs]
        assert len(set(results)) == 1, (w, results)


if __name__ == "__main__":
    _child()
