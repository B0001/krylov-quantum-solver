"""
Acceptance gates for specs/SPEC_driver_orbitals_scf_determinism.md (chem-50w).

Claim: reachability._driver_orbitals's own drv.run_pyscf() call -- the one chem-s6c's fix
deliberately left unpinned (specs/SPEC_reachability_scf_determinism.md) -- shares
chem-7cx/chem-s6c's thread-nondeterminism mechanism, and pinning pyscf.lib.num_threads(1) around
it (restored in a finally) removes it, FOR A FIXED mh. Before the fix, holding ambient threads fixed
at 2 and calling the unpinned _driver_orbitals(mh) ten times in a row (square H4, a=1.35,
conv_tol=1e-9, mh held fixed) gave found={yes,yes,yes,yes,NO,NO,yes,yes,yes,yes} with
run_pyscf call counts {3,1,3,1,5,5,2,2,5,5} -- see spec section 1/9 for the captured table.

Every gate here builds its own `mh` under `pyscf.lib.num_threads(1)` BEFORE varying ambient threads
only around the call under test -- this deliberately isolates "does _driver_orbitals's own RHF call
add nondeterminism on top of a fixed comparison target" (this bead's question) from "is the
comparison target (build_molecular_hamiltonian's own free RHF call) itself reproducible" (chem-7cx's
separate, still-open question on this branch -- see spec section 2/7/R1). Building `mh` under
num_threads(1) is itself deterministic (confirmed directly: 8 repeated builds of the same witness
give a bit-identical qubit_hamiltonian and hf_energy), so this is a legitimate fixed test fixture,
not a tautology.

pyscf + qiskit-nature only (imports reachability and hybrid_quantum_solver.molecular_hamiltonian)
-- NOT block2. Must not share a process with any DMRG gate (CLAUDE.md OpenMP isolation note); own
file, per the test_*_spec.py-runs-in-its-own-process convention. G1/G2 use subprocesses regardless,
so this holds even if invoked alongside other files by accident.
"""
import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (a, conv_tol) witnesses -- the same square-H4 family test_eigenstate_reachability_spec.py and
# test_rebuild_gauge_search_spec.py already gate.
WITNESS_1_19 = (1.19, 1e-6)      # chem-0yh's own flaky witness
WITNESS_1_10 = (1.10, 1e-6)
WITNESS_1_35 = (1.35, 1e-9)      # spec section 1/9's captured repeat-call table


def _sq(a):
    return f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0"


def _build_fixed_mh(a, conv_tol):
    """mh built under num_threads(1) -- a fixed, reproducible comparison target, isolating this
    gate from build_molecular_hamiltonian's own separate (chem-7cx, unmerged-here) nondeterminism."""
    from pyscf import lib

    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    ambient = lib.num_threads()
    lib.num_threads(1)
    try:
        return build_molecular_hamiltonian(atom=_sq(a), conv_tol=conv_tol)
    finally:
        lib.num_threads(ambient)


def _child_single_call():
    """One fresh, counted _driver_orbitals(mh) call at this process's current ambient thread
    count, for the witness named on argv. Prints '<found> <calls>'."""
    from pyscf import lib
    from qiskit_nature.second_q.drivers import PySCFDriver

    import reachability

    a, conv_tol = float(sys.argv[2]), float(sys.argv[3])
    ambient_threads = int(sys.argv[4])

    mh = _build_fixed_mh(a, conv_tol)

    lib.num_threads(ambient_threads)
    real_run_pyscf = PySCFDriver.run_pyscf
    calls = []

    def counting_run_pyscf(self, *a_, **kw):
        calls.append(1)
        return real_run_pyscf(self, *a_, **kw)

    PySCFDriver.run_pyscf = counting_run_pyscf
    found = reachability._driver_orbitals(mh)
    mo_repr = repr(found[1].round(10).tolist()) if found is not None else "None"
    print(f"{found is not None} {len(calls)} {mo_repr}")


def _call_in_subprocess(a, conv_tol, ambient_threads):
    out = subprocess.run(
        [sys.executable, __file__, "--child", str(a), str(conv_tol), str(ambient_threads)],
        check=True, capture_output=True, text=True, cwd=REPO_ROOT)
    line = out.stdout.strip().splitlines()[-1]
    found_str, calls_str, mo_repr = line.split(" ", 2)
    return found_str == "True", int(calls_str), mo_repr


# --- G1: DEFINITION OF DONE -- cross-ambient-thread determinism, fixed mh -------------------------

def test_G1_deterministic_across_ambient_thread_counts():
    """For a fixed mh, _driver_orbitals's found/MO-coefficient outcome is the same across ambient
    thread counts in {1, 2, 4, 8} -- before the fix, a=1.19/conv_tol=1e-6 flipped found itself."""
    for a, conv_tol in (WITNESS_1_19, WITNESS_1_10, WITNESS_1_35):
        results = {n: _call_in_subprocess(a, conv_tol, n) for n in (1, 2, 4, 8)}
        outcomes = {n: (found, mo) for n, (found, _calls, mo) in results.items()}
        assert len(set(outcomes.values())) == 1, (a, conv_tol, results)


# --- G2: repeat-call determinism at one fixed ambient setting (regression floor) ------------------

def test_G2_deterministic_across_repeats_separate_processes():
    """Ten fresh calls, each its own process, ambient threads fixed at 2 (spec section 1/9's
    captured table's own setting), fixed mh, are pairwise bit-identical."""
    a, conv_tol = WITNESS_1_35
    results = [_call_in_subprocess(a, conv_tol, 2) for _ in range(10)]
    outcomes = {(found, mo) for found, _calls, mo in results}
    assert len(outcomes) == 1, results


# --- G3: same-process repeat-call determinism, thread count restored ------------------------------

def test_G3_deterministic_same_process_and_threads_restored():
    from pyscf import lib

    import reachability

    a, conv_tol = WITNESS_1_35
    mh = _build_fixed_mh(a, conv_tol)

    ambient_before = lib.num_threads()
    results = []
    for _ in range(10):
        found = reachability._driver_orbitals(mh)
        assert lib.num_threads() == ambient_before, lib.num_threads()
        results.append(None if found is None else found[1].round(10).tobytes())
    assert len(set(results)) == 1, results
    assert results[0] is not None, "witness should match on this host (see spec section 9)"


# --- G4: the retry loop stops firing when the only reason it fired was thread noise ---------------

def test_G4_exactly_one_run_pyscf_call_across_ambient_thread_counts():
    """The explicit before/after retry-count measurement chem-s6c's notes asked for: for a fixed
    mh, the retry loop makes exactly one run_pyscf call regardless of ambient thread count, for
    every witness this module's own gate file exercises. Before the fix this ranged 1-5 (spec
    section 1/9's captured table)."""
    for a, conv_tol in (WITNESS_1_19, WITNESS_1_10, WITNESS_1_35):
        calls_by_threads = {n: _call_in_subprocess(a, conv_tol, n)[1] for n in (1, 2, 4, 8)}
        assert set(calls_by_threads.values()) == {1}, (a, conv_tol, calls_by_threads)


# --- G5: no regression on SPEC_rebuild_gauge_search.md's own mechanics ----------------------------

def test_G5_rebuild_gauge_search_gates_still_pass():
    """tests/test_rebuild_gauge_search_spec.py's G3-G5 (cheap-when-unneeded call counting, the
    energy pre-filter's retry-not-false-accept behavior, the oversized-active-space refusal) are
    untouched by this fix -- run as a regression check in its own subprocess."""
    out = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "--no-header",
         "tests/test_rebuild_gauge_search_spec.py"],
        cwd=REPO_ROOT, capture_output=True, text=True)
    assert out.returncode == 0, out.stdout[-3000:] + out.stderr[-2000:]


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--child":
        _child_single_call()
