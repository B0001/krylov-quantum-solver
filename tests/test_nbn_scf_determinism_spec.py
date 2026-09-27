"""
Acceptance gates for specs/SPEC_nbn_scf_determinism.md.

Claim: benchmark_nbn.ground_state_mf (via _tight_scf's Newton/SOSCF stability-following loop) is
now deterministic across ambient pyscf.lib.num_threads() settings. Before the chem-y6w fix, this
was false: chem-owr (sandbox-handoffs/chem-owr.md) measured 3/4 default-thread trials on the real
NbN geometry landing on one locally-stable UHF minimum and 1/4 on a different one, 0.548 mHa apart
in SCF energy -- a coin flip from floating-point non-associativity in multi-threaded BLAS
reductions inside pyscf's Davidson-based stability() call.

pyscf + qiskit only (benchmark_nbn imports hybrid_quantum_solver.molecular_hamiltonian, which
imports qiskit/qiskit_nature) -- NOT block2. Must not share a process with any DMRG gate (CLAUDE.md
OpenMP isolation note); own file, per the test_*_spec.py-runs-in-its-own-process convention.
"""
import os

import pytest
from pyscf import lib

import benchmark_nbn as bn

CIF = "specs/nbn_mp-2634.cif"
needs_cif = pytest.mark.skipif(not os.path.exists(CIF), reason=f"{CIF} missing")


def _fresh_e_tot(tmp_path, tag, ambient_threads):
    """One from-scratch SCF (distinct chkfile, so neither run restores another's cache)."""
    lib.num_threads(ambient_threads)
    bn.CHKFILE = str(tmp_path / f"determinism_{tag}.chk")
    mf = bn.ground_state_mf(CIF)
    return float(mf.e_tot)


@needs_cif
def test_G1_deterministic_across_ambient_thread_counts(tmp_path):
    """DEFINITION OF DONE: fresh runs at ambient thread counts 8 and 2 give bit-identical e_tot.
    Before the fix this varied run to run at a fixed ambient count (see G2) -- varying the ambient
    setting explicitly here makes the check itself reproducible rather than relying on getting
    unlucky in CI."""
    e_8 = _fresh_e_tot(tmp_path, "ambient8", 8)
    e_2 = _fresh_e_tot(tmp_path, "ambient2", 2)
    assert e_8 == e_2, (e_8, e_2)


@needs_cif
def test_G2_deterministic_across_repeats(tmp_path):
    """Regression floor: three fresh runs at the SAME (default) ambient thread count are pairwise
    bit-identical. In miniature, this is chem-owr's own 4-trial experiment (which originally found
    3/4 vs 1/4) -- now made to pass every time, not most of the time."""
    lib.num_threads(lib.num_threads())  # no-op; keeps intent explicit (ambient = whatever it was)
    energies = [_fresh_e_tot(tmp_path, f"repeat{i}", 8) for i in range(3)]
    assert len(set(energies)) == 1, energies
