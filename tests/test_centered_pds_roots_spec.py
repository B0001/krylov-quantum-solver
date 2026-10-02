"""
Acceptance gates G1-G3 for specs/SPEC_centered_pds_roots.md (chem-dcz): does the centered moment
frame remove the raw-frame PDS excited-root bound violation on stretched H4?

The raw violation (SPEC_pds_excited_roots R3) looked "bimodal across fresh processes" even with BLAS
threads pinned. The hidden variable is PYTHONHASHSEED: hash-randomized set/dict iteration upstream
changes the Pauli-term order of the qubit Hamiltonian, which changes the last bits of the moments,
which at cond(M) ~ 1e17 selects the path. PYTHONHASHSEED must be set before the interpreter starts,
so each configuration runs in a child process (this file as __main__) with EVERYTHING pinned:
hash seed, BLAS/OpenMP threads, geometry, basis, K range. Given those, the result is deterministic.

PySCF/qiskit only (no block2); `make gates` runs it in its own process.
"""
import hashlib
import json
import os
import subprocess
import sys

import numpy as np

GEOMETRY = "H 0 0 0; H 0 0 2.0; H 0 0 4.0; H 0 0 6.0"  # SPEC_pds_excited_roots' stretched H4, sto3g
KS = range(3, 9)
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
VIOLATING_SEED, CLEAN_SEED = 2, 0  # measured: seeds {2,4,5,6,8,10,14} violate, the other 9 of 0..15 don't


def _child():
    """Run in the pinned child: per K, min_j(root_j - reachable_j) raw and centered, + cond(M_raw)."""
    from device_odmd import centered_frame
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
    from moment_expansion import hamiltonian_moments, pds_roots

    mh = build_molecular_hamiltonian(atom=GEOMETRY)
    w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    hf = np.asarray(mh.hf_state().data, dtype=complex)
    reach = np.sort(w[np.abs(V.conj().T @ hf) ** 2 > 1e-8].real + mh.energy_offset)
    raw, off = hamiltonian_moments(mh, 2 * max(KS))
    cen, off_c = hamiltonian_moments(centered_frame(mh)[0], 2 * max(KS))

    def margin(mom, o, K):
        r = pds_roots(mom, K, o)
        n = min(len(r), len(reach))
        return float(np.min(r[:n] - reach[:n]))

    def cond(mom, K):
        return float(np.linalg.cond(
            np.array([[mom[2 * K - i - j] for j in range(1, K + 1)] for i in range(1, K + 1)])))

    print(json.dumps({
        "raw_hash": hashlib.sha1(raw.tobytes()).hexdigest(),
        "raw": {K: margin(raw, off, K) for K in KS},
        "cen": {K: margin(cen, off_c, K) for K in KS},
        "cond_raw": {K: cond(raw, K) for K in KS},
        "cond_cen": {K: cond(cen, K) for K in KS},
    }))


_cache = {}


def _run(seed, threads=1):
    """Pinned configuration -> child's JSON result (keys K come back as strings)."""
    if (seed, threads) not in _cache:
        env = dict(os.environ, PYTHONHASHSEED=str(seed), **{v: str(threads) for v in THREAD_VARS})
        out = subprocess.run([sys.executable, __file__], env=env, check=True, capture_output=True,
                             text=True, cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        _cache[seed, threads] = json.loads(out.stdout.strip().splitlines()[-1])
    return _cache[seed, threads]


def test_G1_raw_violation_reproduces_in_a_pinned_configuration():
    """Pinned (PYTHONHASHSEED=2, 1 thread): raw PDS(7) puts a root below its reachable target by
    mHa (measured -3.548 mHa on macOS/Accelerate; SPEC_pds_excited_roots measured -7.49 mHa on its
    container's BLAS -- the magnitude is LAPACK-dependent, the violation is not), at cond(M) > 1e14.
    Re-running the same configuration gives bit-identical moments and margins."""
    a = _run(VIOLATING_SEED)
    assert a["raw"]["7"] < -1e-3, a["raw"]
    assert a["cond_raw"]["7"] > 1e14, a["cond_raw"]
    assert _cache.pop((VIOLATING_SEED, 1)) == _run(VIOLATING_SEED) == a   # fresh process, same bits


def test_G1b_the_hidden_variable_is_the_hash_seed_not_threads():
    """The 'bimodality' is PYTHONHASHSEED: seed 0 builds different moment bits and its raw PDS(7)
    respects the bound (+0.026 mHa); seed 2 at 4 threads reproduces seed 2 at 1 thread exactly."""
    v, c = _run(VIOLATING_SEED), _run(CLEAN_SEED)
    assert v["raw_hash"] != c["raw_hash"]
    assert c["raw"]["7"] >= -1e-9, c["raw"]
    assert _run(VIOLATING_SEED, threads=4) == v


def test_G2_centered_roots_stay_above_reachable_targets_to_K8_on_both_paths():
    """THE CLAIM: on both raw paths (violating and clean seed), every centered root_j >= reachable_j
    - 1e-9 Ha at K=3..8 (measured min +0.0042 mHa at K=8), and the centered margins agree across the
    two paths to < 1e-8 Ha -- centering makes the result insensitive to the bits that split raw."""
    v, c = _run(VIOLATING_SEED), _run(CLEAN_SEED)
    for r in (v, c):
        assert all(m >= -1e-9 for m in r["cen"].values()), r["cen"]
    assert all(abs(v["cen"][k] - c["cen"][k]) < 1e-8 for k in v["cen"]), (v["cen"], c["cen"])


def test_G3_centered_frame_is_well_inside_float64_where_raw_is_not():
    """Why G2 holds: at K=7 raw cond(M) > 1e14 on both paths, centered cond(M) < 1e10."""
    for r in (_run(VIOLATING_SEED), _run(CLEAN_SEED)):
        assert r["cond_raw"]["7"] > 1e14 and r["cond_cen"]["7"] < 1e10, (r["cond_raw"], r["cond_cen"])


if __name__ == "__main__":
    sys.path.insert(0, os.getcwd())
    _child()
