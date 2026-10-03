"""
build_molecular_hamiltonian is bit-for-bit independent of PYTHONHASHSEED (chem-8tl).

The fermionic operator's key order follows Python's randomized hashing, and the Jordan-Wigner
mapper sums colliding Pauli terms in that order -- so before map_canonical, stretched H4 built the
same 185 labels with a seed-dependent order and seed-dependent low coefficient bits, which at
cond(M) ~ 1e17 flipped raw PDS(7) across the variational bound for 7 of 16 seeds. The hash seed is
fixed at interpreter start, so each seed builds in its own child process (this file as __main__).
BLAS/OpenMP threads are pinned to 1: multi-threaded OpenBLAS reductions in the SCF are a SEPARATE
source of run-to-run bit noise (measured on Linux: same seed, 4 threads, 3 different digests).
"""
import hashlib
import os
import subprocess
import sys

GEOMETRY = "H 0 0 0; H 0 0 2.0; H 0 0 4.0; H 0 0 6.0"  # stretched H4, sto3g (SPEC_pds_excited_roots)
SEEDS = range(16)
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")


def _child():
    import numpy as np

    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    q = build_molecular_hamiltonian(atom=GEOMETRY).qubit_hamiltonian
    digest = hashlib.sha1(repr([p.to_label() for p in q.paulis]).encode())
    digest.update(np.ascontiguousarray(q.coeffs).tobytes())
    print(len(q), digest.hexdigest())


def _build(seed):
    env = dict(os.environ, PYTHONHASHSEED=str(seed), **{v: "1" for v in THREAD_VARS})
    out = subprocess.run([sys.executable, __file__], env=env, check=True, capture_output=True,
                         text=True, cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    return out.stdout.strip().splitlines()[-1]


def test_h4_qubit_hamiltonian_is_identical_across_hash_seeds():
    """Same term count, labels, order and coefficient bits for PYTHONHASHSEED 0..15."""
    results = {seed: _build(seed) for seed in SEEDS}
    assert len(set(results.values())) == 1, results
    assert results[0].split()[0] == "185", results[0]


if __name__ == "__main__":
    _child()
