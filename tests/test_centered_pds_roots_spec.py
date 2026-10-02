"""
Acceptance gates G1-G3 for specs/SPEC_centered_pds_roots.md (chem-dcz, revised for chem-8tl): does
the centered moment frame remove the raw-frame PDS excited-root bound violation on stretched H4?

The raw violation is a function of the LAST BITS of the qubit Hamiltonian: the Jordan-Wigner mapper
sums colliding Pauli terms in the fermionic operator's key order, and at cond(M) ~ 1e17 those bits
pick the path. Before chem-8tl that order followed PYTHONHASHSEED; build_molecular_hamiltonian is
now canonical, so the paths are pinned instead by mapping the same fermionic terms in a seeded
SHUFFLED order (SHUFFLES; None = the canonical build). Each configuration runs in a child process
(this file as __main__) with hash seed and BLAS/OpenMP threads pinned; given those it is
deterministic.

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
SHUFFLES = range(16)
VIOLATING = "15"  # measured (Linux/OpenBLAS): shuffles 4 and 15 violate (-1.33, -5.91 mHa), 14 don't


def _child():
    """Per path (canonical + each shuffle): raw moment hash, per-K min_j(root_j - reachable_j) raw
    and centered, and cond(M) raw and centered."""
    import dataclasses

    from qiskit_nature.second_q.drivers import PySCFDriver
    from qiskit_nature.second_q.mappers import JordanWignerMapper
    from qiskit_nature.second_q.operators import FermionicOp

    from device_odmd import centered_frame
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
    from moment_expansion import hamiltonian_moments, pds_roots

    mh = build_molecular_hamiltonian(atom=GEOMETRY)
    w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    hf = np.asarray(mh.hf_state().data, dtype=complex)
    reach = np.sort(w[np.abs(V.conj().T @ hf) ** 2 > 1e-8].real + mh.energy_offset)
    fermionic = PySCFDriver(atom=GEOMETRY, basis="sto3g").run().hamiltonian.second_q_op()
    items = sorted(fermionic.items())

    def margin(mom, o, K):
        r = pds_roots(mom, K, o)
        n = min(len(r), len(reach))
        return float(np.min(r[:n] - reach[:n]))

    def cond(mom, K):
        return float(np.linalg.cond(
            np.array([[mom[2 * K - i - j] for j in range(1, K + 1)] for i in range(1, K + 1)])))

    def path(m):
        raw, off = hamiltonian_moments(m, 2 * max(KS))
        cen, off_c = hamiltonian_moments(centered_frame(m)[0], 2 * max(KS))
        return {
            "raw_hash": hashlib.sha1(raw.tobytes()).hexdigest(),
            "raw": {K: margin(raw, off, K) for K in KS},
            "cen": {K: margin(cen, off_c, K) for K in KS},
            "cond_raw": {K: cond(raw, K) for K in KS},
            "cond_cen": {K: cond(cen, K) for K in KS},
        }

    out = {"canonical": path(mh)}
    for s in SHUFFLES:
        perm = np.random.default_rng(s).permutation(len(items))
        op = FermionicOp(dict(items[i] for i in perm), num_spin_orbitals=fermionic.num_spin_orbitals)
        out[str(s)] = path(dataclasses.replace(mh, qubit_hamiltonian=JordanWignerMapper().map(op)))
    print(json.dumps(out))


_cache = {}


def _run(seed, threads=1):
    """Pinned configuration -> child's JSON result (keys come back as strings)."""
    if (seed, threads) not in _cache:
        env = dict(os.environ, PYTHONHASHSEED=str(seed), **{v: str(threads) for v in THREAD_VARS})
        out = subprocess.run([sys.executable, __file__], env=env, check=True, capture_output=True,
                             text=True, cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        _cache[seed, threads] = json.loads(out.stdout.strip().splitlines()[-1])
    return _cache[seed, threads]


def test_G1_raw_violation_reproduces_in_a_pinned_configuration():
    """Pinned (shuffle 15, hash seed 2, 1 thread): raw PDS(7) puts a root below its reachable target
    by mHa (measured -5.91 mHa on Linux/OpenBLAS; the magnitude is LAPACK-dependent, the violation
    is not), at cond(M) > 1e14. A fresh process gives bit-identical moments and margins."""
    a = _run(2)[VIOLATING]
    assert a["raw"]["7"] < -1e-3, a["raw"]
    assert a["cond_raw"]["7"] > 1e14, a["cond_raw"]
    _cache.pop((2, 1))
    assert _run(2)[VIOLATING] == a   # fresh process, same bits


def test_G1b_the_hidden_variable_is_summation_order_not_the_hash_seed():
    """The canonical build is the same bits under hash seeds 0 and 2 (chem-8tl), and so is every
    shuffle; the summation order alone splits raw PDS(7) into violating and bound-respecting paths."""
    a, b = _run(2), _run(0)
    assert a == b
    k7 = [a[s]["raw"]["7"] for s in a]
    assert min(k7) < -1e-3 and max(k7) >= -1e-9, k7
    assert len({a[s]["raw_hash"] for s in a}) > 1


def test_G2_centered_roots_stay_above_reachable_targets_to_K8_on_every_path():
    """THE CLAIM: on every path (canonical + 16 shuffles, violating ones included), every centered
    root_j >= reachable_j - 1e-9 Ha at K=3..8, and the centered margins agree across paths to
    < 1e-8 Ha -- centering makes the result insensitive to the bits that split raw."""
    r = _run(2)
    ref = r["canonical"]["cen"]
    for p in r.values():
        assert all(m >= -1e-9 for m in p["cen"].values()), p["cen"]
        assert all(abs(p["cen"][k] - ref[k]) < 1e-8 for k in ref), (p["cen"], ref)


def test_G3_centered_frame_is_well_inside_float64_where_raw_is_not():
    """Why G2 holds: at K=7 raw cond(M) > 1e14 on every path, centered cond(M) < 1e10."""
    for p in _run(2).values():
        assert p["cond_raw"]["7"] > 1e14 and p["cond_cen"]["7"] < 1e10, (p["cond_raw"], p["cond_cen"])


if __name__ == "__main__":
    sys.path.insert(0, os.getcwd())
    _child()
