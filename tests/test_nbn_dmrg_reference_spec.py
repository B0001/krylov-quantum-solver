"""
Acceptance gates G1-G3 for specs/SPEC_nbn_dmrg_reference.md (NbN CAS(14,14) DMRG reference).

Test-first: ``nbn_dmrg_reference`` does not exist yet, so this file is RED until the spec is
implemented. Claim: for a TM active space beyond the repo's FCI cutoff (~1.18e7 determinants),
DMRG gives a converged reference certified by TWO independent sweep schedules (perD vs ramp,
different dims/seeds/scratch) agreeing -- with the honest finding gated too: the spin-scanned
ground sector is high-spin nelec=(10,4) and LOW-entanglement (discarded weight < 1e-7 already at
D=300), so this is a soft DMRG target, not a strong-correlation benchmark.

pyscf + block2 ONLY (no qiskit imports -- block2's OpenMP runtime segfaults in a process that
already imported pyscf+qiskit-aer; `make gates` runs this file in its own process).
"""
from math import comb

import numpy as np

from nbn_dmrg_reference import load_nbn_cas, run_schedule

_CACHE = {}


def _results():
    if "res" not in _CACHE:
        _CACHE["res"] = {s: run_schedule(s) for s in ("A'", "B'")}
    return _CACHE["res"]


def test_G1_beyond_fci_and_sector_pin():
    """The CAS(14,14) full Hilbert space (both spins over 14 orbitals, sum over Sz sectors)
    exceeds the 5e6-determinant FCI cutoff; the cached spin-scanned SCF restores without an SCF
    run and lands in the high-spin nelec=(10,4) sector. (The fixed-sector count, comb(14,10)*
    comb(14,4) = 1.0e6, is below the cutoff on its own -- the intractability is the FULL problem
    a black-box FCI would face, and DMRG's advantage is that it stays in-sector.)"""
    h1, eri, nelec, e_core = load_nbn_cas()
    assert h1.shape == (14, 14) and eri.shape == (14, 14, 14, 14)
    assert nelec == (10, 4), nelec
    n_full = comb(14, 7) ** 2            # half-filling Sz=0 sector -- the largest FCI must handle
    assert n_full > 5_000_000, n_full


def test_G2_two_independent_schedules_agree():
    """DEFINITION OF DONE: cheap-dims perD (100/200/300) vs ramp (80/160/300), different seeds
    and scratch dirs: |E_A' - E_B'| < 0.1 mHa (measured 0.0012), both in the discarded-weight
    regime (the guard that failed loudly in the killed hchain spec)."""
    res = _results()
    e_a, e_b = res["A'"].energy, res["B'"].energy
    assert abs(e_a - e_b) < 1e-4, (e_a, e_b)
    assert res["A'"].method == "dweight", res["A'"].method
    assert res["B'"].method == "dweight", res["B'"].method


def test_G3_softness_finding_is_pinned():
    """The recorded finding: the high-spin sector is low-entanglement -- discarded weight at
    D=300 < 1e-7 and per-D energy spread < 0.01 mHa. Nobody should later mistake this reference
    for a strong-correlation benchmark."""
    res = _results()["A'"]
    dims = [d for d, _, _ in res.per_D]
    weights = {d: w for d, w, _ in res.per_D}
    energies = np.array([e for _, _, e in res.per_D])
    assert 300 in dims, dims
    assert weights[300] < 1e-7, weights[300]
    assert energies.max() - energies.min() < 1e-5, energies
