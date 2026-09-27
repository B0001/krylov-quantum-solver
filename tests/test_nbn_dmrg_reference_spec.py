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
import os
from math import comb

import numpy as np
import pytest

from hybrid_quantum_solver.dmrg_reference import DISCARD_WEIGHT_FLOOR, dmrg_available
from nbn_dmrg_reference import ensure_reference_data, load_nbn_cas, run_schedule

# data/ is gitignored, so a fresh clone has no chkfile -- but chem-bbi vendored the real checkpoint
# and its source CIF at specs/nbn_scf_reference.chk / specs/nbn_mp-2634.cif; materialize them into
# the data/ paths load_nbn_cas expects. block2 is a genuinely optional extra; its absence is still
# a SKIP, not a FAIL (SPEC_extrap_regime.md R-note).
ensure_reference_data()
_CHK = "data/nbn_scf.chk"
needs_chk = pytest.mark.skipif(not os.path.exists(_CHK),
                                reason=f"{_CHK} missing and vendored copy absent too")
needs_dmrg = pytest.mark.skipif(not dmrg_available(), reason="block2 not installed")

_CACHE = {}


def _results():
    if "res" not in _CACHE:
        _CACHE["res"] = {s: run_schedule(s) for s in ("A'", "B'")}
    return _CACHE["res"]


@needs_chk
def test_G1_beyond_fci_and_sector_pin():
    """The CAS(14,14) full Hilbert space (both spins over 14 orbitals, sum over Sz sectors)
    exceeds the 5e6-determinant FCI cutoff; the cached spin-scanned SCF restores without an SCF
    run and lands in the high-spin nelec=(10,4) sector. (That is the UHF ground spin, not the CAS
    ground: S=1 is 1.538 mHa lower on the real, vendored geometry -- G3 below; see chem-bbi's
    correction in SPEC_nbn_dmrg_reference.md, which supersedes the 15.1 mHa figure recorded by
    chem-dc7 against a since-explained wrong-bond-length reconstruction.) (The fixed-sector
    count, comb(14,10)* comb(14,4) = 1.0e6, is below the cutoff on its own -- the intractability
    is the FULL problem a black-box FCI would face, and DMRG's advantage is that it stays
    in-sector.)"""
    h1, eri, nelec, e_core = load_nbn_cas()
    assert h1.shape == (14, 14) and eri.shape == (14, 14, 14, 14)
    assert nelec == (10, 4), nelec
    n_full = comb(14, 7) ** 2            # half-filling Sz=0 sector -- the largest FCI must handle
    assert n_full > 5_000_000, n_full


@needs_chk
@needs_dmrg
def test_G2_two_independent_schedules_agree():
    """DEFINITION OF DONE: cheap-dims perD (100/200/300) vs ramp (80/160/300), different seeds
    and scratch dirs: |E_A' - E_B'| < 0.1 mHa (measured 0.0012), both in the discarded-weight
    regime (the guard that failed loudly in the killed hchain spec)."""
    res = _results()
    e_a, e_b = res["A'"].energy, res["B'"].energy
    assert abs(e_a - e_b) < 1e-4, (e_a, e_b)
    # Was `method == "dweight"`. That assertion was the reason the CI dims had to be chosen
    # DOWNWARD: SPEC_nbn_dmrg_reference.md:55 records this spec's own headline run (D=400/800/1200,
    # weights 1e-9..1e-13) landing in "invD" precisely because it converged. The gate now excludes
    # only uncontrolled truncation, so raising D can no longer fail it.
    # No accuracy reference here (the agreement check above is A' vs B'), so corroborate a
    # "converged" verdict against the weights rather than trusting the label.
    for tag in ("A'", "B'"):
        assert res[tag].regime != "uncontrolled", (tag, res[tag].regime)
        dws = [w for _, w, _ in res[tag].per_D]
        assert res[tag].regime != "converged" or max(dws) <= DISCARD_WEIGHT_FLOOR, (tag, dws)


@needs_chk
@needs_dmrg
def test_G3_softness_finding_is_pinned():
    """The recorded finding: the high-spin (10,4) sector is itself low-entanglement --
    discarded weight at D=300 < 1e-7 and per-D energy spread < 0.01 mHa. That characterizes
    the *sector* G1/G2 pin, not the CAS ground state -- see G3b below."""
    res = _results()["A'"]
    dims = [d for d, _, _ in res.per_D]
    weights = {d: w for d, w, _ in res.per_D}
    energies = np.array([e for _, _, e in res.per_D])
    assert 300 in dims, dims
    assert weights[300] < 1e-7, weights[300]
    assert energies.max() - energies.min() < 1e-5, energies


@needs_chk
@needs_dmrg
def test_G3b_true_ground_is_pinned():
    """chem-bbi correction: the committed (10,4)/S=3 sector (G1/G2/G3 above) is not the CAS
    ground. Exact FCI on the real, vendored geometry (Ms=0 sector, unconstrained, all 1.18e7
    determinants, 1010.6s wall on 8 cores) finds S=1 (nelec=(8,6)) lower: E=-110.04756504307636
    Ha, <S^2>=2.0000095761529892 (results/nbn_low_spin/runs.jsonl, 2026-09-26T15:07:42). DMRG A'
    on the Sz-projected (8,6) sector must agree with that FCI number to sub-mHa and sit below the
    (10,4) sector's own energy."""
    e_fci_s1 = -110.04756504307636  # exact FCI, Ms=0 sector; results/nbn_low_spin/runs.jsonl
    e_s1 = run_schedule("A'", nelec=(8, 6), scratch_tag="_g3b_s1").energy
    e_s3 = _results()["A'"].energy
    assert abs(e_s1 - e_fci_s1) < 5e-4, (e_s1, e_fci_s1)
    assert e_s1 < e_s3, (e_s1, e_s3)
