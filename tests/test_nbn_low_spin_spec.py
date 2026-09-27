"""
Acceptance gates for specs/SPEC_nbn_low_spin.md: NbN CAS(14,14) low-spin sector.

Originally chem-dc7; corrected chem-g1i (2026-09-26) after chem-bbi found the committed
data/nbn_scf.chk describes a Nb-N separation of 3.7255 A, not the 2.25 A reconstruction chem-dc7
gated against. All energies/margins below are for the REAL, vendored geometry.

Pre-registered (specs/BACKLOG.md "Is NbN's flagged 'hard multireference benchmark' actually
hard?"): KILL the hard-benchmark framing if the (7,7) sector shows dw(D=300) < 1e-7 AND per-D
spread < 1e-5 Ha; CONFIRM it if dw > 1e-5 OR |E_A' - E_B'| > 0.1 mHa. On the real geometry, at
cheap D (<=300) NEITHER condition fires -- a genuine gap in the pre-registered design (G4 below).
The definitive verdict (soft, converging by D~400-1200) is in SPEC_nbn_low_spin.md Table 2, not
re-run every cycle. Which sector is the converged ground? The answer (exact FCI, 1.18e7
determinants) is NEITHER the committed high-spin (10,4) nor the (7,7) singlet but S=1 -- gated
here with DMRG so it re-runs in minutes (G5).

block2 runs in SU(2) mode with spin = na - nb: nelec=(7,7) is S=0, (8,6) is S=1, (10,4) is S=3.

pyscf + block2 ONLY (no qiskit imports; own process under `make gates` / run_gates.sh).
"""
import os

import numpy as np
import pytest

from hybrid_quantum_solver.dmrg_reference import dmrg_available
from nbn_dmrg_reference import run_schedule

_CHK = "data/nbn_scf.chk"
pytestmark = [
    pytest.mark.skipif(not os.path.exists(_CHK), reason=f"{_CHK} missing (gitignored; "
                       "regenerate with `uv run python nbn_low_spin.py cif`)"),
    pytest.mark.skipif(not dmrg_available(), reason="block2 not installed"),
]

_CACHE = {}


def _run(tag, nelec):
    key = (tag, nelec)
    if key not in _CACHE:
        # Own scratch per (schedule, sector): run_gates.sh runs this file alongside
        # test_nbn_dmrg_reference_spec.py, whose A'/B' would otherwise share the directory.
        suffix = "_lowspin_" + ("scf" if nelec is None else f"{nelec[0]}{nelec[1]}")
        _CACHE[key] = run_schedule(tag, n_threads=2, nelec=nelec, scratch_tag=suffix)
    return _CACHE[key]


def test_G4_low_spin_cheap_D_is_pre_registered_inconclusive():
    """chem-g1i correction: on the real geometry, the pre-registered cheap-D (<=300) KILL and
    CONFIRM conditions BOTH fail to fire -- a genuine gap in the pre-registered design, not a
    threshold nudge (measured dw(300) ~1.05e-7 for both schedules, just outside the 1e-7 KILL
    floor; per-D spread ~3-6e-5 Ha, also just outside the 1e-5 Ha KILL floor; |E_A'-E_B'| ~4.2e-6
    Ha, far under the 1e-4 Ha CONFIRM floor). The definitive verdict -- soft, converging by
    D~400-1200 -- comes from the headline schedules recorded in SPEC_nbn_low_spin.md Table 2, not
    re-run here (~9-17 min each); this gate pins the cheap-D near-miss as a regression check."""
    a, b = _run("A'", (7, 7)), _run("B'", (7, 7))
    w_a = {d: w for d, w, _ in a.per_D}
    w_b = {d: w for d, w, _ in b.per_D}
    e_a = np.array([e for _, _, e in a.per_D])
    killed = w_a[300] < 1e-7 and (e_a.max() - e_a.min()) < 1e-5
    confirmed = w_a[300] > 1e-5 or abs(a.energy - b.energy) > 1e-4
    assert not killed, (w_a, e_a)
    assert not confirmed, (w_a[300], a.energy, b.energy)
    # order-of-magnitude regression pins on the measured near-miss, not re-derived "safe" values
    assert 1e-8 < w_a[300] < 1e-6, w_a[300]
    assert 1e-8 < w_b[300] < 1e-6, w_b[300]
    assert abs(a.energy - b.energy) < 1e-4, (a.energy, b.energy)


def test_G5_committed_sector_is_not_the_cas_ground():
    """chem-g1i correction: on the real geometry the S=1 sector (nelec=(8,6)) lies BELOW the
    committed S=3 (10,4) reference by 1.538 mHa (DMRG A'; exact FCI agrees at 1.537 mHa) -- not
    the 15.1 mHa measured against chem-dc7's wrong-bond-length reconstruction. Margin lowered from
    5 mHa to 1 mHa: ~35% headroom below the measured gap, ~1000x above the ~1 uHa DMRG/FCI
    cross-check noise floor (SPEC_nbn_low_spin.md Table 1). The singlet (7,7) still lies above
    both."""
    e_s1 = _run("A'", (8, 6)).energy
    e_s3 = _run("A'", None).energy          # None = the SCF's own split, (10, 4)
    e_s0 = _run("A'", (7, 7)).energy
    assert e_s1 < e_s3 - 1e-3, (e_s1, e_s3)
    assert e_s0 > e_s3, (e_s0, e_s3)
