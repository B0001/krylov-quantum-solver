"""
Acceptance gates for specs/SPEC_nbn_low_spin.md (bead chem-dc7): NbN CAS(14,14) low-spin sector.

Pre-registered (specs/BACKLOG.md "Is NbN's flagged 'hard multireference benchmark' actually
hard?"): KILL the hard-benchmark framing if the (7,7) sector shows dw(D=300) < 1e-7 AND per-D
spread < 1e-5 Ha; CONFIRM it if dw > 1e-5 OR |E_A' - E_B'| > 0.1 mHa. Then: which sector is the
converged ground? The answer found (exact FCI, 1.18e7 determinants) is NEITHER the committed
high-spin (10,4) nor the (7,7) singlet but S=1 -- gated here with DMRG so it re-runs in minutes.

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
        _CACHE[key] = run_schedule(tag, n_threads=2, nelec=nelec)
    return _CACHE[key]


def test_G4_low_spin_is_hard_preregistered():
    """(7,7) singlet, cheap schedules A'/B': the pre-registered KILL condition must NOT hold and
    the CONFIRM condition must (measured dw(300) = 4.8e-5 for both, spread 2.3 mHa)."""
    a, b = _run("A'", (7, 7)), _run("B'", (7, 7))
    w_a = {d: w for d, w, _ in a.per_D}
    e_a = np.array([e for _, _, e in a.per_D])
    killed = w_a[300] < 1e-7 and (e_a.max() - e_a.min()) < 1e-5
    assert not killed, (w_a, e_a)
    assert w_a[300] > 1e-5 or abs(a.energy - b.energy) > 1e-4, (w_a[300], a.energy, b.energy)


def test_G5_committed_sector_is_not_the_cas_ground():
    """The S=1 sector (nelec=(8,6)) lies well BELOW the committed S=3 (10,4) reference: exact FCI
    says by 15.1 mHa; the gate asks DMRG A' for > 5 mHa, and that the singlet lies above both."""
    e_s1 = _run("A'", (8, 6)).energy
    e_s3 = _run("A'", None).energy          # None = the SCF's own split, (10, 4)
    e_s0 = _run("A'", (7, 7)).energy
    assert e_s1 < e_s3 - 5e-3, (e_s1, e_s3)
    assert e_s0 > e_s3, (e_s0, e_s3)
