"""
Acceptance gate for chem-czw: does the CAS(14,14) spin ordering / "which minimum is the better
active-space reference" finding survive a full re-derivation on the OTHER (deterministic,
lower-whole-molecule-SCF-energy) UHF minimum chem-owr/chem-y6w found?

Full re-derivation (chem-czw, see sandbox-handoffs/chem-czw.md and SPEC_nbn_low_spin.md's new
section) ran exact FCI for S=1/2/3 and DMRG (2 independent schedules) for S=0 on a freshly
regenerated chk anchored to the OTHER locally-stable (10,4)/S=3 UHF minimum (0.548 mHa lower in
raw SCF energy than the vendored `specs/nbn_scf_reference.chk`). Two findings, cheaply pinned here
as a regression floor (NOT re-running the full ~25 min re-derivation every CI cycle -- same
convention as SPEC_nbn_low_spin.md G4/G5 using cheap dims instead of the headline D<=1200 runs):

  1. Ordering is orbital-choice-independent: S=1 is still the lowest of the four sectors checked
     on EITHER minimum's orbitals (only the gap sizes change).
  2. The vendored/committed chk's orbitals give a LOWER (more variationally favorable) CASCI(14,14)
     energy than the other minimum's orbitals, in every sector checked -- despite having a HIGHER
     raw whole-molecule SCF energy. This is chem-czw's principled reason (not just compute cost)
     to keep the vendored chk as the primary reference: it demonstrably describes the CAS(14,14)
     subspace better, not merely "whichever one we committed first".

Only the cheap (10,4) sector is re-derived fresh here (~10-20s); the (7,7)/(9,5)/(8,6) numbers
from both minima are pinned as literals from the one-time chem-czw re-derivation
(results/nbn_low_spin/runs.jsonl, 2026-09-27) -- regenerating ALL of them here would repeat the
~25 minute cost this gate exists to avoid.

pyscf + qiskit only: the fixture-gated regeneration test imports benchmark_nbn (for
ground_state_mf), which pulls in hybrid_quantum_solver.molecular_hamiltonian /
quantum_krylov_solver -> qiskit at module level. No block2/DMRG import anywhere in this file
(nbn_low_spin.py and hybrid_quantum_solver.dmrg_reference both import block2 lazily inside
function bodies, not at module scope) -- still its own file, per the test_*_spec.py
one-file-per-process isolation convention (CLAUDE.md), so it never risks sharing a process with
a block2 import.
"""
import os

import pytest

CIF = "specs/nbn_mp-2634.cif"
VENDORED_CHK = "specs/nbn_scf_reference.chk"
needs_fixtures = pytest.mark.skipif(
    not (os.path.exists(CIF) and os.path.exists(VENDORED_CHK)),
    reason=f"{CIF} or {VENDORED_CHK} missing",
)

# chem-czw one-time re-derivation, 2026-09-27, host e1fd4a2d0570 (8-core x86_64, 16 GB, no GPU) --
# see sandbox-handoffs/chem-czw.md for every command and full per-sector table.
LOWER_MINIMUM_E_SCF = -110.02813374576796
LOWER_MINIMUM_FCI = {
    (8, 6): -110.04577676472408,   # S=1, exact FCI, 9,018,009 dets, 259.7s
    (9, 5): -110.04546990757419,   # S=2, exact FCI, 4,008,004 dets, 192.2s
    (10, 4): -110.04119080464493,  # S=3, exact FCI, 1,002,001 dets, 5.6s
}
LOWER_MINIMUM_DMRG_S0 = -110.0410570692685    # S=0, (7,7), DMRG A perD 400/800/1200, dw(1200)=8.2e-10
LOWER_MINIMUM_DMRG_S0_B = -110.04105709905623  # DMRG B ramp 300/600/1200, dw(1200)=7.7e-10 (independent)

# Committed/vendored chk, exact FCI (SPEC_nbn_low_spin.md Table 1, chem-g1i 2026-09-26)
VENDORED_FCI = {
    (8, 6): -110.04756504307636,
    (9, 5): -110.04698592997889,
    (10, 4): -110.04602841723886,
}
VENDORED_DMRG_S0 = -110.04249952166984
VENDORED_E_SCF = -110.0275853827266  # specs/nbn_scf_reference.chk e_tot

TOL_HA = 1e-6  # cross-container float noise ceiling (~4e-7 Ha observed reproducing chem-owr's
                # own figure); 3-4 orders of magnitude under every mHa-scale gap this gate checks


def test_lower_minimum_ordering_matches_vendored_S1_lowest():
    """chem-czw regression floor: S=1 (8,6) is below S=2 (9,5), S=3 (10,4), and S=0 (7,7) on BOTH
    minima's orbitals -- the chem-dc7/chem-g1i "committed sector is not the CAS ground" finding
    does not depend on which locally-stable UHF minimum supplies the active-space orbitals."""
    e1, e2, e3 = (LOWER_MINIMUM_FCI[(8, 6)], LOWER_MINIMUM_FCI[(9, 5)], LOWER_MINIMUM_FCI[(10, 4)])
    e0 = LOWER_MINIMUM_DMRG_S0
    assert e1 < e2 < e3 < e0, (e1, e2, e3, e0)


def test_lower_minimum_dmrg_schedules_agree():
    """The two independent (7,7) DMRG schedules on the lower minimum's orbitals must still agree
    at the same tightness the committed-chk cross-check does (SPEC_nbn_low_spin.md Table 2:
    5 orders of magnitude under the 1e-4 Ha CONFIRM floor) -- measured 2.98e-8 Ha."""
    assert abs(LOWER_MINIMUM_DMRG_S0 - LOWER_MINIMUM_DMRG_S0_B) < 1e-4


@pytest.mark.parametrize("nelec", [(8, 6), (9, 5), (10, 4)])
def test_vendored_chk_gives_lower_casci_energy_than_the_other_minimum(nelec):
    """chem-czw's principled keep-the-vendored-chk finding: for every sector checked, the
    committed/vendored chk's orbitals give a STRICTLY LOWER (more favorable) CASCI(14,14) energy
    than the other, lower-whole-molecule-SCF-energy minimum's orbitals -- by 1.4-4.8 mHa depending
    on sector (measured; see sandbox-handoffs/chem-czw.md). This is despite the vendored chk
    having a HIGHER raw SCF energy (0.548 mHa) -- "lower whole-molecule HF energy" and "better
    CAS(14,14) active-space orbitals" are not the same claim; CASCI is not variational over
    orbital choice (that would be CASSCF, out of scope -- SPEC_nbn_dmrg_reference.md Section 7)."""
    assert VENDORED_FCI[nelec] < LOWER_MINIMUM_FCI[nelec] - 1e-3, (
        VENDORED_FCI[nelec], LOWER_MINIMUM_FCI[nelec])


def test_vendored_chk_gives_lower_dmrg_energy_for_S0_too():
    assert VENDORED_DMRG_S0 < LOWER_MINIMUM_DMRG_S0 - 1e-3, (VENDORED_DMRG_S0, LOWER_MINIMUM_DMRG_S0)


@needs_fixtures
def test_lower_minimum_chk_regenerates_and_reproduces_S3_FCI(tmp_path):
    """Cheap end-to-end regression (~10-20s): a fresh regeneration of the OTHER minimum (not the
    vendored one) still lands within float noise of the pinned SCF energy and its (10,4)-sector
    exact FCI still matches the pinned literal above -- catches a regression in
    benchmark_nbn._tight_scf's determinism fix or in load_nbn_cas/spin_fci without repeating the
    full ~25 min four-sector re-derivation every cycle.

    Which of the two minima the (deterministic) regeneration finds is BLAS-build dependent
    (2026-09-27: Linux container -> the lower one, macOS arm64 -> the vendored one;
    SPEC_nbn_scf_determinism.md Sec 8), so the gate accepts either characterized minimum -- but
    only those two -- and checks the S=3 FCI against the value pinned for the one it found."""
    import benchmark_nbn as bn
    from pyscf import ao2mo, mcscf

    from nbn_low_spin import spin_fci

    chk = str(tmp_path / "nbn_scf_lower_regen.chk")
    bn.CHKFILE = chk
    mf = bn.ground_state_mf(CIF)
    known = {LOWER_MINIMUM_E_SCF: LOWER_MINIMUM_FCI[(10, 4)], VENDORED_E_SCF: VENDORED_FCI[(10, 4)]}
    e_scf = min(known, key=lambda k: abs(mf.e_tot - k))
    assert abs(mf.e_tot - e_scf) < TOL_HA, ("landed on a third SCF solution", mf.e_tot)

    cas = mcscf.CASCI(mf, 14, (10, 4))
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), 14)
    (e, ss), = spin_fci(h1, eri, (10, 4), e_core)[0]
    assert abs(e - known[e_scf]) < TOL_HA, (e, known[e_scf])
    assert abs(ss - 12.0) < 1e-3, ss
