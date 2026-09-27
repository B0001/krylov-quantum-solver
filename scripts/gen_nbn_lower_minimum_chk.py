#!/usr/bin/env python3
"""
One-off, chem-czw: materialize the deterministic (pyscf.lib.num_threads(1)-pinned) LOWER-energy
UHF minimum benchmark_nbn._tight_scf actually finds on a fresh regeneration from specs/nbn_mp-2634.cif
-- the "other" locally-stable (10,4)/S=3 solution chem-owr/chem-y6w found, distinct from the
vendored specs/nbn_scf_reference.chk. Writes the resulting chk to the path given on argv[1] (does
NOT touch data/nbn_scf.chk or specs/nbn_scf_reference.chk) and asserts the energy matches the
value pinned by tests/test_nbn_scf_determinism_spec.py's G1/G2 gates, as a sanity check that this
run landed on the same, already-characterized minimum rather than some third solution.

Run in its own process (imports benchmark_nbn -> qiskit-nature; do not follow with a block2 import
in the same process -- CLAUDE.md OpenMP isolation note).
"""
import sys

import benchmark_nbn as bn

EXPECTED_E_TOT = -110.02813374576796  # tests/test_nbn_scf_determinism_spec.py G1/G2; chem-owr
# G1/G2 pin bit-identical e_tot only WITHIN one container/BLAS build (SPEC_nbn_scf_determinism.md
# R1); across a fresh uv sync on a different container, expect last-bit noise, not exact equality.
# 1e-9 Ha is 5 orders of magnitude tighter than the 0.548 mHa SCF gap between the two minima --
# generous headroom to distinguish "same minimum, float noise" from "landed on a third solution".
TOL_HA = 1e-9
# NOTE (2026-09-27): only the Linux container's BLAS build lands on this minimum; macOS arm64's
# deterministic path lands on the VENDORED one (-110.02758538272664) and this assert fires there.
# See SPEC_nbn_scf_determinism.md Sec 8.

if __name__ == "__main__":
    out_chk = sys.argv[1]
    bn.CHKFILE = out_chk
    mf = bn.ground_state_mf("specs/nbn_mp-2634.cif")
    print(f"e_tot = {mf.e_tot!r}")
    diff = abs(mf.e_tot - EXPECTED_E_TOT)
    assert diff < TOL_HA, (
        f"landed on a DIFFERENT minimum than chem-owr/chem-y6w's characterized one: "
        f"got {mf.e_tot!r}, expected {EXPECTED_E_TOT!r} (diff {diff:.3e} Ha >= {TOL_HA:.0e} Ha)"
    )
    print(f"OK: matches the deterministic lower-energy minimum (diff {diff:.3e} Ha), "
          f"written to {out_chk}")
