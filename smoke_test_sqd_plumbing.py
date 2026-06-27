"""Smoke test: exercise the REAL shipped functions on systems with known answers."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # find the script beside this file
import run_nbn_sqd_sweep as R   # the actual deliverable
from pyscf import gto, scf

def build_mf(atom, basis="sto-3g", spin=0):
    mol = gto.M(atom=atom, basis=basis, spin=spin)
    mf = scf.RHF(mol) if spin == 0 else scf.ROHF(mol)
    mf.verbose = 0; mf.kernel()
    return mf

rng = np.random.default_rng(0)
print("="*70)

# --- Test 1: H2, CAS(2,2) = full space. SQD over all determinants == CASCI. ---
mf = build_mf("H 0 0 0; H 0 0 0.74")
hcore, eri, e_core, nelec, casci = R.integrals_for_spin(mf, 2, 2, 0)
norb = hcore.shape[0]
ba = R.generate_bit_array_uniform(20_000, norb*2, rand_seed=rng)
best, max_dim, iters, _ = R.run_sqd_for_sector(hcore, eri, e_core, nelec, norb, ba,
                                               samples_per_batch=100, spin_sq=0.0, rng=rng)
err = abs(best - casci)*1e3
print(f"[H2  CAS(2,2)] CASCI={casci:.8f}  SQD={best:.8f}  |Δ|={err:.4f} mHa  "
      f"frame_ok={err < 1e-3}")
assert err < 1e-3, "FRAME MISALIGNED: SQD != CASCI on full space"
assert best >= casci - 1e-6, "VARIATIONAL VIOLATION"
print(f"            verdict={R.validate_row(best, casci)}  (e_core={e_core:.6f} Ha added once)")

# --- Test 2: deliberately BREAK the frame (drop e_core) -> verdict must catch it ---
broken = best - e_core   # what you'd get if e_core were discarded (the 507 Ha class of bug)
print(f"\n[frame-break sim] dropped e_core -> 'energy'={broken:.6f}  "
      f"verdict={R.validate_row(broken, casci)}")
assert R.validate_row(broken, casci) in ("FRAME_ERROR", "ABOVE_TOL")

# --- Test 3: H4 chain, CAS(4,4). Sweep subspace, check the trend gate on REAL output. ---
mf = build_mf("H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0")
hcore, eri, e_core, nelec, casci = R.integrals_for_spin(mf, 4, 4, 0)
norb = hcore.shape[0]
ba = R.generate_bit_array_uniform(40_000, norb*2, rand_seed=rng)
rows = []
print(f"\n[H4  CAS(4,4)] CASCI={casci:.8f} Ha")
for spb in [10, 30, 80, 200]:
    best, max_dim, iters, _ = R.run_sqd_for_sector(hcore, eri, e_core, nelec, norb, ba,
                                                   samples_per_batch=spb, spin_sq=0.0, rng=rng)
    d = abs(best - casci)*1e3
    rows.append({"delta_mHa": d, "subspace_dim": max_dim})
    print(f"   spb={spb:>3}  SQD={best:.8f}  Δ={d:7.3f} mHa  dim={max_dim}  "
          f"[{R.validate_row(best, casci)}]")
    assert best >= casci - 1e-6, "VARIATIONAL VIOLATION on H4"
print(f"   -> sector trend: {R.analyze_sector_trend(rows)}")
print("="*70)
print("ALL SMOKE TESTS PASSED")
