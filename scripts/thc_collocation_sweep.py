#!/usr/bin/env python3
"""
chem-5oj corroboration sweep for specs/SPEC_thc_collocation.md.

Independently reproduces (not just trusts) the numbers quoted in the G2/G4 docstrings of
tests/test_thc_collocation_spec.py and in SPEC_thc_collocation.md Sec 8 R2:

  1. A 12-seed x 8000-LM-evaluation (restarts=1 per seed, so each seed gets the full eval budget
     instead of splitting it across 4 restarts) reconstruction-error sweep for `fit_thc` at
     LiH/STO-3G's matched rank M=21 -- corroborates the G2 "fit_thc does not reach <1e-6
     reconstruction" finding beyond the single CI-gated seed=0/restarts=4/max_nfev=4000 point.
  2. A 5-seed beats_df_lambda / ratio_to_random spot check at fit_thc's CI-gate defaults
     (restarts=4, max_nfev=4000) -- corroborates the G4 "beats df_lambda at 4/5 seeds, >=5x below
     random at all 5" claim.
  3. A single, NOT CI-gated norb=7 (H2O/STO-3G full space) run at fit_thc's CI-gate defaults, to
     report (not gate) the ~600-parameter cost case this bead's cost note names.

Pinned seeds throughout (never claim global optimality from a single LM run -- see
specs/SPEC_thc_collocation.md Sec 2 and Sec 8 R2). PySCF/NumPy/SciPy only, no block2 -- safe to run
in the same process as the rest of this sweep.
"""
import time

import numpy as np
from pyscf import ao2mo, gto, mcscf, scf

from df_factorization import double_factorize, df_lambda
from thc_factorization import reconstruct_thc, tensor_hypercontraction, thc_lambda, thc_rank
from lambda_ladder import fit_thc


def build_system(atom, basis="sto-3g"):
    mol = gto.M(atom=atom, basis=basis)
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    norb = mol.nao_nr()
    na = nb = mol.nelectron // 2
    cas = mcscf.CASCI(mf, norb, (na, nb))
    cas.verbose = 0
    cas.kernel()
    h1, ecore = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), norb)
    return dict(h1=h1, eri=eri, norb=norb, nelec=(na, nb), ecore=ecore, e_fci=cas.e_tot)


def sweep_reconstruction_error(eri, norb, M, seeds, max_nfev):
    """One LM run per seed, restarts=1 (full eval budget per seed, no internal best-of-restarts)."""
    out = []
    for s in seeds:
        t0 = time.time()
        X, Z = fit_thc(eri, norb, M, restarts=1, seed=s, max_nfev=max_nfev, return_factors=True)
        dt = time.time() - t0
        err = np.linalg.norm(reconstruct_thc(X, Z) - eri)
        out.append((s, err, dt))
    return out


def spot_check_beats_df(eri, norb, h1, M, seeds, lam_df, lam_r, restarts, max_nfev):
    out = []
    for s in seeds:
        t0 = time.time()
        X, Z = fit_thc(eri, norb, M, restarts=restarts, seed=s, max_nfev=max_nfev,
                        return_factors=True)
        dt = time.time() - t0
        err = np.linalg.norm(reconstruct_thc(X, Z) - eri)
        lam = thc_lambda(X, Z, h1)
        out.append((s, err, lam, lam / lam_r, lam < lam_df, dt))
    return out


def main():
    print("=" * 90)
    print("1) LiH/STO-3G norb=6 M=21 -- 12-seed x 8000-nfev (restarts=1) reconstruction-error sweep")
    print("   corroborates test_thc_collocation_spec.py::test_G2 (fit_thc fails <1e-6 precondition)")
    print("=" * 90)
    lih = build_system("Li 0 0 0; H 0 0 1.6")
    M = thc_rank(lih["norb"])
    assert lih["norb"] == 6 and M == 21

    t_sweep0 = time.time()
    sweep = sweep_reconstruction_error(lih["eri"], lih["norb"], M, seeds=range(12), max_nfev=8000)
    t_sweep_total = time.time() - t_sweep0
    errs = [e for _, e, _ in sweep]
    for s, e, dt in sweep:
        print(f"  seed={s:2d}  recon_err={e:.4e}  time={dt:.1f}s")
    print(f"  min={min(errs):.4e}  max={max(errs):.4e}  mean={np.mean(errs):.4e}")
    print(f"  total wall time = {t_sweep_total:.1f}s ({t_sweep_total/60:.2f} min)")
    all_above_1e6 = all(e > 1e-6 for e in errs)
    print(f"  ALL 12 seeds stay above 1e-6 (G2 precondition never met): {all_above_1e6}")

    print()
    print("=" * 90)
    print("2) LiH/STO-3G norb=6 M=21 -- 5-seed beats_df_lambda / ratio_to_random spot check")
    print("   (fit_thc at CI-gate defaults: restarts=4, max_nfev=4000)")
    print("   corroborates test_thc_collocation_spec.py::test_G4")
    print("=" * 90)
    leaves, _, full_rank = double_factorize(lih["eri"], lih["norb"])
    lam_df = df_lambda(leaves, lih["h1"], lih["norb"])
    chi_r, zeta_r = tensor_hypercontraction(lih["eri"], lih["norb"], n_thc=M, seed=0)
    lam_r = thc_lambda(chi_r, zeta_r, lih["h1"])
    print(f"  df_lambda={lam_df:.4f}  random_lambda={lam_r:.4f} (seed=0)")

    spot = spot_check_beats_df(lih["eri"], lih["norb"], lih["h1"], M, seeds=range(5),
                                lam_df=lam_df, lam_r=lam_r, restarts=4, max_nfev=4000)
    n_beats = 0
    for s, err, lam, ratio, beats, dt in spot:
        n_beats += int(beats)
        print(f"  seed={s}  recon_err={err:.4e}  lambda={lam:.4f}  ratio_to_random={ratio:.4f}"
              f"  beats_df={beats}  time={dt:.1f}s")
    print(f"  beats_df_lambda at {n_beats}/5 seeds")
    ratios = [r for _, _, _, r, _, _ in spot]
    print(f"  ratio_to_random range: [{min(ratios):.4f}, {max(ratios):.4f}]  "
          f"(all <0.2, i.e. >=5x below random: {all(r < 0.2 for r in ratios)})")

    print()
    print("=" * 90)
    print("3) H2O/STO-3G full space, norb=7 -- single NOT-CI-gated run at fit_thc CI-gate defaults")
    print("   (this bead's cost note: '~600 Levenberg-Marquardt parameters at norb=7')")
    print("=" * 90)
    h2o = build_system("O 0 0 0.117; H 0 0.757 -0.467; H 0 -0.757 -0.467")
    M7 = thc_rank(h2o["norb"])
    n_params = h2o["norb"] * M7 + M7 * (M7 + 1) // 2
    print(f"  norb={h2o['norb']}  M={M7}  n_LM_params={n_params}")
    leaves7, _, full_rank7 = double_factorize(h2o["eri"], h2o["norb"])
    lam_df7 = df_lambda(leaves7, h2o["h1"], h2o["norb"])
    chi_r7, zeta_r7 = tensor_hypercontraction(h2o["eri"], h2o["norb"], n_thc=M7, seed=0)
    lam_r7 = thc_lambda(chi_r7, zeta_r7, h2o["h1"])
    t0 = time.time()
    X7, Z7 = fit_thc(h2o["eri"], h2o["norb"], M7, restarts=4, seed=0, max_nfev=4000,
                      return_factors=True)
    dt7 = time.time() - t0
    err7 = np.linalg.norm(reconstruct_thc(X7, Z7) - h2o["eri"])
    lam7 = thc_lambda(X7, Z7, h2o["h1"])
    print(f"  df_lambda={lam_df7:.4f}  random_lambda(seed=0)={lam_r7:.4f}")
    print(f"  fit_thc seed=0 restarts=4 max_nfev=4000: time={dt7:.1f}s  recon_err={err7:.4e}  "
          f"lambda={lam7:.4f}  ratio_to_random={lam7/lam_r7:.4f}  beats_df={lam7 < lam_df7}")
    print("=" * 90)


if __name__ == "__main__":
    main()
