# SPEC: Hₙ to larger n, done right (adequate-D ramp + bulk per-site estimator)

**Status:** PARTIAL. CI gates G1–G3 PASS (`make gates`, 2026-06-29); `bulk_per_site_energy`
is merged. The headline large-n ramp (D=400/800/1600, 4 threads, canonical RHF orbitals) **was run**
(`logs/run_hchain_headline.sh`, 2026-07-04): it took ~21 h to finish n=8..20 (04:50 on 07-04 to 01:53
on 07-05, `logs/hchain_pipeline.log` on the laptop, not tracked), and the n=22/24 stage was killed
after ~19.5 h without a row. Its largest points are out of regime (see bead chem-oxm), so the
definition of done in §5 is still **not met**. An earlier attempt (2026-06-29, n=16..30, 8 threads)
aborted when block2 ran out of stack memory before the first row. **§10 (2026-09-24):** the cost was
the orbital basis. In Loewdin site orbitals, n=20 converges at D=50 in 20 s to within 5.4 µHa of the
D=400 result; canonical orbitals at D=400 are still 6.4 mHa above it after 9.5 min.

> Supersedes [`SPEC_hchain_largen.md`](SPEC_hchain_largen.md), which was **killed**: the cheap
> ramp at D=100/200/400 truncated too hard as chain entanglement grew — discarded-weight stderr
> ballooned to ~5 mHa, the `dweight` extrapolation fell back to `invD` by n=30, and leave-one-out
> reached 1.07 mHa/atom (gate < 0.1). The protocol was fine; the bond dims were not.

---

## 1. Goal

The Hₙ thermodynamic-limit per-atom energy can be pinned to a **leave-one-out shift < 0.1 mHa/atom**
at large n by either (a) the single-ramp DMRG protocol at **adequate bond dimension** (D ≈ 400/800/1600),
or (b) a **bulk per-site estimator** `e_bulk = (E(n) − E(n−Δ))/Δ` that cancels the open-chain surface
term directly. The claim is false if neither route reaches 0.1 mHa/atom while staying in the
discarded-weight regime (`method == "dweight"`).

## 2. Background and honest framing

- Builds directly on the validated `protocol="ramp"` extrapolation ([`SPEC_singleramp.md`]) and the
  TDL fit ([`SPEC_hchain_tdl.md`], e∞ = −0.539967 ± 0.000107 Ha/atom from n ≤ 16).
- **What we can claim if gates pass:** a tighter, better-controlled TDL estimate for the
  minimal-basis Hₙ chain, with two independent extrapolation routes agreeing — and an honest
  account of the bond dimension required to stay in regime at large n.
- **What we cannot claim:** novelty. This is a minimal-basis (sto-6g), open-boundary **model** at a
  fixed geometry that reproduces benchmark physics (cf. Motta et al., PRX 7, 031059, 2017). The
  quantum/statevector solver plays no role at these n. No quantum advantage.

## 3. Approach

For each n, get `E(n)` from `dmrg_energy_extrapolated(..., protocol="ramp", bond_dims=(400,800,1600))`,
which must remain in the discarded-weight regime (`method == "dweight"`). Two extrapolations to the
bulk:

1. **Surface-term fit** — `thermodynamic_limit_fit(ns, E(n)/n)` → e∞ (existing).
2. **Bulk per-site estimator** — `e_bulk = (E(n) − E(n−Δ))/Δ`. For an open chain
   `E(n) ≈ n·e∞ + c` (a constant two-end surface term), so the difference quotient removes `c`
   exactly; e_bulk → e∞ from a different direction and is the cross-check on the fit.

**Reference:** exact FCI where tractable (n ≤ 14); the two extrapolation routes referee each other
at large n where no FCI exists. Driver: `benchmark_hchain_tdl.py` (already resumable, already
supports `--protocol ramp --ns --bond-dims`).

## 4. Public interface

Prefer composing validated primitives; the only new code is the bulk estimator + a summary line.

```
hybrid_quantum_solver.dmrg_reference.bulk_per_site_energy(ns, totals, *, step=None) -> float
    # (E(n) − E(n−Δ))/Δ from the two largest available n (Δ = step or the largest gap); Ha/atom
benchmark_hchain_tdl.py  --protocol ramp --ns 8,...,30 --bond-dims 400,800,1600
    # -> data/hchain_tdl.csv + printed e∞ (1/n fit), e_bulk, and leave-one-out shift
```

## 5. Acceptance criteria (validation gates)

In `tests/test_hchain_largen2_spec.py` (test-first). G1 is pure/instant; G2–G3 need block2 and run
in their own process via `make gates`.

- **G1 — estimator inverts the surface term (instant, no DMRG).** On synthetic open-chain totals
  `E(n) = n·e∞ + c`, `bulk_per_site_energy` returns `e∞` to < 1e-9, and agrees with
  `thermodynamic_limit_fit` on the same ideal data to < 1e-9.
- **G2 — two routes agree on real chains (cheap CI proxy of "done").** For H_n at converged dims
  `(80,160,300)`, n ∈ {8,10,12,14}: `|e_bulk − e∞_fit| < 1e-3` Ha (1 mHa/atom).
  **Finding (revised from 0.1 mHa/atom):** at n ≤ 14 the bulk difference quotient still sits on the
  curved part of `E(n)`, so the two routes agree only to **0.51 mHa/atom** (`e_bulk = −0.540501`
  vs `e∞_fit = −0.539987`, the latter matching the n ≤ 16 TDL reference −0.539967 to 0.02 mHa).
  0.1 mHa/atom agreement between the routes is therefore a **large-n** claim, verified by the driver
  (below), not a small-n one. G2 is loosened to 1 mHa/atom and kept as a "estimator is sane and on
  the right side" guard; the headline number lives in the definition of done.
- **G3 — regime guard.** At n=12, converged dims, the extrapolation stays
  `res.method == "dweight"` (no `invD` fallback — the failure mode that killed the cheap spec).

**Definition of done.** The gate that actually settles the claim is **leave-one-out shift < 0.1
mHa/atom** on the *large-n* ramp run (n up to ~30, D=400/800/1600). That run is minutes-to-an-hour,
so it is **driver-level**, recorded in the PR with the CSV — not a CI unit gate. G1–G3 are the cheap
falsifiers that gate CI; G2 is the small-n proxy for the headline number.

> If G2 proves unsatisfiable at small n (the two routes may only converge to 0.1 mHa/atom at larger
> n), **loosen G2's tolerance and record the n at which they agree** — that crossover is the finding
> (cf. `SPEC_hchain_tdl.md` G1).

## 6. Implementation plan (test-first)

1. Write `tests/test_hchain_largen2_spec.py` encoding G1–G3 (initially RED — `bulk_per_site_energy`
   does not exist yet).
2. Add `bulk_per_site_energy` to `dmrg_reference.py` (≈5 lines) and a summary line to
   `benchmark_hchain_tdl.py`. Reuse `dmrg_energy_extrapolated`/`thermodynamic_limit_fit` unchanged.
3. `make gates` to green. Then run the large-n driver and record the LOO shift in the PR.

## 7. Out of scope

- Larger basis sets / CBS extrapolation (a separate spec).
- Periodic boundary conditions; 2D/3D lattices.
- Any quantum-hardware or statevector path at these n.

## 8. Caveats and risks

- **R1 — adequate D is still not enough at very large n.** Mitigation: the bulk estimator (route b)
  degrades more gracefully than the 1/n fit; if both miss 0.1 mHa/atom, that is a recorded finding
  (D needed grows with n), not a silent pass — exactly how the cheap spec was killed.
  **Update (§10):** on canonical orbitals this was largely an artefact of the orbital basis. In Loewdin
  site orbitals, n=20 converges at D≈100–200.
- **R2 — cost.** D=1600 ramp at n≈30 is ~an hour; the driver is resumable and the CI gates stay at
  small n. Large-n runs belong on a faster/GPU node.
- Honest limitation: minimal-basis open-chain model; reproduces, does not extend, known physics.

## 9. Deliverables

- `hybrid_quantum_solver/dmrg_reference.py` — `bulk_per_site_energy`.
- `benchmark_hchain_tdl.py` — bulk estimate + leave-one-out shift in the printed summary.
- `tests/test_hchain_largen2_spec.py` — gates G1–G3.
- Results summary (with §2/§7 caveats + the large-n LOO number) in the PR description.
- §10 (chem-rhw): `benchmark_hchain_tdl.py --localize`, `probe_hchain_localize.py`, localized gates
  in `tests/test_hchain_tdl_spec.py`, and the `fci_energy` convergence check.

## 10. Result: localized site orbitals (bead chem-rhw, 2026-09-24)

**Hypothesis (inferred from the code, not measured before this):** every H-chain DMRG run above used
canonical, delocalized RHF orbitals (`integrals()` passed `mf.mo_coeff` straight to block2 with no
reordering). For a 1D chain that makes the MPS entanglement grow with n, which would explain the ~21 h
cost and R1's "D needed grows with n".

**Change.** `benchmark_hchain_tdl.py --localize` uses Loewdin site orbitals C = S^(-1/2) over the full
STO-6G basis (one AO per H, so already in chain order). This is a unitary rotation of the full space,
so FCI is invariant. The default is still canonical, and localized rows go to a separate CSV.

**Validation against FCI** (`tests/test_hchain_tdl_spec.py`, new gates next to G1):

| n | E_FCI canonical (Ha) | E_FCI localized (Ha) | difference |
|---|---|---|---|
| 8  | −4.345079400958 | −4.345079400919 | 3.9e-11 |
| 10 | −5.424385375303 | −5.424385375235 | 6.8e-11 |
| 12 | −6.504226955890 | −6.504226955828 | 6.2e-11 |

- Gate `test_G1_localized_fci_invariant`: \|ΔE\| < 1e-8 at n = 8, 10, 12. PASS.
- Gate `test_G1_localized_extrapolation_matches_fci`: CONV_DIMS (80,160,300) extrapolation at
  n = 10, 12 has regime ≠ "uncontrolled" and \|E − E_FCI\| < G1_TOL = 2e-4. It also checks that the
  largest-D energy is within 1e-6 of FCI. PASS. Measured: regime = "converged" at both n (every
  discarded weight is below the 1e-8 floor). The extrapolated error is 3.5e-10 (n=10) and 1.9e-9
  (n=12). Canonical orbitals at the same D give 5.3e-6 and 1.8e-4, with regime = "truncation".
- Discarded weight at D=80: 1.2e-10 localized vs 2.1e-4 canonical (n=10), and 8.7e-10 vs 6.3e-4 (n=12).
- **Converged (not extrapolated) localized DMRG at n=12, D=1000** (schedule 250/500/1000, 20 sweeps):
  E − E_FCI = **+2.2e-11 Ha**.

**Finding on the way (a latent bug).** `fci_energy` called PySCF's Davidson with its defaults
(100 iterations) and **never checked convergence**. In the site basis the diagonal initial guess is
poor, so it silently returned energies 0.28 mHa (n=10) and **91 mHa** (n=12) too high. It now uses
`max_cycle=1000` and raises if it does not converge. Without that fix, the invariance gate fails.
That is the falsifier doing its job. Other scripts still call `fci.direct_spin1.kernel` directly
without a convergence check (`lambda_ladder.py`, `df_factorization.py`, `thc_factorization.py`,
`shift_both_sides.py`, and two spec tests). They are canonical-basis today, but they carry the same
risk.

**Cost probe: n = 20, R = 1.8 bohr, STO-6G** (`probe_hchain_localize.py`). Both bases ran with an
identical schedule: D held fixed for up to 20 sweeps, noise 1e-4×4 / 1e-5×4 / 1e-6×4 / 0×8,
Davidson threshold 1e-10, and block2's default early stop (\|ΔE\| < 1e-8 on a noise-free sweep).
Other settings: random initial MPS seeded 1234, 4 threads, 6 GB block2 stack, one subprocess per run
(so peak RSS is that run's alone), a 90 min cap (not hit), and runs made one at a time on an idle
machine.
Hardware: Linux 6.18.44 x86_64 cloud VM, Intel Xeon @ 2.10 GHz, `nproc` = 4, 15 GiB RAM, no GPU.
n = 20 has no FCI (C(20,10)² ≈ 3.4e10 determinants), so **ΔE is measured against the lowest
(variational) energy obtained, localized D=400.** It is therefore a lower bound on each run's true
error.

| basis | D | E (Ha) | ΔE vs best (Ha) | discarded weight (final sweep) | sweeps | wall time | peak RSS |
|---|---|---|---|---|---|---|---|
| canonical | 400 | −10.819460843 | +6.4e-3 | 7.0e-4 | 20 (did not meet the stop criterion; ΔE/sweep still −3e-6) | 572 s | 2.65 GB |
| localized | 400 | −10.825878868 | 0 (reference) | 8.5e-14 (5e-11 on the noisy sweeps) | 13 | 290 s | 2.16 GB |
| localized | 200 | −10.825878868 | +5.3e-10 | 2.6e-11 | 13 | 84 s | 0.97 GB |
| localized | 100 | −10.825878811 | +5.7e-8 | 4.4e-9 | 15 | 38 s | 0.58 GB |
| localized | 50  | −10.825873428 | +5.4e-6 | 3.3e-7 | 15 | 20 s | 0.49 GB |

**Verdict: yes, localization cuts the D needed by at least 8×.** Localized D=50 (20 s) is already
~1000× more accurate than canonical D=400 (572 s, unconverged). At fixed accuracy, the time drop is
therefore at least ~29× at n=20 (572 s vs 20 s). At the same D=400, the discarded weight drops by 7 to 10 orders of
magnitude and the energy by 6.4 mHa. The accuracy-matched D and time ratios are lower bounds: the
canonical run never reached the localized runs' accuracy, so its true D for 1e-6 Ha was not found.
R1 ("D needed grows with n") and the ~21 h cost were largely an orbital-basis artefact. chem-oxm
should run with `--localize`, and does not need a bigger machine. This measures n = 20 at R = 1.8 bohr
only. Whether the needed D stays flat at n = 40, and how it depends on R, are chem-oxm and chem-1uu.
No TDL extrapolation is done here.
