# SPEC: Hₙ to larger n, done right (adequate-D ramp + bulk per-site estimator)

**Status:** DONE (2026-09-24, bead chem-oxm; §11). The definition of done is met in Loewdin site
orbitals at n ≤ 40: leave-one-out = **0.037 mHa/atom** (< 0.1), with every point in regime.
**Headline e∞ = −0.540353 ± 0.000192 Ha/atom.** The bar is systematic-inclusive (the 8-fit
envelope plus the largest stderr). Motta's −0.540493 lies **inside** it (gap 0.140 mHa/atom =
0.73 bar). The bulk-vs-fit cross-check **fails** as a recorded finding: 0.232 mHa/atom vs the
0.1 target. The a + b/n form is biased at these n; the four a + b/n + c/n² fits all give
−0.540492 to −0.540496. The canonical-orbital history is kept below and in §10: the n ≤ 16
headline −0.539967 ± 0.000107 is **superseded**. The ~21 h canonical ramp to n=20 (2026-07-04,
killed at n=22/24) and the 2026-06-29 stack-memory abort were the orbital basis, not the
hardware (§10).

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
  TDL fit ([`SPEC_hchain_tdl.md`], e∞ = −0.539967 ± 0.000107 Ha/atom from n ≤ 16; **superseded**
  by §11: that ± was a fit stderr blind to finite-size bias).
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

## 11. Thermodynamic limit in localized orbitals (bead chem-oxm)

### 11.1 Pre-registered analysis (committed before any large-n number existed)

The analysis is code, `hchain_tdl_analysis.py`, and was committed before the run finished. It must
not change after the data are seen.

- **Data:** `benchmark_hchain_tdl.py --localize --protocol ramp --bond-dims 100,200,400`, at
  n ∈ {8, 10, 12, 16, 20, 24, 28, 32, 40} as far as the budget allows. Any n whose ladder is
  under-converged is rerun at 200/400/800. Every point records its `regime`.
- **Fits (8):** e(n) = a + b/n and e(n) = a + b/n + c/n², each over n_min ∈ {8, 12, 16, 20}, with
  n_max = the largest n reached.
- **Systematic-inclusive bar:** the envelope of `a` over the 8 fits. Headline = midpoint of
  [min a, max a], and bar = half the envelope width plus the largest single-fit stderr.
- **LOO:** the driver's definition: the all-n a + b/n fit minus the same fit with the largest n
  dropped. Gate: < 0.1 mHa/atom.
- **Bulk:** `bulk_per_site_energy` on the two largest n, vs the all-n a + b/n fit. Gate:
  < 0.1 mHa/atom. A miss means the fit form is wrong at these n. That is reported, not tuned away.
- **Reference:** Motta et al., PRX 7, 031059: −0.540493 Ha/atom. The result is stated as inside or
  outside the headline bar, with the gap in mHa/atom and in units of the bar. No other window or
  fit form is substituted to close a gap. Motta's value is itself an N → ∞ extrapolation, so a small
  disagreement is a finding about two extrapolations, not proof that either energy is wrong.

### 11.2 Result (2026-09-24)

**Run.** `benchmark_hchain_tdl.py --localize --protocol ramp --bond-dims 100,200,400 --threads 4
--stack-mem-gb 6` was run at n = 8, 10, 12, 16, 20, 24, 28, 32, 40. The tracked per-n table is
[`hchain_tdl_localized_table.csv`](hchain_tdl_localized_table.csv), which gate G4 in
`tests/test_hchain_largen2_spec.py` reads. Hardware: `Linux vm 6.18.44-fc-v37 x86_64`, Intel Xeon
@ 2.10 GHz, `nproc` = 4, 15 GiB RAM (container limit), no GPU.

| n | E_total (Ha) | E/atom (Ha) | extrap. stderr | dw at D_max | regime | D ladder | E_extrap − E(D_max) | vs FCI | wall (s) |
|---|---|---|---|---|---|---|---|---|---|
| 8  | −4.345079401 | −0.5431349 | 1.2e-14 | 3.7e-20 | converged | 100/200/400 | +2.9e-11 | 5.1e-11 | 19 |
| 10 | −5.424385375 | −0.5424385 | 2.2e-11 | 5.0e-20 | converged | 100/200/400 | −6.4e-11 | 1.3e-10 | 2 |
| 12 | −6.504226956 | −0.5420189 | 2.6e-10 | 5.9e-17 | converged | 100/200/400 | −4.2e-10 | 4.9e-10 | 6 |
| 16 | −8.664761408 | −0.5415476 | 4.1e-09 | 9.9e-15 | converged | 100/200/400 | −6.2e-09 | — | 18 |
| 20 | −10.825878902 | −0.5412939 | 2.2e-08 | 1.2e-13 | converged | 100/200/400 | −3.4e-08 | — | 41 |
| 24 | −12.987288247 | −0.5411370 | 1.8e-10 | 5.8e-13 | truncation | 100/200/400 | −2.2e-10 | — | 79 |
| 28 | −15.148861967 | −0.5410308 | 6.4e-10 | 1.8e-12 | truncation | 100/200/400 | −7.1e-10 | — | 142 |
| 32 | −17.310536221 | −0.5409543 | 9.6e-10 | 7.2e-12 | truncation | 100/200/400 | −1.2e-09 | — | 308 |
| 40 | −21.634061781 | −0.5408515 | 5.9e-08 | 1.5e-13 | converged | **200/400/800** | −9.1e-08 | — | 1317 |

Wall time is the DMRG time only (4 threads), and excludes FCI at n ≤ 12.

**n=40 needed the 200/400/800 rerun.** On the 100/200/400 ladder, the first ramp stage (D=100, 4
sweeps from a random MPS) stalled at −20.924 Ha, **0.71 Ha** above the D=200/400 energies. It
still reported a small, monotone discarded weight, so `regime` read "truncation". The
discarded-weight fit then extrapolated 0.82 mHa *below* its own D=400 energy (stderr 8.0e-4, 800 s).
**The regime label did not catch this.** Weights alone cannot see a stage that never converged.
The flag was the stderr, 4–5 orders above every other point, and E_extrap − E(D_max). The rerun's
D=400 energy (−21.6340616891) matches the stalled run's D=400 energy to every printed digit, and
its D=800 point moves it by 1e-9 Ha. G4 now also checks \|E_extrap − E(D_max)\| < 1e-6 Ha for every
vendored point. That check was added after this failure, and says so. (Since fixed at the source: `truncation_regime` now reads the energies too and rejects this ladder, reproduced -- see SPEC_regime_stalled_stage.md.) The stalled row is not in
the table. **n=48 and n=56** (200/400/800) were attempted and **OOM-killed** by the 15 GiB container
limit, with 13.9 GB resident after ~1 h each. Neither produced a row, so nothing half-finished is
in the fit.

**Pre-registered analysis** (`uv run python hchain_tdl_analysis.py`, n_max = 40):

| form | n_min | points | a (Ha/atom) | stderr |
|---|---|---|---|---|
| a + b/n | 8  | 9 | −0.540209 | 4.8e-05 |
| a + b/n | 12 | 7 | −0.540320 | 2.7e-05 |
| a + b/n | 16 | 6 | −0.540373 | 1.7e-05 |
| a + b/n | 20 | 5 | −0.540403 | 1.2e-05 |
| a + b/n + c/n² | 8  | 9 | −0.540492 | 1.6e-06 |
| a + b/n + c/n² | 12 | 7 | −0.540496 | 4.9e-07 |
| a + b/n + c/n² | 16 | 6 | −0.540495 | 5.8e-07 |
| a + b/n + c/n² | 20 | 5 | −0.540494 | 4.6e-07 |
| bulk (E(40) − E(32))/8 | — | 2 | −0.540441 | — |

- Envelope of a: [−0.540496, −0.540209], width 0.288 mHa. Largest stderr: 0.048 mHa.
- **Headline: e∞ = −0.540353 ± 0.000192 Ha/atom.**
- **LOO** (all-n a + b/n, dropping n=40): −0.540209 → −0.540172, shift **0.037 mHa/atom**. PASS (< 0.1).
- **Bulk vs fit:** \|−0.540441 − (−0.540209)\| = **0.232 mHa/atom**. FAIL (target < 0.1). This is
  reported, not tuned away. The a + b/n form is wrong at n ≤ 40. The linear fits drift monotonically
  as n_min rises (−0.540209 → −0.540403), which is the finite-size bias the old ± missed. The bulk
  quotient on (32, 40) still carries its own O(1/n²) curvature.
- **Motta −0.540493: INSIDE** the headline bar. The gap is +0.140 mHa/atom = 0.73 of the bar.

**Honest reading.** The bar is wide because the pre-registered envelope includes the biased
a + b/n fits. It contains Motta's value; the old ± (0.107 mHa, a fit stderr) did not, and missed
it by ~5×. A post-hoc observation, *not* the headline: the four a + b/n + c/n² fits agree with each
other to 4 µHa (−0.540492 to −0.540496), and with Motta to ≤ 3 µHa. Promoting that to the headline
would mean choosing the fit form after seeing the data, which §11.1 forbids. A future spec can
pre-register it. Motta's value is itself an N → ∞ extrapolation from larger N, so agreement or
disagreement at this level is a statement about two extrapolations. Scope: STO-6G, open chain,
R = 1.8 bohr. This is a model number that reproduces published physics; it does not extend it.
