# SPEC: Hₙ to larger n, done right (adequate-D ramp + bulk per-site estimator)

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

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
- **R2 — cost.** D=1600 ramp at n≈30 is ~an hour; the driver is resumable and the CI gates stay at
  small n. Large-n runs belong on a faster/GPU node.
- Honest limitation: minimal-basis open-chain model; reproduces, does not extend, known physics.

## 9. Deliverables

- `hybrid_quantum_solver/dmrg_reference.py` — `bulk_per_site_energy`.
- `benchmark_hchain_tdl.py` — bulk estimate + leave-one-out shift in the printed summary.
- `tests/test_hchain_largen2_spec.py` — gates G1–G3.
- Results summary (with §2/§7 caveats + the large-n LOO number) in the PR description.
