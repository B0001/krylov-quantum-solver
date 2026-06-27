# SPEC: Hₙ thermodynamic limit at larger n (ramp protocol)

**Status:** CLOSED — gate **G1 FAILED** (see §5); recorded as a finding, not faked into a pass.
A follow-up (larger D, or a bulk per-site estimator) is needed. Builds on `SPEC_hchain_tdl.md` +
`SPEC_singleramp.md`.

---

## 1. Goal

Extend the Hₙ per-atom thermodynamic-limit estimate to larger chains (n up to ~40) using the
now-validated `protocol="ramp"` (≈2× cheaper), and check that the larger range **tightens** the
n→∞ fit. Claim: with more, longer chains the extrapolation is more stable.

## 2. Honest framing

Still the minimal-basis (sto-6g) open-chain MODEL at fixed R; reproduces benchmark physics, not a
new result. No new code physics — this exercises the validated extrapolation + TDL machinery at a
larger n range via the cheaper protocol.

## 3. Approach

`benchmark_hchain_tdl.py --protocol ramp --ns 8,10,12,16,20,30,40 --bond-dims 100,200,400`
(resumable; one uniform protocol across all n for a consistent fit).

## 4. Interface

Driver flags only: `--protocol`, `--ns`, `--bond-dims` (no library changes). The extrapolation and
TDL fit are unchanged and already gated (`test_hchain_tdl_spec`, `test_singleramp_spec`).

## 5. Acceptance criteria (result-level gate)

- **G1 — tighter, stable TDL fit.** Over the extended n range, the **leave-one-out shift of e_∞ is
  < 0.1 mHa/atom** (vs the 0.20 mHa/atom from the n≤16 study) — i.e. dropping the largest n barely
  moves the estimate. This is the "done" gate, checked against the study's CSV output.
- (The per-point soundness — ramp ≈ perD ≈ FCI — is already gated elsewhere; not re-tested here.)

### RESULT: G1 FAILED — and the failure is the finding.

Ran `--protocol ramp --ns 8,10,12,16,20,30,40 --bond-dims 100,200,400`. Per-point quality
**degrades sharply with n** at this (cheap) bond dimension:

| n | E/atom | stderr | method |
|---|---|---|---|
| 8–12 | −0.5431 … −0.5420 | 1e-16 … 3e-5 | dweight (clean, monotone) |
| 16 | −0.541678 | 1e-3 | dweight |
| 20 | −0.541899 | 5e-3 | dweight (E/atom non-monotone) |
| 30 | −0.539852 | 4e-3 | **invD** (discarded weights no longer usable) |

**Leave-one-out shift = 1.07 mHa/atom ≫ 0.1** → gate fails.

**Why:** a 1D chain's entanglement grows with length, so D=400 truncates too hard at large n —
the discarded weight stops being small/monotone (forcing the cruder `invD` fallback at n=30) and
the per-atom energy becomes unreliable (non-monotone at n=20, badly off at n=30). The **ramp
protocol is not at fault** (gated against per-D/FCI in `test_singleramp_spec`); the **bond
dimension is**. The hoped-for "cheap large-n" does not exist: reliable large-n Hₙ needs larger D
(≈ 400/800/1600), costing back most of the ramp speedup — consistent with the original per-D study
stalling at n=20.

**Revision (follow-up spec):** the gate is only achievable with adequate D at each n. Two honest
paths: (a) larger bond dims at large n (real compute), or (b) a *bulk* per-site estimator
`(E(n) − E(n−Δ))/Δ` from long chains, which cancels the open-boundary surface term and is far less
sensitive to total-energy noise. The committed n≤16 per-D result (e_∞ = −0.539967 ± 0.000107,
leave-one-out 0.20 mHa/atom) remains the best validated number.

## 6. Out of scope

- n up to 60 (kept as a stretch; gated only if n=40 completes in tractable time on this hardware).
- Complete-basis-set / periodic / equation-of-state extensions (separate specs).

## 7. Caveats and risks

- **R1:** even with the ramp protocol, large n at D=400 is the runtime pole; the driver is
  resumable, so a capped run still yields a valid fit from the completed n.
- Mixing protocols in one CSV is avoided by regenerating all n with `ramp`.

## 8. Deliverables

- `benchmark_hchain_tdl.py` — `--protocol/--ns/--bond-dims` flags.
- `data/hchain_tdl.csv` — uniform ramp study; e_∞ ± stderr + leave-one-out in the summary.
- `BACKLOG.md` — mark item done.
