# SPEC: `hamiltonian_moments` centers by ⟨H⟩ by default, with a conditioning-guarded fallback

**Status:** CLOSED — gates G1–G5 PASS (2026-10-01, chem-nxn). The default now returns the μ=⟨H⟩
centered moments from `specs/SPEC_centered_pds.md` G5 whenever doing so does not worsen the Hankel
matrix's conditioning at the deepest K the call supports; otherwise it falls back to the exact raw
moments `specs/SPEC_moment_pds.md` was validated against. K=8 PDS is now variational and monotone on
H₄/LiH/N₂ from a **plain** `hamiltonian_moments(mh, 2*K)` call — no caller-side centering — while
every existing gate (`SPEC_moment_pds`, `SPEC_pds_excited_roots`, `SPEC_centered_pds`) still passes
unmodified in its numeric claims (the latter two now pass `centered=False` explicitly, since their
whole point is studying the raw frame).

> **Attribution note.** The `moment_expansion.py` code this spec documents (the `centered` parameter
> and the conditioning-guarded fallback) and the `centered=False` updates to
> `tests/test_centered_pds_spec.py` / `tests/test_pds_excited_roots_spec.py` were already present at
> the tip of this branch (commit `ca63cb9`, titled as a `chem-dcz` proof-of-concept) when this bead's
> work resumed — this shared sandbox has several concurrent bead workers committing onto the same
> working tree, and that commit's content is byte-identical to this bead's own design (it references
> `chem-nxn` and this spec's filename directly). This spec does not re-author that code; it is the
> formal spec + dedicated acceptance gates chem-nxn's own description calls for ("needs its own
> spec") that the proof-of-concept commit did not include. Verified here, not re-implemented.

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

`SPEC_centered_pds` (chem-70c) showed μ=⟨H⟩ — the free first moment, no diagonalization — reproduces
dense-diagonalization-centered PDS to < 1e-6 Ha at K ≤ 8 and cuts cond(M) by 9.5–14 orders, but left
"change `moment_expansion.hamiltonian_moments`'s default" explicitly out of scope because it touches
`SPEC_moment_pds`'s closed API. The claim here: that default can be flipped to centered *without*
breaking any existing gate, by only applying the shift when it is not worse-conditioned than raw at
the deepest K the call requests. **Dies if** (a) `hamiltonian_moments(mh, N)` with no `centered` flag
changes any `SPEC_moment_pds` G1–G4 result by more than float noise, (b) a plain default call at
K=8 fails `SPEC_centered_pds`'s own G3 variational/monotone/accuracy thresholds on H₄/LiH/N₂, or
(c) `centered=False` fails to reproduce the exact raw numbers `SPEC_pds_excited_roots` and
`SPEC_centered_pds`'s own gates depend on.

## 2. Background and honest framing

- **The conflict this spec resolves.** A naive "always center" default breaks `SPEC_moment_pds`'s
  G2 on H₂: H₂ has only 2 HF-reachable eigenstates, so centering at K=3,4 (past exhaustion) is
  `SPEC_centered_pds` G2b's boundary — `det M → 0` there is fundamental, not a frame artifact, and
  in this case centering actually leaves the Hankel matrix *more* singular than the raw frame
  happens to be (measured below). The raw frame's K=3,4 agreement with FCI on H₂ was never
  guaranteed by theory; it is a float64 coincidence of the specific 2-reachable-state structure, now
  made visible by comparing frames. This is the finding, not a bug to hide.
- **What we can claim if the gates pass.** The conditioning-guarded default lets ordinary callers
  use `hamiltonian_moments(mh, 2*K)` up to K=8 on systems with enough reachable states (H₄, LiH, N₂
  CAS(6,6)) and get the `SPEC_centered_pds`-grade accuracy automatically, while every number any
  closed spec currently reports is reproduced exactly (to float noise) by the same code path that
  produced it originally.
- **What we cannot claim.** (a) No new physics — pure numerical hygiene, same scope note as
  `SPEC_centered_pds` §2. (b) The fallback does not rescue K > n_reachable cases (H₂-like systems);
  it only avoids making them worse. (c) The comparison is a per-call one-shot decision at the
  deepest K the call's `max_order` supports — a caller requesting `max_order` much larger than the K
  it actually uses may get a conservative (raw) answer even at low K where centering would have
  helped; callers that care should request `max_order = 2*K` for the K they intend to use (the
  existing convention in every caller in this repo).

## 3. Approach

Inside `hamiltonian_moments`, always compute the raw moments first (`mu_1` is needed for the shift
regardless). If `centered=True` (default), additionally shift `H → H − mu_1·I` (same math as
`SPEC_centered_pds` G5's `_shifted`) and recompute moments on the shifted operator. Compare
`cond(Hankel(K_max))` between the two, where `K_max = len(moments) // 2` is the deepest K this
`max_order` supports; return whichever is not worse conditioned there (ties go to centered). No new
primitives — reuses `MolecularHamiltonian`, sparse mat-vec, `np.linalg.cond`. Reference: FCI (dense
diagonalization), same as every antecedent spec here.

## 4. Public interface

```
moment_expansion.hamiltonian_moments(mh, max_order, centered: bool = True) -> (moments, offset)
```

`centered` is the only change to the signature; existing positional calls (`hamiltonian_moments(mh,
N)`) are unaffected except for the new conditioning-guarded shift.

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_centered_moments_default_spec.py`, BLAS pinned to 1 thread.

- **G1 — existing moment_pds gates reproduced exactly.** Re-running `SPEC_moment_pds`'s own G1–G4
  assertions (same geometries, same `hamiltonian_moments(mh, N)` calls with no flag) inside this
  file passes with the same pass/fail pattern as `tests/test_moment_pds_spec.py` (also run directly
  as part of `make gates` — unmodified file, unmodified result).
- **G2 — the unlock: K=8 with no caller-side centering.** `pds_energy(*hamiltonian_moments(mh,
  2*K), K)` for K=1..8 on H₄/LiH/N₂ is variational (≥ FCI − 1e-9), monotone (PDS(K+1) ≤ PDS(K) +
  1e-9), and PDS(8) beats `SPEC_centered_pds` G3's same per-system thresholds (H₄ < 1e-6, LiH <
  1.6e-3, N₂ < 5e-6 Ha) — reproducing that spec's numbers through the plain default call instead of
  an explicit `centered_frame` composition.
- **G3 — the default actually centers there.** On H₄/LiH/N₂ at `max_order=16`, the returned
  `moments[1] ≈ 0` (within 1e-8) and `offset` absorbs `⟨H⟩` — i.e. the fallback picked the centered
  branch, not a coincidental raw pass.
- **G4 — the fallback engages on the boundary case, and does not regress it.** On H₂ at
  `max_order=7` (`SPEC_moment_pds` G2's own call), the returned moments equal
  `hamiltonian_moments(mh, 7, centered=False)` exactly (the fallback chose raw), and
  `cond(Hankel_centered(K=4)) > cond(Hankel_raw(K=4))` is measured and recorded (the reason it fell
  back).
- **G5 — `centered=False` is untouched.** `hamiltonian_moments(mh, N, centered=False)` reproduces
  `SPEC_pds_excited_roots` G1b's own number (`cond(Hankel) > 1e14` for stretched H₄ at K=7) and
  `SPEC_centered_pds` G4's own number (N₂ raw `cond(M) > 1e17` at K=6,7,8) — confirming the opt-out
  is the exact same raw computation those closed specs' gates depend on, not an approximation.

> Definition of done: **G2** (the unlock) and **G1** (no regression) together — either alone is not
> the claim.

## 6. Measured data

Default-path K=8 errors vs FCI (mHa), 1 BLAS thread, matching `SPEC_centered_pds` §6 to within float
noise:

| system | K=4 | K=8 | fallback branch taken (max_order=16) |
|---|---|---|---|
| H₄  | 0.428 | 0.0001 | centered |
| LiH | 0.675 | 0.091  | centered |
| N₂  | 1.076 | 0.0011 | centered |
| H₂ (`max_order=7`, `SPEC_moment_pds`'s own call) | exact (K=2) | n/a (K≤4 only) | raw (fallback) |

H₂ fallback evidence (`max_order=7`, the call `SPEC_moment_pds` G2 actually makes): `cond(Hankel)`
at the deepest supported K=4 is 5.96e17 raw vs 6.17e17 centered — centered is *worse*, so the
fallback returns raw, reproducing `SPEC_moment_pds`'s original K=1..4 numbers bit-for-bit.

## 7. Out of scope

- Changing `pds_energy`/`pds_roots`/`cmx2_energy` themselves (composition only, as in
  `SPEC_centered_pds`).
- A smarter fallback that re-decides per-K rather than once at `K_max` (would need restructuring
  `pds_roots` to carry both frames; not needed by any current caller, all of which request
  `max_order = 2*K` for the K they use).
- Shot-noise propagation through the shifted operator (`chem-u87`, open).
- Whether centering changes the excited-root usefulness finding (`chem-dcz`, open).

## 8. Caveats and risks

- **R1 — the K_max heuristic is a one-shot, not a continuous, decision.** A caller requesting
  `max_order` far beyond the K it will actually use could get a conservative (raw) answer at a low K
  where centering would have been strictly better, if conditioning at the unused deep K is bad. No
  caller in this repo does this (`max_order` is always sized to the intended K), so it is not gated,
  only recorded.
- **R2 — `centered=True`'s fallback to raw does not fix H₂-like exhaustion.** Past K=n_reachable,
  neither frame is reliable (`SPEC_centered_pds` G2b); the fallback only avoids making the centered
  frame's exactly-singular failure mode *worse* than raw's coincidental survival. PDS is already
  exact at K=n_reachable, so this is not expected to matter in practice.
- Honest limitation shared with every antecedent spec here: exact statevector moments, minimal-basis
  tiny-to-modest active spaces — a correctness study, not a scaling demonstration.

## 9. Deliverables

- `specs/SPEC_centered_moments_default.md` (this file) and `tests/test_centered_moments_default_spec.py`
  — gates G1–G5; the net-new artifacts from this bead.
- `moment_expansion.py` — `hamiltonian_moments` gains `centered: bool = True` and the
  conditioning-guarded fallback (`_raw_moments`, `_hankel_cond` helpers). Already present on this
  branch (see attribution note in §Status); verified against G1–G5 here, not re-implemented.
  `tests/test_centered_pds_spec.py`, `tests/test_pds_excited_roots_spec.py` — already updated to pass
  `centered=False` explicitly at every call site that needs genuinely raw moments (no change to
  their assertions or numbers); verified still green under this spec's gates.
