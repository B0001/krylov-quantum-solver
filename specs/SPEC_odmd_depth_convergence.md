# SPEC: A depth-convergence check on ODMD's eps(K) — and the GapBracket premise flag it pairs with

**Status:** IMPLEMENTED — gates G1–G3 green (chem-36m, follow-up to chem-z3h /
`SPEC_odmd_temple_repair_map.md` §7, which named "a depth-convergence check on ε(K)" as explicit
out-of-scope future work).

---

## 1. Goal

Two gaps left by `SPEC_odmd_temple_repair_map.md`:

1. `temple_bounds.EnergyBracket` has an oracle-free `premise_refuted` flag (ε > θ₁ ⟹ the premise
   ε ≤ E₁ is proven false); `certified_gaps.GapBracket` rests on the identical premise (ε₁ ≤ E₁,
   checked against the same θ₁ by Cauchy interlacing) with no equivalent flag.
2. The measured map in `SPEC_odmd_temple_repair_map.md` §10 shows ODMD's ε(K) is **not monotone in
   K**: on H₄ the premise is repaired at K=16–26 and violated again at K=28 (+0.74 mHa); on N₂ it
   is repaired at K=12–16 and violated again at K=18–28 (+14.8…+35.4 mHa). Treating "the premise
   currently holds at this K" as a green light to substitute ε(K) is unsound — the next K can flip
   it back. A depth-convergence check (ε(K) stable over a trailing window of already-measured K,
   no look-ahead) is the cheap, oracle-free guard that spec named and deferred.

Falsifiable: (1) `GapBracket.premise_refuted` kills the claim if it ever fires without a genuine
ε₁ > E₁ overshoot, or fires in self mode. (2) the depth-convergence check kills the claim if any
K it declares "converged" coincides with a measured containment ESCAPE anywhere in the gated
(M, K) map — i.e. if requiring convergence does not buy any actual safety over using ε(K) raw.

## 2. Background and honest framing

- `SPEC_odmd_temple_repair_map.md` already built and gated the (M, K) map this spec reuses
  (`odmd_temple_repair_map.k_vs_validity_map`, `K_DEPTHS`, `TOL`); no new Krylov/ODMD machinery.
- **What we can claim if the gates pass:** `GapBracket.premise_refuted` is sound (never fires
  without a real overshoot, never in self mode) by the same Cauchy-interlacing argument as
  `EnergyBracket`'s; and requiring depth-convergence (stability over a trailing window of K)
  before trusting ε(K) as a Temple floor never green-lights a K that the gated map shows produces
  a containment escape, on any of the four non-trivial systems (H₄, N₂, LiH CAS(2,5), full-space
  LiH) at any base depth M.
- **What we cannot claim:** depth-convergence is **necessary, not sufficient**. ε(K) can stabilize
  at a value that still violates the premise — N₂'s ε(K) is stable within 1 mHa over K=24–28 while
  sitting +34…+35 mHa *above* E₁ the whole time (a consistently-biased fixed point, not a
  transient). The check catches the *documented* failure mode (an unstable repair that was about
  to flip back); it is blind to a stable-but-wrong one, exactly the model-misspecification blind
  spot `SPEC_odmd_uq` G4 already recorded for single-signal resampling. Pairs with, does not
  replace, the oracle-free `premise_refuted` flag and an eventual Lehmann/Kato rigorous floor
  (still out of scope, as it was before).
- Thresholds (window=3 grid points, tol=2 mHa) were chosen by inspecting the already-gated
  `SPEC_odmd_temple_repair_map.md` §10 numbers (committed before this spec; re-derived here, not
  re-measured blind) — the same provenance discipline that spec used for its own pre-registration.

## 3. Approach

- **Reference:** the already-gated `odmd_temple_repair_map.k_vs_validity_map` cells (exact
  reachable E₀/E₁ by dense diagonalization feeds `summarize`'s `escapes` set — an oracle, used only
  to check the new flag/check, never to compute them).
- `certified_gaps.GapBracket.premise_refuted` — `eps1 > theta1 + 1e-9` (identical form and margin
  to `temple_bounds.EnergyBracket.premise_refuted`; `GapBracket.theta1` gains the same shot-noise
  guard `EnergyBracket.theta1` already has, so the new flag cannot look sound and not be — no
  current caller passes a noisy solver to `gap_bracket`, so this is a no-op for every existing use).
- `odmd_depth_convergence.depth_converged(odmd, k_grid, k, window=3, tol=2e-3)` — True iff the
  `window` deepest grid points at or below `k` are all resolved (not `None`) and span `<= tol` Ha.
  No look-ahead: only K's at or below the query depth are used, matching how this would run live
  (increase K until the estimate stops moving, then stop).
- `odmd_depth_convergence.safe_eps1(mh, solver, k, ...)` — `odmd_temple_repair_map.odmd_eps1`
  gated by `depth_converged`; returns `None` (abstain) rather than an unconverged estimate. Feeds
  either `temple_bounds.krylov_bracket(eps=...)` or `certified_gaps.gap_bracket(e1=...)` — both
  accept an externally supplied E₁ floor, so one guard protects both brackets' substitution path.

## 4. Public interface

```
certified_gaps.GapBracket.premise_refuted            # NEW property, mirrors EnergyBracket's
odmd_depth_convergence.depth_converged(odmd, k_grid, k, window=3, tol=2e-3) -> bool
odmd_depth_convergence.converged_depths(odmd, k_grid, window=3, tol=2e-3) -> list[int]
odmd_depth_convergence.safe_eps1(mh, solver, k, k_grid=K_DEPTHS, window=3, tol=2e-3) -> float|None
```

## 5. Acceptance criteria (validation gates)

`tests/test_odmd_depth_convergence_spec.py`. Systems: H₄, N₂ CAS(6,6), LiH CAS(2,5), full-space
LiH (`SPEC_odmd_temple_repair_map`'s four non-trivial systems; H₂ excluded, no escapes to check
against). K grid and TOL as in that spec.

- **G1 — `GapBracket.premise_refuted` is sound (DEFINITION OF DONE, part 1).** Never fires without
  `eps1 > E1 + TOL` against the oracle; never fires in self mode on any (system, M) in the grid; a
  deliberately-wrong oracle `e1` set above `theta1` does fire it; a noisy solver leaves `theta1` at
  `+inf` so the flag stays `False` regardless of the (now-unverifiable) `eps1`.
- **G2 — depth-convergence buys real safety (DEFINITION OF DONE, part 2).** Across H₄, N₂,
  LiH CAS(2,5), full-space LiH: no K that `depth_converged` declares True coincides with any K in
  that system's measured containment-escape set (any M). One counterexample kills the claim that
  the check is a genuine (if partial) safety improvement over substituting ε(K) raw.
- **G3 — necessary, not sufficient (the honest boundary).** At least one (system, K) exists where
  `depth_converged` is True yet the premise is still violated against the oracle (`eps > E1 + TOL`)
  — demonstrating convergence does not imply correctness. Measured on N₂ at K ∈ {26, 28}.

## 6. Implementation plan (test-first)

1. `tests/test_odmd_depth_convergence_spec.py` encoding G1–G3 (G1 red: no `premise_refuted` on
   `GapBracket`; G2/G3 red: no `odmd_depth_convergence` module).
2. `certified_gaps.py`: add `GapBracket.premise_refuted` + the `theta1` noise guard.
3. `odmd_depth_convergence.py`: `depth_converged`, `converged_depths`, `safe_eps1`, reusing
   `odmd_temple_repair_map`'s map/ladder.
4. `make gates-odmd_depth_convergence` and `make gates-certified_gaps` until green.

## 7. Out of scope

- Making the substitution unconditionally safe (Lehmann/Kato rigorous E₁ floors) — still deferred,
  as `SPEC_odmd_temple_repair_map.md` §7 already said.
- Centered-frame ODMD, shot noise on the ODMD signal itself, K or M off the grid.
- Wiring `safe_eps1` into a live pipeline default — this spec delivers the guard and measures its
  coverage; using it by default in `pipeline.py` is a follow-up product decision, not a gate here.
- Tightening `depth_converged`'s window/tol beyond what G2/G3 need — a smaller tol buys nothing
  once G2 already holds at window=3, tol=2 mHa; a parameter sweep is future work if the boundary
  in G3 needs to move.

## 8. Caveats and risks

- **R1 — G3 is a feature, not a bug, but it means `depth_converged=True` must never be read as
  "this ε(K) is correct."** Document this at every call site; `safe_eps1` returning a value means
  "stopped moving," not "verified."
- **R2 — brittle exact K-sets.** Like the parent spec, the specific converged/escaping K's are a
  measurement on this machine; §10-style provenance notes in the test docstrings, not re-derived
  thresholds, carry the "why these numbers" burden if a cell flips elsewhere.
- Sector-restricted, exact statevector, uncentered ODMD frame — identical scope to
  `SPEC_odmd_temple_repair_map.md`; this spec adds no new physics, only a guard and a flag.

## 9. Deliverables

- `certified_gaps.py` — `GapBracket.premise_refuted`, `theta1` noise guard.
- `odmd_depth_convergence.py` — `depth_converged`, `converged_depths`, `safe_eps1`.
- `tests/test_odmd_depth_convergence_spec.py` — gates G1–G3.
- `specs/SPEC_odmd_depth_convergence.md` — this spec.
- `specs/BACKLOG.md` — entry recording the result.
