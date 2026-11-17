# chem-36m — `GapBracket.premise_refuted` + a depth-convergence check on ODMD's ε(K)

Status: **closed, both halves of the bead delivered**. Acceptance criteria (my own
pre-registration, in `specs/SPEC_odmd_depth_convergence.md`) met; all touched/adjacent gates green.

## What changed

- `certified_gaps.py` — `GapBracket` gains a `premise_refuted` property, byte-for-byte mirroring
  `temple_bounds.EnergyBracket.premise_refuted`: `eps1 > theta1 + 1e-9` (oracle-free proof the
  Weinstein/self premise `eps1 <= E_1` failed, via Cauchy interlacing `theta1 >= E_1`). While adding
  it I found `gap_bracket`'s `theta1` field was *not* guarded against shot noise the way
  `EnergyBracket.theta1` already is (`temple_bounds.krylov_bracket` sets `theta1=np.inf` when
  `solver.noise_sigma` is truthy, since under noise θ₁ is no longer a variational bound — so an
  unguarded `theta1` would make the new flag unsound under noise). Fixed by adding the identical
  guard to `gap_bracket`'s return. Confirmed via repo-wide grep that no current caller passes a
  noisy solver to `gap_bracket` (`certchem/core.py::_estimate_gap` always builds a fresh noiseless
  solver; `gap_selfcheck_noise.py` constructs `GapBracket` directly, bypassing `gap_bracket`
  entirely) — so this is a no-op for every existing call site, purely a latent-soundness fix.
- `odmd_depth_convergence.py` (new) — `depth_converged(odmd, k_grid, k, window=3, tol=2e-3)`:
  True iff the `window` deepest grid points at or below `k` are all resolved (not `None`) and span
  `<= tol` Ha (no look-ahead — only K's already measured at or below the query depth). 
  `converged_depths(odmd, k_grid, ...)`: all K in `k_grid` where that holds.
  `safe_eps1(mh, solver, k, ...)`: `odmd_temple_repair_map.odmd_eps1(mh, solver, k)` gated by
  `depth_converged` — returns `None` (abstain) rather than an unconverged estimate; the function
  other code should call instead of `odmd_eps1` directly before substituting into
  `temple_bounds.krylov_bracket(eps=...)` or `certified_gaps.gap_bracket(e1=...)` (both accept an
  externally supplied E₁ floor).
- `specs/SPEC_odmd_depth_convergence.md` (new) — full spec: goal, background/honest framing,
  approach, public interface, G1–G3 acceptance gates, implementation plan, out of scope, caveats.
- `tests/test_odmd_depth_convergence_spec.py` (new) — gates G1–G3 plus two unit tests for the
  convergence/abstention logic itself. Reuses `odmd_temple_repair_map`'s already-gated (M, K) map
  (chem-z3h) — no new Krylov/ODMD machinery, same oracle (exact reachable E₀/E₁ via dense
  diagonalization, used only to check the flag/check, never to compute them).
- `specs/SPEC_certified_gaps.md` — added a one-line interface entry + a short §10 addendum noting
  the new `premise_refuted` property and pointing at the new spec, for discoverability. No gates
  changed in that spec; G1–G4 there are untouched.
- `specs/BACKLOG.md` — new `[x] CLOSED (chem-36m, ...)` entry before the `## Killed` section,
  following the existing closed-entry convention (bold one-line claim, measured findings, spec/test
  links).

## Pre-registered acceptance criteria (`specs/SPEC_odmd_depth_convergence.md` §5) and pass/fail

**G1 — `GapBracket.premise_refuted` is sound (definition of done, part 1).** ✅ Never fires without
`eps1 > E1 + TOL` against the oracle, on any of H₄/N₂/LiH CAS(2,5)/full-space LiH at M∈{2,4,6,8};
never fires in self mode; a deliberately-wrong oracle `e1` set above `theta1` DOES fire it; under
shot noise `theta1` is left at `+inf` so the flag stays `False` even though the (same bad) `eps1` is
then unverifiable — soundness over usefulness, matching `EnergyBracket`'s existing contract exactly.

**G2 — depth-convergence buys real safety (definition of done, part 2).** ✅ Across H₄, N₂ CAS(6,6),
LiH CAS(2,5), full-space LiH: no K that `converged_depths` returns coincides with any K in that
system's measured containment-escape set (`lower > E0 + TOL`, chem-z3h's "escapes" field — the
*actual* bracket break, not mere premise violation) at any base depth M on the grid. E.g. full-space
LiH escapes at M=6, K∈{12,14,16}; `depth_converged` is False at those K because the trailing-3-point
window spans 174–698 mHa there, nowhere near the 2 mHa tolerance (ε(K) is still swinging wildly).

**G3 — necessary, not sufficient (the honest boundary).** ✅ At least one (system, K) has
`depth_converged=True` while the premise is still violated against the oracle. Measured: N₂'s
ε(28) sits +34.5 mHa above E₁ while the trailing window (K=24,26,28: 35.38, 34.97, 34.53 mHa
overshoot) spans only 0.84 mHa — well under the 2 mHa tolerance. A stable-but-wrong fixed point,
not a transient. (N₂ at K=26 does *not* converge under window=3 — its own trailing triple
(22,24,26) spans 2.26 mHa, just over tolerance — so the actual converged-and-violated set is `{28}`
on this grid, not `{26,28}` as I'd guessed from the parent spec's prose before running the numbers;
G3 only requires the set be non-empty, which it is.)

## Commands run and verbatim results

```
$ uv run pytest -q tests/test_odmd_depth_convergence_spec.py
.....                                                                    [100%]
5 passed, 10 warnings in 160.32s (0:02:40)

$ uv run pytest -q tests/test_certified_gaps_spec.py        # regression check, pre-existing gate
....                                                                     [100%]
4 passed, 3 warnings in 460.41s (0:07:40)

$ uv run pytest -q tests/test_odmd_temple_repair_map_spec.py  # regression check, pre-existing gate
....                                                                     [100%]
4 passed, 6 warnings in 321.69s (0:05:21)

$ uv run ruff check .
All checks passed!
```

All three pytest invocations run in their own process (block2-isolation convention; none of these
three import `block2`/DMRG, but this also matches how `make gates` would run each
`tests/test_*_spec.py` individually). Environment: `uv sync --extra dmrg --extra test --extra dev`
already in place from a prior session step (not re-run here since `uv.lock`/`pyproject.toml` were
untouched by this bead).

## Calibration provenance (window=3, tol=2 mHa)

These thresholds were chosen by hand-inspecting the real eps(K)−E₁ ladders for all four non-trivial
systems (computed via `odmd_temple_repair_map.k_vs_validity_map`) before writing
`odmd_depth_convergence.py`, then confirmed by the actual pytest run above (not re-measured blind
after the fact). Representative N₂ ladder (mHa, `E1 = -106.9151343690812` Ha):

```
K:    4       5       6       7       8      10       12       14       16      18      20      22      24      26      28
eps-E1: 1529.71 1481.28 1026.31 1027.13 544.04  342.04 -117.88 -332.32 -283.05  14.84   19.72   33.12   35.38   34.97   34.53
```

Trailing-3 span at K=28 (24,26,28): 35.38−34.53 = 0.84 mHa ≤ 2 mHa → converged, still +34.5 mHa
above E₁ → G3's witness. Full ladders for H₄, LiH, LiH CAS(2,5) are reproducible with
`uv run python odmd_depth_convergence.py` (prints `converged K` per system; not committed as a data
file — cheap to regenerate, ~2–3 min per system on this machine, PySCF/qiskit statevector path).

## What I decided not to do, and why

- **Did not wire `safe_eps1` into any live pipeline path** (e.g. `pipeline.py`,
  `odmd_temple_repair_map.py`'s own cells). The spec's §7 explicitly scopes this out: delivering the
  guard and measuring its coverage is this bead; making it the pipeline default is a separate
  product decision with its own tradeoffs (e.g. the extra ODMD sweep cost at every K up to the
  target depth) that wasn't asked for here.
- **Did not attempt a Lehmann/Kato rigorous E₁ floor.** Both the parent spec (chem-z3h, §7) and this
  one explicitly defer it — it needs a genuine lower bound on E₂, which real-time Krylov data alone
  doesn't give you; out of scope for this bead's ask.
- **Did not sweep `window`/`tol`** beyond what G2/G3 needed. Once window=3/tol=2mHa satisfies both
  gates on all four systems, tightening further buys nothing without a reason — recorded as explicit
  out-of-scope in the spec (§7) rather than silently skipped.
- **Did not re-run the full `make gates` (~90 spec files).** Ran the one new gate plus the two
  gates whose *production* code (`certified_gaps.py`) or *adjacent* module
  (`odmd_temple_repair_map.py`, unmodified but closely coupled) this bead touched. A full run would
  burn a large share of the shared compute envelope on ~90 gates unrelated to this change; `make
  lint`/`ruff check .` (repo-wide, cheap) did run clean.
- **Did not touch `odmd_temple_repair_map.py` itself.** It's read-only infrastructure this bead
  reuses (`odmd_eps1`, `K_DEPTHS`, `K_vs_validity_map`, `summarize`); its own gate (re-run above)
  confirms it's unaffected.

## What I could not verify

- **Whether `window=3, tol=2mHa` generalizes beyond this spec's five-system, fixed-geometry,
  fixed-`K_DEPTHS` grid.** Like the parent spec (chem-z3h, R2), the specific converged/escaping K
  sets are a measurement on particular systems and this machine's float path, not a derived
  universal constant — flagged explicitly in the spec (§8, R2) and the BACKLOG entry.
- **Shot-noise / hardware behavior of the convergence check itself.** Everything here is
  noiseless statevector (inherited scope from `odmd_temple_repair_map.md`); out of scope, listed in
  §7.
- **Did not independently re-verify the ad-hoc manual spot-check I ran earlier in this session**
  (constructing `GapBracket` with a bad oracle and a noisy solver via a throwaway script) against
  the final committed code — the actual `test_G1_...` pytest gate covers the same cases
  end-to-end and passed, which supersedes that earlier manual check.

## Scope note

No work outside this bead's scope was discovered that needed a new bead — the two asks (flag +
depth-convergence check) were self-contained once the parent spec's already-gated map was reused.
