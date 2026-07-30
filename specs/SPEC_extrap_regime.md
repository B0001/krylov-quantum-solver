# SPEC: the DMRG extrapolation regime label conflates "converged" with "uncontrolled"

**Status:** IMPLEMENTED once gates green. Backlog hypothesis: *"The `dweight` regime guard is
inverted — it flags the best DMRG points as out-of-regime"* (specs/BACKLOG.md, Open → Scale).

---

## 1. Goal

`hybrid_quantum_solver/dmrg_reference.py` labels every bond-dimension extrapolation `method =
"dweight" | "invD"`, and **four gate files assert `method == "dweight"` as a quality guard**. The
label is produced by

```python
usable = (len(per_D) >= 2 and not np.allclose(dws, 0.0)
          and np.all(np.diff(dws) <= 1e-12))
if not usable: x, method = 1.0 / Ds, "invD"
```

`np.allclose(dws, 0.0)` reduces to `max(dws) <= 1e-8` (numpy's default `atol`; the `rtol` term
vanishes against zero). So `"invD"` fires for **two opposite reasons** — the schedule *converged*
(nothing left to extrapolate) or the truncation is *uncontrolled* (non-monotone weights) — and the
four gates cannot tell them apart. **The falsifiable claim: any schedule converged enough that all
discarded weights fall below 1e-8 fails a gate written to demand quality.**

## 2. Background and honest framing

This repo has already recorded the phenomenon without naming it as a defect.
`SPEC_nbn_dmrg_reference.md:55` records the headline D=400/800/1200 run as
`invD (weights ~1e-9–1e-13: nothing left to extrapolate)`, and `:88-91` says the CI gate uses cheap
dims "**unlike the converged headline dims**" — i.e. the gate's bond dimensions were chosen
*downward* to keep the label. Same defect class as the recorded
`SPEC_validate_and_cost_composition` finding (one threshold quietly governing two boundaries).

**What we can claim if the gates pass:** the label is a two-into-three conflation; a three-state
regime separates them; and four gates were asserting a property their own specs show is not the
property they wanted.

**What we cannot claim — and this is the correction that matters.** The *axis fallback itself is
numerically justified and is NOT a bug.* On the recorded NbN converged weights the discarded-weight
axis has `cond(vander) = 2.1e9` against `1.4e3` for `1/D` — six orders worse conditioned, because
near-identical x-values make the slope meaningless. Switching axes there is correct. **Only the
label and the gates built on it are defective.** An earlier framing of this hypothesis implied the
fallback was itself wrong; that framing is falsified by the conditioning measurement above and is
withdrawn here. No energy, `stderr`, or extrapolated number changes in this spec.

**Not claimed:** that any recorded result is wrong. Every energy in the repo stands.

## 3. Approach

Extract the classification into a **pure predicate** over `per_D` triples so it is testable with
**no block2** and no DMRG run:

```
truncation_regime(per_D) -> "converged" | "truncation" | "uncontrolled"
```

Ordering is the one real design decision: the **floor test runs before the monotonicity test**, so a
converged schedule whose weights are non-monotone by float noise (e.g. `[1e-12, 3e-12, 1e-12]`)
classifies as `converged`, not `uncontrolled`.

`method` keeps its exact current value for every possible input — the new `regime` field sits
*alongside* it. That makes the migration provably numerics-free, and is verified by a bit-for-bit
oracle gate (G3) rather than by review.

**Reference:** the legacy expression itself (G3 oracle) plus the recorded weight profiles from
`SPEC_nbn_dmrg_reference.md:55` (converged) and `SPEC_hchain_largen.md` (the killed non-monotone
run).

**The floor constant.** `DISCARD_WEIGHT_FLOOR = 1e-8` is kept, now *named and justified* rather than
inherited from a numpy default. Near convergence `dE ≈ C·δ`; the recorded NbN ladder gives `C ≈ 30`
(dw 1e-9 → 1e-13 moves E by 3e-8 Ha), so a 1e-8 floor implies `|dE| ~ 3e-7 Ha ≈ 0.0003 mHa` —
three-plus orders below the tightest DMRG gate in the repo (`SPEC_singleramp` 0.1 mHa). **Honest:
the value is defensible but it arrived by accident.** That an implicit library default landed in the
right range is why nobody noticed; `C ≈ 30` is one system, so the constant is order-of-magnitude
justified, not derived.

## 4. Public interface

```
hybrid_quantum_solver.dmrg_reference.DISCARD_WEIGHT_FLOOR : float = 1e-8
hybrid_quantum_solver.dmrg_reference.REGIMES : tuple[str, ...]
hybrid_quantum_solver.dmrg_reference.truncation_regime(per_D, *, floor=...) -> str
ExtrapResult.regime : str        # NEW field; ExtrapResult.method unchanged
```

## 5. Acceptance criteria (validation gates)

`tests/test_extrap_regime_spec.py` — **pure, no block2, no pyscf**, so it runs even where the four
DMRG gates skip.

- **G1 — the three regimes separate (DEFINITION OF DONE).** The recorded NbN converged profile
  `[1e-9, 1e-11, 1e-13]` → `"converged"` (today: demoted, indistinguishable from failure); a healthy
  ladder → `"truncation"`; a non-monotone ladder above the floor → `"uncontrolled"`; `len < 2` →
  `"uncontrolled"`.
- **G2 — the floor boundary is pinned explicitly.** `max(dws) = 2e-8` → `"truncation"`;
  `= 2e-9` → `"converged"`. The physical threshold is asserted, never implied by a library default.
- **G3 — bit-for-bit `method` compatibility (the no-numerics-changed oracle).** For every fixture,
  the new `method` equals the legacy expression recomputed inline. One mismatch kills the claim that
  this change is label-only.
- **G4 — ordering: converged beats non-monotone.** `[1e-12, 3e-12, 1e-12]` (below floor *and*
  non-monotone) → `"converged"`. Kills a predicate that tests monotonicity first.
- **G5 — no gate compares against a string outside `REGIMES`.** Scan the four DMRG gate files for
  literals compared against `.regime`; every one must be in the imported `REGIMES`. This is the
  guard against shipping a green, dead assertion in a test that skips locally.

## 6. Implementation plan (test-first)

1. `tests/test_extrap_regime_spec.py` encoding G1–G5 (initially failing — no predicate).
2. `dmrg_reference.py`: add the constant, `REGIMES`, `truncation_regime`, the `regime` field; replace
   lines 193-198 with a call. No change to the fit or the returned energy.
3. Migrate the four assertions from `method == "dweight"` to `regime != "uncontrolled"`.
4. `make gates`.

## 7. Out of scope

- Changing any extrapolated energy, `stderr`, or fit. (Explicitly guarded by G3.)
- Retiring `method`. It answers a real question — which axis was fitted — and is recorded in spec
  tables and the `extrap_method` CSV column; rewriting recorded findings to fix a predicate is the
  wrong trade.
- A per-system or relative floor. `C ≈ 30` is one system; a defensible relative criterion needs the
  dE/dw slope measured across several systems — a follow-up.
- The pre-existing gate-hygiene bug found en route (below). Recorded, not fixed here.

## 8. Caveats and risks

- **R1 — a migrated assertion could be silently vacuous in a test that does not run.**
  `regime != "invD"` (a plausible slip: that value moved to `method`) passes forever, as does a
  typo. Mitigated by G5 plus importing `REGIMES` so a stale name is an `ImportError`, not a skip.
- **R2 — `regime != "uncontrolled"` admits one world the old assertion rejected** (the converged
  one), so it is strictly weaker as a set of passing worlds. That is the point — the admitted world
  is the better one — but the two gates with no independent accuracy assertion pair it with a
  positive check that `converged` is corroborated by the weights, not accepted on the label alone.
- **R3 — `max(dws)` not `dws[-1]`.** Using the last weight would label a normally-converging ladder
  (1e-3, 1e-5, 1e-9) `converged`, abandon the δ-fit and *silently change energies*. G3 catches it.
- **Pre-existing, found en route, not caused by this spec:** `tests/test_nbn_dmrg_reference_spec.py`
  lacks the `skipif(not dmrg_available())` guard its three siblings have, so without block2 it
  **fails** rather than skips. Recorded for a follow-up.
- **Environment, found en route:** `uv sync` failed for *every* extra on this machine because the
  `gpu` extra's `cupy-cuda12x` has no macOS wheels while `[tool.uv] required-environments` pins an
  x86_64-darwin resolution. Fixed with a platform marker so `--extra dmrg` can install at all.

## 9. Deliverables

- `hybrid_quantum_solver/dmrg_reference.py` — `DISCARD_WEIGHT_FLOOR`, `REGIMES`,
  `truncation_regime`, `ExtrapResult.regime`.
- `tests/test_extrap_regime_spec.py` — G1–G5, pure, runs without block2.
- Four migrated gate assertions.
- `pyproject.toml` — the cupy platform marker that unblocks `uv sync --extra dmrg`.
