# SPEC: Does an HF-weighted band-edge gap predict the Krylov depth the chem-pm3 η-bound needs? (chem-j71)

**Status:** DONE 2026-10-15 (pre-registered 2026-10-14; this section and §4–§5 were written, and the
new-witness search was kicked off, before any correlation number below was computed; the new-witness
*search* log was inspected to decide which geometries survive at all, which is a selection step, not
a peek at the correlation result it is used to test). **Result (§10): G1 PASS; G2 SUPPORTS the
weighted-gap hypothesis (ρ_WG=−0.867 vs ρ_overshoot=−0.758, n=10); G3 finds the moment-bound LP is
~100% unsound under measured shot noise at every tested shot count up to 1e8, at a scale matching
chem-pm3's own analytic worst case rather than beating it.**
**Depends on:** `SPEC_pm3_eta_bound.md` (chem-pm3 — the bound `_build_moment_bound_lp` and the three
original escape witnesses), `SPEC_subspace_floor_resolvability.md` (the escape-witness mechanism).

---

## 1. Goal

`SPEC_pm3_eta_bound.md` §10 found that the self-mode floor's overshoot `β − E_d` does not order its
three escape witnesses by how hard the η-bound is to satisfy: the asym witness overshoots by 5× R=1.1's
overshoot yet is the easier of the two (survives at M ≤ 24; R=1.1 needs M ≈ 60). This spec asks the
bead's two follow-up questions about that same bound:

- **(A)** does a band-edge gap that also accounts for the HF weight sitting at it predict the bound's
  M-dependence (required depth / achieved bound) better than the overshoot does, across a wider witness
  set than chem-pm3's three?
- **(B)** how much does the bound actually move under realistic shot noise on the moments `S_0k`, as
  opposed to the analytic worst case chem-pm3 §7 derived but did not measure
  (≈ (2M−1)·cap·max|ΔS_0k| ≈ 4.7e4·max|ΔS_0k| at M=24, cap 1e3)?

## 2. Background and honest framing

- **What's already established** (chem-pm3, not re-derived here): `_build_moment_bound_lp` is a sound
  upper bound on η for any measure on `[e_min, e_max]`, for any coefficient box, monotone in M and in
  the box. The three witnesses' verdicts (R=1.0 survives, R=1.1 killed, asym survives) are fixed facts
  this spec reuses, not reruns.
- **Why overshoot is the wrong variable (already shown).** Overshoot sets the *width* of the band
  `[E_d, β]` the bound must empty out, but says nothing about how much weight sits at the band's edges
  or how close those edges are to the band itself — and it is exactly that edge structure a degree-(M−1)
  trigonometric polynomial has to resolve.
- **The candidate predictor (pre-registered, one specific operationalization — not a claim that this is
  *the* right functional form from first principles).** Define, per witness:
  - `gap_lo = E_d − E_{d−1}` and `gap_hi = (first reachable level ≥ β) − β` (already computed by
    `run_witness`),
  - `η` = the exact HF weight on reachable levels in `[E_d, β)` (already computed by `run_witness`),
  - **weighted gap** `WG = min(gap_lo, gap_hi) / η`. Larger `WG` (more room relative to the weight that
    must be excluded) is predicted to mean an easier bound (survives at smaller M); smaller `WG` means
    harder (killed, or killed at a larger required M).
  - This is one simple, cheap, falsifiable candidate — not a derivation from the Markov–Krein extremal
    theory of the LP (that theory is not worked out here; see §7).
- **What you can claim if the gates pass:** that `WG` orders / correlates with the bound's difficulty
  better than the overshoot does, on this witness set, and by how much; or that neither predicts it
  well, which is itself the honest answer to "does a cheap scalar explain this."
- **What you cannot claim:** a general law relating gap/weight to required Krylov depth (this is single
  predictor on O(10) witnesses from one molecule family, not a derived asymptotic result); an
  oracle-free predictor (`E_d`, `β`, and the support are still oracle inputs, exactly as in chem-pm3).

## 3. Approach

**Part A — wider witness set.** Reuse `scripts/spec_pm3_subspace_eta_bound.py`'s `_chain`, `_asym_chain`,
`_ritz_floor`, and `run_witness` unchanged. New candidate geometries (linear H₆ at
R ∈ {0.85, 0.90, 0.95, 1.05, 1.15, 1.20, 1.25} Å, plus two asym variants
`[0.9,0.9,2.0,0.9,0.9]` and `[0.9,0.9,2.5,0.9,0.9]`, all d=3) are searched over
M ∈ {6, 8, 10, 12, 16, 20, 24} for the **smallest M at which the disjoint-interval guard passes AND the
self-mode floor overshoots the true E_d by > 1e-6 Ha** (the same "guard-passing, invalid-floor" escape
condition `SPEC_subspace_floor_resolvability.md` defines) — a deterministic, pre-committed selection
rule applied *before* any witness's moment-bound number exists, so it cannot be tuned post hoc to flatter
either predictor. A candidate that never escapes in that M range contributes no witness (recorded, not
discarded silently). `pyscf.lib.num_threads(1)` is pinned around each SCF/integral build (bd memory
`chem-7cx-root-cause`: free-RHF-branch nondeterminism on near-degenerate geometries is otherwise an
OpenMP thread-count artifact, not physics).

Each escaping candidate plus the three original witnesses gets the full `run_witness` treatment
(M_SWEEP = 8,12,16,20,24; extended sweep to find `m_pass` if killed). `WG` and the overshoot are computed
for all of them, and compared against `b_min` (and `m_pass` for killed witnesses) by Spearman rank
correlation.

**Part B — measured noise sensitivity.** On the R=1.1 witness (the one chem-pm3's analytic estimate was
computed for) at M=24, cap=1000: perturb `moments[1:24]` (the measured `S_0k`, k ≥ 1; `moments[0] = 1`
is held exact, mirroring `QuantumKrylovSolver._hermitian_shot_noise`'s `unit_diagonal=True` treatment of
a self-overlap) with iid complex Gaussian noise of std `noise.shot_noise_sigma(shots)` per real/imaginary
part, for `shots` ∈ {1e2, 1e3, 1e4, 1e5, 1e6, 1e8} and 40 seeded realizations per shot count. Recompute
the bound each time; report the empirical mean/std/[min,max] of the bound and the fraction of
realizations where the (no-longer-exact-moment) bound falls below the true η — the soundness proof
assumed exact moments, so this is the honest thing to check, not an assumption to carry over. Compare
the empirical std against the chem-pm3 worst-case slope `(2M−1)·cap·σ` and against a typical-case
estimate `‖x*‖_2 · σ · sqrt(2(M−1))` (`x*` = the optimal LP coefficient vector at the exact moments).

## 4. Public interface

```
scripts/spec_pm3_depth_predictor.py
  find_escape_witnesses(candidates, d, m_search) -> list[dict]   # Part A geometry search
  weighted_gap(result) -> float                                  # WG = min(gap_lo, gap_hi) / eta
  rank_correlation_report(results) -> dict                       # Spearman(WG, b_min), Spearman(overshoot, b_min)
  noise_sensitivity_sweep(witness, M, cap, shots_list, n_trials) -> dict
  main()
```

## 5. Acceptance criteria — pre-registered

- **G1 — new witnesses reuse the validated primitives.** Every new witness's G-a and G-b-soundness
  checks pass exactly as `test_pm3_eta_bound_spec.py`'s do (bound ≥ η, bound non-increasing in M, lower
  ≤ bound) — this is a regression check that the new geometries are not somehow breaking the LP, not a
  new claim.
- **G2 — the predictor comparison is recorded either way.** Spearman ρ(WG, b_min) and
  ρ(overshoot, b_min) computed over all witnesses (original 3 + new escapes) and reported. **Supports
  the hypothesis** if `|ρ(WG, b_min)| ≥ 0.7` and exceeds `|ρ(overshoot, b_min)|`; anything else
  (including a clean kill) is recorded as the finding, not reworked until it passes.
- **G3 — noise sensitivity is measured, not just asserted.** For every `shots` value: report mean, std,
  [min, max] of the bound over 40 trials, and the fraction of trials where the bound is unsound
  (`< η_exact`). Report both analytic comparators (worst-case slope, typical-case `‖x*‖_2` estimate)
  alongside the empirical std at each shot count.
- **Descriptive, not gated:** whether the fraction of unsound trials is ever non-negligible at
  shot counts realistic for near-term hardware (1e4–1e6) — if the LP bound is frequently unsound under
  noise, that is a load-bearing caveat for ever using it outside the noiseless setting chem-pm3 scoped
  it to, and must be stated plainly in §10, not softened.

Done = G1–G3 green (G1 a true pass/fail; G2/G3 "green" means *computed and recorded*, since both can
honestly resolve either way) and §10 filled in.

## 6. Implementation plan

1. `tests/test_pm3_depth_predictor_spec.py`: G1 regression gates on the new witnesses; G2/G3 as
   "ran and recorded" checks (structure/shape/monotonicity assertions that must hold regardless of which
   way the substantive finding goes — e.g., the noise sweep's mean bound must still be ≥ a 0-noise lower
   sanity floor, Spearman must be computed on the actually-recorded witness count, not a hardcoded one).
2. `scripts/spec_pm3_depth_predictor.py`: the search, the predictor, the correlation report, the noise
   sweep, importing `run_witness`/`_build_moment_bound_lp`/etc. from
   `scripts/spec_pm3_subspace_eta_bound.py` unchanged.
3. Run, record §10, correct `specs/BACKLOG.md`'s chem-j71 entry.

## 7. Out of scope

- Deriving `WG`'s functional form from the LP's actual extremal (Markov–Krein / Fejér-kernel) theory —
  §2 is explicit that this is one cheap candidate, not a derivation.
- A corrected, sound bound under noise (Part B only measures how unsound/how volatile the existing,
  noiseless-only bound is; it does not fix it).
- Oracle-free band edges — unchanged from chem-pm3.
- Any witness family other than linear/asym H₆ chains (same scope chem-pm3 used).

## 8. Caveats and risks

- **Small-n statistics.** Likely single-digit to low-teens witness count; a Spearman ρ here is
  suggestive, not a law. Report the exact n alongside ρ.
- **Selection-rule dependence.** "Smallest escaping M" is one deterministic rule; a different rule
  (e.g., largest escaping M in range) could change which witnesses exist and shift ρ. Stated once,
  applied uniformly, not re-chosen after seeing results.
- **Noise model fidelity.** The per-`S_0k` iid complex Gaussian with std `shot_noise_sigma(shots)`
  mirrors `noise.py`'s documented scale but is not bit-identical to
  `QuantumKrylovSolver._hermitian_shot_noise`'s Hermitian-averaged perturbation (that halves the
  variance via symmetrization); stated as a modeling choice, not hidden.

## 9. Deliverables

- `scripts/spec_pm3_depth_predictor.py` — search + predictor + noise sweep driver.
- `tests/test_pm3_depth_predictor_spec.py` — gates.
- §10 results; the chem-j71 entry in `specs/BACKLOG.md` corrected to point here.

## 10. Results

**Run 2026-10-15, Linux x86_64 (8 vCPU / 16 GB RAM container).** Commands and raw logs:

```
uv run python scripts/spec_pm3_depth_predictor.py   # full numbers below (/tmp/j71_main_run.log)
uv run pytest tests/test_pm3_depth_predictor_spec.py -v   # 3 passed, 19 warnings in 1218.16s
```
Regression check that the shared, unmodified `scripts/spec_pm3_subspace_eta_bound.py` primitives
still behave as chem-pm3 recorded them:
```
uv run pytest tests/test_pm3_eta_bound_spec.py tests/test_subspace_floor_resolvability_spec.py -q
# 47 passed, 8 warnings in 1267.18s (ran concurrently with the driver above, hence the wall time)
```

### Part A — witness search and the predictor comparison

`find_escape_witnesses()` over the 9 pre-registered candidates found **7 escapes** (2 candidates,
R=0.85 and R=1.25, never escaped in M ∈ {6,8,10,12,16,20,24} and contribute no witness — recorded,
not discarded silently). Combined with chem-pm3's 3 original witnesses, n = 10 witnesses got the
full `run_witness` treatment (M_SWEEP = 8,12,16,20,24; `m_extended=()`, `caps=()` — see deviation
note below):

| Witness | WG = min(gap_lo,gap_hi)/η | overshoot (Ha) | b_min (M=24) | survives (≤0.05)? |
|---|---:|---:|---:|---|
| linear H₆ R=1.0 *(chem-pm3)* | 2.4247 | +0.2821 | 0.0396 | True |
| linear H₆ R=1.1 *(chem-pm3)* | 0.3624 | +0.2155 | 0.0594 | **False** |
| asym H₆ [0.9,0.9,2.2,0.9,0.9] *(chem-pm3)* | 0.6624 | +1.0548 | 0.0359 | True |
| linear H₆ R=0.90 | 8.3542 | +0.3070 | 0.0217 | True |
| linear H₆ R=0.95 | 4.8622 | +0.2927 | 0.0291 | True |
| linear H₆ R=1.05 | 0.0522 | +0.2765 | 0.0522 | **False** |
| linear H₆ R=1.15 | 0.2575 | +0.2269 | 0.0747 | **False** |
| linear H₆ R=1.20 | 0.3069 | +0.1537 | 0.0812 | **False** |
| asym H₆ [0.9,0.9,2.0,0.9,0.9] | 2.4612 | +0.9945 | 0.0357 | True |
| asym H₆ [0.9,0.9,2.5,0.9,0.9] | 1.3517 | +1.1242 | 0.0368 | True |

**Spearman ρ(WG, b_min) = −0.867 (p = 0.001); Spearman ρ(overshoot, b_min) = −0.758 (p = 0.011).**
Both predictors anti-correlate with `b_min` in the expected direction (larger WG or larger overshoot
→ lower/easier `b_min`), but `|ρ(WG, b_min)| = 0.867 ≥ 0.7` **and** exceeds `|ρ(overshoot, b_min)| =
0.758`.

**G2 verdict: SUPPORTS the weighted-gap hypothesis** under the pre-registered rule. On this
10-witness set, the HF-weighted band-edge gap orders witness difficulty noticeably better than the
raw overshoot — consistent with, and a quantitative strengthening of, chem-pm3's own observation that
overshoot alone gets the ordering of its three witnesses wrong in places (R=1.1 is the smallest
overshoot yet the hardest witness; R=1.05 and R=1.15/R=1.20, added here, are further points the
overshoot-only story does not explain well but WG mostly does).

**Deviation from §3's stated plan (recorded honestly, not hidden):** §3 said each witness would get
"extended sweep to find `m_pass` if killed." The actual run passed `m_extended=()` to every
`run_witness` call — `m_pass` was **not** recomputed for the new witnesses, for compute-budget
reasons (an extended sweep to M≈60-64 means re-extending the Krylov basis past 24 for every one of
6 killed witnesses here, on top of the ~33 min already spent on the base sweep + search). G2's
*formal* acceptance criterion (§5) is stated in terms of `b_min`, which was computed in full, so this
deviation does not weaken the gate itself — it only means the richer "required depth" comparison
sketched in §3's prose was not run. Filed as follow-up scope, not done here.

### Part B — measured noise sensitivity (linear H₆ R=1.1, M=24, cap=1000, η_exact=0.0378)

40 trials per shot count, iid complex Gaussian perturbation of `moments[1:24]` at
`std = shot_noise_sigma(shots)`:

| shots | σ | mean(bound) | std(bound) | min | max | unsound % | worst-case `(2M−1)·cap·σ` | typical-case `‖x*‖₂σ√(2(M−1))` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1e2 | 0.1000 | −5028.41 | 837.77 | −7374.93 | −3501.66 | **100.0%** | 4700.00 | 3613.21 |
| 1e3 | 0.0316 | −1518.68 | 184.06 | −1966.50 | −1111.68 | **100.0%** | 1486.27 | 1142.60 |
| 1e4 | 0.0100 | −476.11 | 70.85 | −622.62 | −339.41 | **100.0%** | 470.00 | 361.32 |
| 1e5 | 0.0032 | −144.30 | 25.34 | −205.71 | −109.28 | **100.0%** | 148.63 | 114.26 |
| 1e6 | 0.0010 | −46.08 | 7.43 | −59.27 | −27.62 | **100.0%** | 47.00 | 36.13 |
| 1e8 | 0.0001 | −4.28 | 0.71 | −6.07 | −2.49 | **100.0%** | 4.70 | 3.61 |

**This is the headline, load-bearing finding of Part B, and it is worse than chem-pm3's analytic
estimate implied.** At *every* shot count tested — including 1e8 shots, far beyond any near-term
or even medium-term hardware budget — **100% of the 40 trials produced an unsound bound** (the
perturbed-moment LP value fell below the true η, frequently going negative, i.e. not merely loose
but outright non-physical as an upper bound on a quantity that must be ≥ 0). The empirical mean and
std track chem-pm3's own analytic worst-case slope `(2M−1)·cap·σ` to within roughly a factor of 1-2
at every shot count (e.g. at 1e8 shots: empirical mean magnitude 4.28 vs. worst-case 4.70 vs.
typical-case 3.61) — so that worst-case estimate, far from being a conservative/loose upper bound,
is a *reasonable estimate of the typical case*. The soundness proof for `_build_moment_bound_lp`
assumed exact moments; this measurement shows that assumption is not a formality to wave past —
under the measured noise model, with `coeff_cap=1000` and `M=24`, the bound has **no regime in the
tested range where it is usable**, and shrinking shot noise by six orders of magnitude (1e2 → 1e8)
only shrinks the typical swing proportionally, never closing the gap enough to recover soundness
with any of the margin this witness needs (η̂ would need the perturbed bound to land within
~0.01-0.02 of η_exact = 0.0378 to be useful, and even at 1e8 shots the swing is still ~100× that).
**Caveat stated plainly, not softened:** this bound, as currently built (this `coeff_cap`, this `M`),
is a noiseless-regime-only certificate; extending it to any noisy setting is open work, explicitly
out of scope here (§7).

### Overall verdict

- **G1:** PASS — all 10 witnesses (3 original + 7 new) satisfy bound ≥ η and non-increasing-in-M,
  exactly as chem-pm3's own gates require. Reused `_build_moment_bound_lp` unchanged.
- **G2:** computed and recorded; **substantively SUPPORTS** the weighted-gap hypothesis
  (ρ_WG = −0.867 vs. ρ_overshoot = −0.758, n=10).
- **G3:** computed and recorded; the honest finding is that the bound is **essentially always
  unsound** under the measured noise model at every tested shot count, 100% of trials, at a scale
  matching chem-pm3's own worst-case analytic estimate rather than beating it.
- `uv run pytest tests/test_pm3_depth_predictor_spec.py -v` → **3 passed** in 1218.16s.
- `uv run pytest tests/test_pm3_eta_bound_spec.py tests/test_subspace_floor_resolvability_spec.py -q`
  → **47 passed** (regression, unaffected).
