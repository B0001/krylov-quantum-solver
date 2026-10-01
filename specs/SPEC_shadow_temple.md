# SPEC: Classical shadows estimate ⟨H⟩ and ⟨H²⟩ from one dataset — but do they buy Temple coverage?

**Status:** CLOSED — gates G1–G4 PASS (2026-10-01). Backlog hypothesis: specs/BACKLOG.md
"Method rungs", *"Classical shadows buy ⟨H⟩ and ⟨H²⟩ from one dataset — cheaper than the repo's λ²
model, but they do not buy coverage"* (chem-loe). **Both halves of the claim confirmed, plus an
unpredicted third finding:** shadows beat the λ² model's variance for both H and H² (G2/G3); shadow-
estimated Temple coverage lands in the same broken ~0.40 regime as the λ²-model, not materially
higher (G4 — the headline claim survives); and the HKP-style `shadow_norm` bound itself, while valid
for H, is empirically **violated** for H² (G3 — unpredicted, see §5/§8).

---

## 1. Goal

`SPEC_certified_noise` found λ_{H²} ≫ λ_H (H4: 62.9 vs 10.3 at the time of that spec; `λ_H2_bridge`'s
later identity-term correction moves both numbers — this spec measures λ_H=7.66, λ_{H²}=55.23,
ratio 7.21 — but the qualitative finding, λ_{H²} ≫ λ_H, is unchanged and in fact sharpens slightly),
making the Temple lower bound the noise-expensive side of the certified arc, and `temple_bounds`
states outright that the ⟨H²⟩
hardware cost is unmodeled. Random-Pauli classical shadows (`classical_shadows.py`, closed spec)
estimate **both** moments from the **same** measurement snapshots. Claim: shadows give a lower
variance bound than the repo's assumed λ² model for both ⟨H⟩ and ⟨H²⟩ — but they do **not** remove
the noisy-Temple coverage break `SPEC_certified_noise` found; feeding shadow-estimated moments into
the Temple bracket still lands coverage near the same ~0.40, because the break is a structural
variational knife-edge (ρ₀ → E₀), not an artifact of the λ²-Gaussian noise model. The claim is false
if shadow-estimated coverage is materially higher than the λ²-model coverage (a real "free lunch"),
or if the shadow variance bound is *not* below the λ² model (no efficiency win at all).

## 2. Background and honest framing

- **Prior art.** Composes two closed specs: `SPEC_classical_shadows` (the unbiased random-Pauli
  estimator, HKP shadow norm) and `SPEC_certified_noise` (the λ²-Gaussian noise model, the ~0.40
  coverage break). No new measurement primitive — this is a bridge, like `SPEC_lambda_h2_bridge`.
- **What we can claim if gates pass.** The shadow estimator is unbiased on the Krylov Ritz state
  (not just HF/FCI) for both ⟨H⟩ and ⟨H²⟩; its variance (bound *and* measured empirically) sits
  below the λ²-model's assumed variance for both moments; and the H²-vs-H cost asymmetry survives
  under shadows (shrunk, not removed). Whether shadow-estimated-moment Temple coverage matches or
  beats the λ²-model coverage is an open, gated question — not asserted in advance.
- **What we cannot claim.** (a) Everything `SPEC_classical_shadows` already disclaims: exact
  statevector simulation, no quantum advantage, random-Pauli only (no grouped/derandomized
  shadows). (b) A converged-depth result only — the Ritz state at other M is out of scope. (c) The
  comparison is single-molecule (H4, the system `SPEC_certified_noise` recorded the λ ratio on);
  generalizing the coverage-doesn't-improve finding to other systems is a follow-up.

## 3. Approach

Reuse `classical_shadows.collect_classical_shadow` / `shadow_energy_samples` / `shadow_norm`
unchanged (they already accept any `SparsePauliOp`, so H² needs no new estimator code — just
`(H @ H).simplify()`, the same expansion `certified_noise.hamiltonian_one_norms` already builds).
For a given Krylov Ritz state |Ψ₀(M)⟩:

1. Collect ONE shadow dataset (bases, signs); estimate ⟨H⟩ and ⟨H²⟩ from the **same** snapshots by
   calling `shadow_energy_samples` once with `H` and once with `H²`.
2. Compare `shadow_norm(H)` / `shadow_norm(H²)` (HKP bound) against `λ_H²` / `λ_{H²}²` (the repo's
   assumed per-shot variance under the λ² model) — **and** measure the empirical single-shot
   variance of the shadow estimator directly (`.var()` over snapshots), so the comparison is not
   two different upper bounds talking past each other.
3. Monte-Carlo over `trials` independent shadow datasets (fresh RNG seed each), forming the noisy
   Temple bracket `[τ₀, ρ₀]` from each trial's shadow-estimated ⟨H⟩, ⟨H²⟩, and checking coverage of
   the exact reachable E₀ (`reachability.reachable_eigenpairs`) — the exact empirical analogue of
   `certified_noise.shot_noise_coverage`, with real shadow sampling in place of the assumed Gaussian.

**Reference:** exact `⟨Ψ₀|H|Ψ₀⟩`, `⟨Ψ₀|H²|Ψ₀⟩` via `temple_bounds.mean_and_variance` (one sparse
matvec, the repo's existing convention) for unbiasedness; exact reachable `E₀` via
`reachability.reachable_eigenpairs` for coverage; `certified_noise.shot_noise_coverage`'s λ²-model
coverage as the baseline the shadow-based coverage is checked against.

## 4. Public interface

```
shadow_temple.h2_operator(mh) -> SparsePauliOp                            # (H@H).simplify(), electronic frame
shadow_temple.ritz_state(mh, m, solver=None) -> (psi, solver)
shadow_temple.shadow_moments(psi, n_qubits, op, op2, n_shots, seed=0) -> (per_shot_H, per_shot_H2)
shadow_temple.shadow_vs_lambda_model(mh, psi, n_shots=16000, seed=0) -> dict  # bound + empirical var, H and H2
shadow_temple.shadow_temple_coverage(mh, m, n_shots, trials=200, seed=0, solver=None) -> dict
```

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_shadow_temple_spec.py` (test-first). Seeded; H4/STO-3G (the
`SPEC_certified_noise` case with the recorded λ_{H²} ≫ λ_H boundary, 8 qubits, H: 185 / H²: 1775
Pauli terms — the O(N⁸) blowup recorded honestly, not hidden). PySCF/qiskit only, no block2.

- **G1 — unbiased on the Krylov Ritz state.** For the M=12 Ritz state (not HF, not FCI — the state
  the certified arc actually runs on), the shadow mean of **both** ⟨H⟩ and ⟨H²⟩ is within
  `4·stderr` of the exact value at 16k snapshots.
- **G2 — shadow bound beats the λ² model.** `shadow_norm(H²) < λ_{H²}²` **and**
  `shadow_norm(H) < λ_H²` (dies if either fails — no efficiency win). Also records the ratio
  `shadow_norm(H²)/shadow_norm(H)` against `(λ_{H²}/λ_H)²`: the claim is that shadows shrink the
  H²-vs-H asymmetry, not remove it — dies if the shadow ratio is not strictly smaller than the
  λ²-based ratio, and dies (as a separate, stronger failure) if the shadow ratio is not itself > 1
  (i.e. shadows accidentally claiming H² is *cheaper* than H would itself be too surprising to
  accept silently). Units note: `shadow_norm` and `λ²` are both variance-scale quantities, so the
  comparison ratio must be `(λ_{H²}/λ_H)²`, not the bare amplitude ratio `λ_{H²}/λ_H` — an earlier
  draft of the gate compared the unsquared ratio and failed as a result (measured `λ`-ratio 7.21,
  `λ²`-ratio 52.0, shadow ratio 29.7 — inside the squared bracket, outside the unsquared one). The
  backlog's own recorded "λ-based ratio of 35.6" is itself a λ²-scale number
  (`(62.9/10.3)² ≈ 37.3` at the scout-probe's pre-correction λ values), confirming the squared
  reading was the one intended.
- **G3 — empirical variance, measured directly (not inferred from bounds alone).** The empirical
  single-shot variance of the shadow estimator (`.var()` over snapshots) is `< λ²` (beats the model
  it is being compared against) for **both** H and H² — the real efficiency claim, and survives.
  **Unpredicted, sharper finding:** the HKP-style additive `shadow_norm` bound itself (`Σ|c_k|²3^{w_k}`,
  validated on H alone by the closed `SPEC_classical_shadows`) holds empirically for H (ratio
  measured/bound 0.69–0.95 over 10 independent seeds) but is **violated, reproducibly, for H²**
  (ratio 1.16–1.49× over the same 10 seeds — not a one-off fluctuation). The naive per-term
  diagonal sum ignores positive cross-correlations between H²'s 1775 heavily support-overlapping
  Pauli terms, which are not negligible once term count and overlap grow this far (H²'s O(N⁸)
  expansion). `shadow_norm` therefore cannot be trusted as a literal shot-noise-variance bound for
  composite operators like H² — recorded here (gate pins the violation itself, not a loosened
  tolerance) rather than smoothed over.
- **G4 — the real test: Temple coverage under shadow-estimated moments.** 200-trial Monte Carlo
  coverage of the exact reachable E₀ using shadow-estimated ⟨H⟩, ⟨H²⟩ (fresh snapshots per trial),
  compared against `certified_noise.shot_noise_coverage`'s λ²-model coverage at a matched shot
  budget. Records whether shadow coverage is materially higher (a real free lunch — would falsify
  the headline claim) or lands within the same broken regime (~0.40, i.e. estimator-independent).

> Definition of done: **G4** — it is the question the bead asks. G1–G3 are necessary preconditions
> (an estimator that is biased, or not actually cheaper, would make G4 moot either way).

## 6. Implementation plan (test-first)

1. Write `tests/test_shadow_temple_spec.py` encoding G1–G4 (initially failing — module absent).
2. Add `shadow_temple.py`, composing `classical_shadows`, `certified_noise`, `temple_bounds`,
   `reachability` — no new measurement primitive.
3. Iterate to green via `uv run pytest tests/test_shadow_temple_spec.py` (own process; pyscf/qiskit,
   no block2). Record actual numbers (shots/trials, timing) in §9 / the PR description.

## 7. Out of scope

- Grouped/derandomized shadows (would shrink the 3^weight cost further — `SPEC_classical_shadows`
  already scoped this out).
- Any system other than H4 (the recorded λ_{H²} ≫ λ_H case); generalizing across bond lengths,
  active spaces, or molecules is a follow-up.
- Inflated (z·se) shadow-based brackets — `SPEC_certified_noise`'s inflation-restores-coverage
  result is not re-derived here; this spec only checks whether the *raw* break moves.
- Real hardware grouped-Pauli measurement covariances (both this spec and `SPEC_certified_noise`
  use the idealized i.i.d. per-operator model).

## 8. Caveats and risks

- **R1 — bound vs. bound is not evidence on its own.** `shadow_norm` and λ² are upper bounds from
  two different protocols; G3's direct empirical-variance measurement is the mitigation the bead
  calls for, not a formality. **It paid off**: the direct measurement is what caught R1's own risk
  materializing — `shadow_norm` is not actually a valid empirical bound for H², so a gate that only
  compared `shadow_norm` to λ² (never measuring the real variance) would have reported an efficiency
  win for H² that systematically understates the true per-shot cost by ~20–50%.
- **R2 — compute cost.** Each shadow snapshot rotates and samples a 256-dim statevector; 200 trials
  × n_shots is the dominant cost of G4. Mitigated by keeping n_shots modest (a few thousand) — the
  claim under test is about the *coverage regime*, not shot-optimality.
- **R3 — H² term-count blowup is real and unhidden.** 1775 Pauli terms vs H's 185 on H4 (O(N⁸));
  this makes shadow estimation of H² itself non-trivially expensive per snapshot, independent of
  the norm/variance comparison — recorded, not smoothed over.
- Honest limitation: single molecule, single Krylov depth, idealized shot noise — a boundary map,
  not a general theorem.

## 9. Deliverables

- `shadow_temple.py` — the four interface functions above.
- `tests/test_shadow_temple_spec.py` — gates G1–G4 (5 tests; all green).
- Results summary, H4/STO-3G, M=12 Krylov Ritz state, `uv run python shadow_temple.py`:
  - H has 185 Pauli terms, H² has 1775 (O(N⁸) blowup, recorded not hidden).
  - `λ_H=7.66, λ_{H²}=55.23` (λ²-ratio 52.0) — post `SPEC_lambda_h2_bridge` correction; see §1 note.
  - H: `shadow_norm=40.18 < λ²=58.63`; empirical variance ≈ 29.7–38.2 over 10 seeds (ratio to bound
    0.69–0.95) — **bound holds**.
  - H²: `shadow_norm=1193.38 < λ²=3050.85`; empirical variance ≈ 1384–1774 over 10 seeds (ratio to
    bound 1.16–1.49) — **bound measurably violated, reproducibly, while still beating λ²**.
  - shadow-norm ratio H²/H = 29.7, vs λ²-ratio 52.0: asymmetry shrinks (not removed), confirming the
    scout probe's qualitative direction (its specific numbers predate the λ_H2_bridge correction).
  - **G4, the real test** (seed 11, 200 trials × 4000 shots, 114 s): shadow-estimated-moment coverage
    `cov_raw=0.395, cov_upper=0.505`; λ²-model coverage at the matched shot budget (4000 trials):
    `cov_raw=0.385, cov_upper=0.487`. Difference 0.01 — **no free lunch**: shadows land in the same
    broken ~0.39-0.40 coverage regime as the λ²-Gaussian model, confirming the structural,
    estimator-independent reading of the variational knife-edge.
  - Hardware/timing: single-process `uv run pytest tests/test_shadow_temple_spec.py`, 5 tests,
    ~120 s wall time (8-core container, no GPU; see handoff for `uname -a`).
