# SPEC: Tensor hypercontraction — a rank-compact ERI factorization with a validated qubitization λ

**Status:** CLOSED — gates G1–G4 PASS (2026-06-29); `thc_factorization.py` merged. Reconstruction
exact (~5e-13), FCI preserved (<1e-7 mHa), `thc_lambda` == `df_lambda` on the structured THC; rank
28<196 (H₂O) / 55<530 (N₂). **Finding:** unoptimized collocation gives λ_THC ≈ 5377 vs λ_DF ≈ 87
(H₂O, ~62×) — the λ advantage needs ISDF/optimized points (out of scope, §7).

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

Add a tensor-hypercontraction (THC) factorization of the active-space ERIs,
`(pq|rs) ≈ Σ_{μν} χ_p^μ χ_q^μ ζ_{μν} χ_r^ν χ_s^ν`, with its qubitization 1-norm `λ_THC` computed in
the **same convention** as the double-factorization `df_lambda`. Two falsifiable claims: (a) a linear
least-squares THC reconstructs the ERIs *exactly* and preserves FCI at THC rank
`M = norb(norb+1)/2`, compressing the rank relative to the DF-derived THC; and (b) the `thc_lambda`
formula is correct — it reproduces `df_lambda` to machine precision on the DF-derived (structured)
THC. The **honest boundary** this spec also pins: with *unoptimized* collocation, `λ_THC` is far
*larger* than `λ_DF`; the literature THC λ advantage requires optimized (ISDF) collocation, which is
out of scope. The claim set is false if reconstruction is not exact at full pair rank, if the λ
formula disagrees with `df_lambda` on the structured THC, or if the rank is not compressed.

## 2. Background and honest framing

- **Prior art / reference.** Lee, Berry, Gidney, Huggins, McClean, Wiebe, Babbush, *Even more
  efficient quantum computations of chemistry through tensor hypercontraction*, PRX Quantum 2,
  030305 (2021); the DF reference already in `df_factorization.py` (von Burg/Lee 2021). THC is the
  asymptotically cheaper qubitization route at scale.
- **What we can claim if the gates pass.** A correct, exact THC factorization of the ERIs; a
  qubitization λ that is validated against the existing DF λ; and a quantified rank compression and
  λ-advantage *boundary* on small molecules.
- **What we cannot claim (stated up front).** (a) **No λ advantage at this scale with this method** —
  unoptimized collocation gives `λ_THC ≫ λ_DF` (measured ≈ 60× on H₂O/STO-3G). The literature λ
  advantage needs ISDF/optimized points (a nonlinear fit), explicitly out of scope (§7). (b)
  Reproduction, not novelty. (c) Minimal-basis active spaces; not an FT resource claim — it feeds
  `ft_resource_estimator.py` only once an optimized-collocation THC exists.

## 3. Approach

- **Structured THC from DF (the anchor).** Substitute each DF leaf
  `L^t_pq = Σ_k w^t_k u^tk_p u^tk_q` into `(pq|rs) = Σ_t g_t L^t_pq L^t_rs`: this *is* a THC with
  `χ` columns `{u^tk}` and `ζ` block-diagonal (`ζ_{(t,k),(t,l)} = g_t w^t_k w^t_l`). Exact by
  construction; `thc_lambda` on it must equal `df_lambda`.
- **Linear-LS THC (the compression).** Pick a random full-rank collocation `χ` (norb × M), solve the
  central matrix `ζ` by linear least squares (`ζ = pinv(Pᵀ) V pinv(P)`, `P_{μ,(pq)} = χ_p^μ χ_q^μ`).
  At `M ≥ norb(norb+1)/2` the pair space is spanned ⇒ exact reconstruction, no nonlinear iteration.
- **λ.** `λ_THC = ¼ Σ_{μν} |ζ_{μν}| ‖v^μ‖² ‖v^ν‖² + ‖h1‖_nuc`, with `‖v^μ‖² = Σ_p (χ_p^μ)²` (each THC
  density is rank-1, so `Σ_k|eig_k| = ‖v^μ‖²`) — the `df_lambda` convention specialized to THC.
- **Reference.** `df_factorization.df_lambda` / `double_factorize` (λ + structure anchor) and PySCF
  FCI on the reconstructed ERIs (energy).

## 4. Public interface

```
thc_factorization.tensor_hypercontraction(eri, norb, n_thc=None, seed=0) -> (chi, zeta)
thc_factorization.thc_from_df(eri, norb, rank=None)                      -> (chi, zeta)
thc_factorization.reconstruct_thc(chi, zeta)                            -> eri
thc_factorization.thc_lambda(chi, zeta, h1)                             -> float
thc_factorization.thc_rank(norb)                                        -> int   (norb(norb+1)/2)
```

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_thc_lambda_spec.py` (test-first). PySCF/NumPy only (no block2). Small CAS so
the dense ERI tensor and FCI are cheap.

- **G1 — exact reconstruction + FCI (definition of done).** On H₂O and N₂ (STO-3G full space),
  `tensor_hypercontraction` at `M = norb(norb+1)/2` reconstructs the ERIs to `< 1e-9` and the FCI
  energy on the reconstructed ERIs matches the exact FCI to `< 1e-6` Ha.
- **G2 — λ machinery validated against DF.** `thc_lambda(thc_from_df(eri, norb), h1)` equals
  `df_lambda(double_factorize(eri, norb), h1)` to `< 1e-9` on H₂O and N₂. (This anchors the THC
  1-norm to already-validated code, independent of any remembered literature constant.)
- **G3 — rank compression.** `thc_rank(norb) = norb(norb+1)/2` is strictly less than the DF-derived
  THC rank `norb × (DF full rank)` (e.g. 28 < 196 for H₂O, 55 < 530 for N₂).
- **G4 — λ correctness + recorded boundary.** `thc_lambda` is invariant under `χ → cχ` (scale
  invariance of a 1-norm, `< 1e-6` relative) **and** the recorded finding holds: with unoptimized
  collocation `λ_THC > λ_DF` (measured ≈ 60×). This locks the boundary; a collocation method that
  *beat* DF would flip it and require revising this gate (which is exactly the finding to record).

> Definition of done: **G1 + G2**. G4's inequality is a *finding lock*, not an aspiration — the
> honest content is "naive THC does not beat DF λ here; the advantage needs ISDF." If a future
> optimized-collocation THC breaks G4, revise it and record the crossover.

## 6. Implementation plan (test-first)

1. Write `tests/test_thc_lambda_spec.py` encoding G1–G4 (initially failing).
2. Add `thc_factorization.py` (reusing `df_factorization.double_factorize` for the anchor).
3. Iterate to green via `make gates` (own process; PySCF/NumPy, no block2).

## 7. Out of scope

- **Optimized / ISDF collocation and nonlinear (ALS/CP) THC refinement** — the route to a genuinely
  *small* λ_THC. This is the hard, research-grade part and the reason THC helps at scale; pinned as
  the boundary, not attempted here.
- Sub-pair-rank compression (`M < norb(norb+1)/2`), which also needs the nonlinear fit.
- Feeding λ_THC into `ft_resource_estimator.py` for a Toffoli/T-gate count (meaningful only once an
  optimized-collocation THC exists).

## 8. Caveats and risks

- **R1 — naive λ is large.** The headline risk *is* the finding: unoptimized collocation gives
  `λ_THC ≫ λ_DF`. Mitigation: state it up front (§2), gate it as a recorded boundary (G4), and keep
  the optimized route out of scope.
- **R2 — random collocation conditioning.** A pathological `χ` could be rank-deficient. Mitigation:
  seed the RNG; `M = norb(norb+1)/2` random Gaussian columns span the pair space with probability 1,
  and G1 would fail loudly otherwise.
- Honest limitation: minimal-basis, full-space tiny molecules; a correctness/■-norm study, not an FT
  resource result.

## 9. Deliverables

- `thc_factorization.py` — `tensor_hypercontraction`, `thc_from_df`, `reconstruct_thc`,
  `thc_lambda`, `thc_rank`.
- `tests/test_thc_lambda_spec.py` — gates G1–G4.
- Results summary (with the §2/§7 caveats and the measured λ_THC/λ_DF ratio) in the PR description.
