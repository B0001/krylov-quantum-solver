# SPEC: Self-mode subspace floor is not rigorous for d ≥ 3 — a discovered failure and a fail-safe guard

**Status:** Specced 2026-07-27. A soundness fix to `hf_overlap_subspace.py`
(merged PR #20). Deferred boundary named in `SPEC_hf_overlap_subspace.md` honest-caveats:
"a general certified resolvability test for d > 1 is deferred."
**Depends on:** `hf_overlap_subspace` (the self-mode E_d floor), `certified_overlap` (the block
certificate it feeds).

## The finding (killable)

`certify_hf_subspace_overlap` in self mode floors E_d by `theta_d − sigma_d` — the Weinstein
floor on the (d+1)-th Ritz value, generalized from the d = 1 `gap_bracket`. **This floor is not
rigorous for d ≥ 3.** On the linear H₆ chain (R = 1.2 Å, d = 3):

| M | β_self = θ_d − σ_d | true reachable E_d | β_self ≤ E_d ? |
|---|---|---|---|
| 8  | −2.25106 | −2.58257 | **NO** |
| 12 | −2.41174 | −2.58257 | **NO** |

β_self **exceeds** the true (d+1)-th reachable level, so it is not a valid lower bound, and the
resulting block certificate γ_min (0.88) came out *above* the oracle-mode γ_min (0.53) — a
near-miss from an outright invalid certificate (it stayed under the exact ‖P_S u‖ = 0.92 only by
slack). The cause is visible in the residuals: the excited Ritz states are unresolved at these M
(σ₁…σ₃ ≈ 0.29, 0.34, 0.45), so `theta_d` has not converged to E_d — it sits near a higher level.

The d = 1 (gap_bracket) and d = 2 (square-H₄ / linear-H₄) paths did not exhibit this in testing;
d = 3 on a genuinely multireference chain is where the unresolved-excited-state regime bites.

## The fix: a checkable resolvability signal, fail-safe

The failure has an **in-band signature**: when the Krylov space has not resolved the cluster, the
Weinstein intervals `I_k = [theta_k − sigma_k, theta_k + sigma_k]` of the first d + 1 Ritz states
**overlap**. When it has resolved them, the intervals are disjoint. Observed:

| case | intervals disjoint | β_self valid |
|---|---|---|
| H₆ R=1.2 d=3 M=8/12 (fails) | **No**  | No  |
| square-H₄ d=2 (good)         | Yes | Yes |
| linear-H₄ d=2 (good)         | Yes | Yes |

`certify_hf_subspace_overlap` self mode now computes all d + 1 Ritz residuals and checks pairwise
disjointness of the Weinstein intervals. **If they overlap, it returns an explicit VACUOUS
certificate** (reason: self-mode floor unresolved; increase M or supply an oracle e_d) rather than
the possibly-unsound positive. Self mode is now *fail-safe* (refuses when unresolved), not
*fail-silent*. Oracle mode is unchanged and remains the rigorous path.

## Honest scope of the guard

Disjoint Weinstein intervals is a **necessary, in-band resolvability signal that is checkable
without an oracle** — the property the existing certified line calls "not self-verifiable" and
gates only empirically at M ≥ 6. It is **not proven sufficient**: a pathological reachable level
of vanishing HF amplitude near the cluster boundary could in principle leave the intervals
disjoint while a level is mislocalized. The guard converts the demonstrated *fail-silent* into a
*fail-safe*; it does not turn self mode into a rigorous certificate. Oracle mode remains the
proof; the empirical soundness sweep (G4) is exactly the standard the d = 1 premise is held to
(zero escapes across a family, honestly labeled empirical).

## Gates

- G1 (the bug, regression, killable): on H₆ (R = 1.2 Å, d = 3, M ∈ {8, 12}) the *raw* self-mode
  floor `theta_d − sigma_d` exceeds the true reachable E_d. Documents the failure exists; if a
  future change makes the raw floor rigorous here, this gate flips and the spec is revisited.
- G2 (guard catches it): the guarded `certify_hf_subspace_overlap` returns a **VACUOUS**
  certificate on H₆ d = 3 (M ∈ {8, 12}) — never the unsound positive.
- G3 (guard does not over-reject): on square-H₄ (a ∈ {1.4, 1.2, 1.0}, d = 2) and linear-H₄
  (R = 1.0, d = 2) the guard passes and the certificate is non-vacuous and valid (γ_min ≤ exact).
  Regression-guards the merged PR #20 demonstration.
- G4 (empirical soundness, zero-tol): across a sweep (square-H₄, linear-H₄ R ∈ {1.0,1.5},
  linear-H₆ R ∈ {1.2,1.8}, d ∈ {2,3}, M ∈ {6,8,12}), **every** case the guard passes has a valid
  self-mode floor (β_self ≤ true reachable E_d) and a valid certificate (γ_min ≤ exact ‖P_S u‖).
  One guard-passes-but-invalid escape kills the guard.

## Out of scope

A proven-sufficient resolvability certificate (the sound sufficient condition, requiring a
certified anchor below the cluster — deferred, hard, the certified_gaps docstring's open
problem). Guarding the d = 1 `gap_bracket` / `hf_overlap_certificate` self-mode paths (empirically
safe at M ≥ 6; unchanged here). Shot noise on the residuals.

## Honest caveats

Exact-statevector only. The guard is a necessary check, not a sufficiency proof (see scope). The
G4 soundness claim is empirical over the tested family, in the same spirit as the d = 1 M ≥ 6
boundary — not a theorem.
