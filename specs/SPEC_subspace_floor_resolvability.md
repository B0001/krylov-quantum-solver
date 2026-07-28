# SPEC: The self-mode subspace floor is not rigorous for d ≥ 2 — a discovered failure, a heuristic pre-filter, and its falsification

**Status:** Updated 2026-07-28 after a parallel falsification sweep. The self-mode floor is a
**heuristic**; the disjoint-interval guard reduces but does **not** eliminate the failure (proven
insufficient below). **Oracle mode is the only rigorous path.** Corrects the earlier "fail-safe"
framing of this same spec.
**Depends on:** `hf_overlap_subspace.py`, `certified_overlap` (the block certificate).

## The finding (killable): the self-mode floor is not rigorous

`certify_hf_subspace_overlap` in self mode floors E_d by `theta_d − sigma_d` — the Weinstein floor
on the (d+1)-th Ritz value, generalized from the d = 1 `gap_bracket`. **This floor can exceed the
true reachable E_d**, so a non-vacuous self-mode certificate is not a proof. Linear H₆ (R = 1.2 Å,
d = 3): β_self = −2.251 > true E_d = −2.583 at M = 8 (−2.412 at M = 12). The excited Ritz states
are unresolved (σ₁…σ₃ ≈ 0.29–0.45), so `theta_d` sits near a higher level.

## The heuristic pre-filter (kept, but NOT sound)

The gross-failure case has an in-band signature: when the Krylov space badly fails to resolve the
cluster, the d + 1 lowest Weinstein intervals `[theta_k − sigma_k, theta_k + sigma_k]` **overlap**.
`certify_hf_subspace_overlap` self mode computes all d + 1 residuals and returns a **VACUOUS**
certificate when they overlap. This rejects the gross failures (in the stress sweep below, 106 of
114 cases), including the H₆ R=1.2 d=3 bug.

## The falsification: the guard is insufficient (killable)

A parallel adversarial sweep (linear H₆, symmetric and asymmetric, d = 3, M ∈ {6,8,10,12,16,20})
found **8 guard-PASSING cases whose self-mode floor is still invalid** (β_self > true E_d).
Representative, deterministic:

| case | intervals disjoint (guard passes) | β_self | true E_d | floor valid |
|---|---|---|---|---|
| linear H₆ R=1.0, d=3, M=16 | Yes | −2.218 | −2.500 | **NO** |
| linear H₆ R=1.1, d=3, M=20 | Yes | −2.369 | −2.584 | **NO** |
| asym H₆ [0.9,0.9,2.2,0.9,0.9], d=3, M=6 | Yes | −1.379 | −2.433 | **NO** |

Two facts make this decisive, not a corner case:
- **The blind spot is common, not pathological.** Every stretched H₆ geometry has a reachable
  level of HF amplitude ~1e-4 near the cluster boundary. Disjoint intervals localize d + 1
  *distinct* eigenvalues but not necessarily the d + 1 *lowest*; a missed low-amplitude level is
  localized as a higher one, with small residuals (disjoint intervals) all the way.
- **More Krylov dimension does not fix it** — escapes occur at M = 16 and M = 20.
- No local residual/disjointness criterion can detect a *missed* level; a sound test needs a
  certified lower anchor below the cluster (the `certified_gaps` open problem).

In the sweep the resulting *certificate* γ_min stayed ≤ exact ‖P_S u‖ in all 8 cases — but only by
unquantified slack, as thin as **0.0045** (asym H₆, M=6). So the library did not emit an outright
invalid number here, yet self-mode rests on slack it does not control. **It is not rigorous.**

## Conclusion / rigor statement

- **Oracle mode (a true E_d floor) is rigorous** — block sin-θ, validated with zero escapes across
  36k adversarial synthetic certificates and the molecular sweep (G4b, and the SPEC-21/21b gates).
- **Self-mode d ≥ 2 is HEURISTIC.** The guard is a pre-filter that rejects gross non-resolution,
  not a soundness proof. Docstrings say so; a non-vacuous self-mode certificate must not be quoted
  as certified. The d = 1 path was already "not self-verifiable" (gated empirically at M ≥ 6); this
  spec establishes the same, more sharply, for d ≥ 2.

## Gates

- G1 (the floor bug, regression, killable): raw `theta_d − sigma_d` > true E_d on H₆ R=1.2 d=3.
- G2 (pre-filter rejects gross failure): guarded self-mode is VACUOUS on H₆ R=1.2 d=3.
- G3 (no over-rejection on resolved cases): square-H₄ / linear-H₄ d=2 pass, non-vacuous & valid
  (regression-guards the PR #20 demonstration, which used resolved cases).
- G4 (**the falsification, killable**): a deterministic guard-PASSING case with an invalid floor
  exists (linear H₆ R=1.0, d=3, M=16) — documents the guard is not sound; flips if a future change
  makes it sound.
- G4b (the rigorous path, zero-tol): oracle-mode certificates are valid across the sweep.

## Out of scope

A proven-sufficient resolvability certificate (needs a certified sub-cluster anchor — deferred,
hard). Making self-mode rigorous. Guarding the d = 1 paths (unchanged; empirically M ≥ 6).

## Honest caveats

Exact-statevector only. The self-mode floor is heuristic; the guard is a heuristic pre-filter,
demonstrably insufficient. The block-bound math and oracle mode are the rigorous, validated core.
