# SPEC: the near-term/FT bridge prices only ⟨H⟩, and one call site never got the identity fix

**Status:** IMPLEMENTED once gates green. Backlog hypothesis: *"The bridge prices only ⟨H⟩, but the
method it represents needs ⟨H²⟩ — and one call site never got the identity fix"* (specs/BACKLOG.md,
Open → Fault-tolerant stack).

---

## 1. Goal

Two claims, one verified defect and one measurement.

**(a) The defect.** `SPEC_lambda_meas_identity` established that the identity term must not enter a
shot-noise 1-norm — it is a constant of zero variance and costs zero shots. That fix landed in
`precision_cost.measurement_lambda`, which now carries a full docstring about it. It **never reached
`certified_noise.hamiltonian_one_norms`**, five files away, which still computes
`np.abs(Hop.coeffs).sum()` on *both* λ_H and λ_{H²}. Measured identity fraction: **21.8–30.1%**.

**(b) The measurement.** `precision_cost.near_term_shots` charges the near-term side using the
1-norm of **H alone**. But the near-term certified method the bridge is *defined as* is the Temple
bracket, which needs ⟨H⟩ **and** ⟨H²⟩, and `SPEC_certified_noise` already recorded λ_{H²} ≫ λ_H
("the Temple lower bound is the noise-expensive side"). So the bridge understates near-term cost —
i.e. **understates FT's advantage**, the opposite direction from this arc's two prior
self-corrections (`SPEC_shift_both_sides`, `SPEC_lambda_meas_identity`, which both moved toward
near-term).

## 2. Background and honest framing

**What we can claim:** the identity fix has a missed call site; the resulting λ's are inflated by
~22–30%; the ⟨H²⟩ undercount is real and moves real `cost_advisor` verdicts.

**What we cannot claim — the naive size law is FALSE.** The backlog entry predicted the undercount
ratio λ_{H²}/λ_H "grows with size". Measured, it does **not** grow with qubit count across
heterogeneous active spaces: LiH CAS(2,3) at 6 qubits has ratio **1.29**, *below* H₂ CAS(2,2) at 4
qubits (**1.74**). The growth is real only within a homogeneous family — the Hₙ/STO-3G series gives
1.74 → 7.09 → 16.69 at 4/8/12 qubits. **The undercount cannot be predicted from system size
alone**, so `SPEC_precision_cost`'s "margin grows with size" headline is *not* cleanly cancelled;
it is confounded by active-space composition. G5 gates this boundary so the killed form of the
claim cannot quietly return.

**What we cannot claim — the shot-budget model.** How λ_{H²} enters the shot count is a **modelling
choice, not a derivation**. The Temple bracket's dependence on ⟨H²⟩ precision is not a plain
additive term. We therefore gate the *direction* and the *existence* of verdict movement under one
clearly-labelled inflation (λ_H + λ_{H²}, the summed-1-norm proxy), and explicitly do **not**
re-derive a crossover. A different allocation model would move the numbers; it would not restore
the missing ⟨H²⟩ cost to zero.

**Not claimed:** grouped-Pauli or commuting-set measurement (which would shrink λ_{H²}); any
statement about quantum advantage.

## 3. Approach

Fix `hamiltonian_one_norms` to exclude identity by default, mirroring `measurement_lambda`'s
existing API (`include_identity=True` reproduces the old inflated value for archaeology). Then
measure, against real SCDF-shifted λ_DF from CASCI integrals, whether pricing ⟨H²⟩ moves any
`cost_advisor` verdict over ρ ∈ [1, 1e6].

**Reference:** `precision_cost.measurement_lambda` is the convention oracle — the two modules must
agree exactly on λ_H (they currently do not). `SPEC_certified_noise`'s recorded λ_{H²} > λ_H
boundary is the finding that must survive the fix.

## 4. Public interface

```
certified_noise.hamiltonian_one_norms(mh, *, include_identity=False) -> (lambda_H, lambda_H2)
```

Signature is additive; the default flips to the corrected convention.

## 5. Acceptance criteria (validation gates)

`tests/test_lambda_h2_bridge_spec.py`.

- **G1 — the two modules agree on λ_H (DEFINITION OF DONE).** `hamiltonian_one_norms(mh)[0] ==
  measurement_lambda(mh)` exactly on every case. This is the defect: today they disagree by the
  identity fraction.
- **G2 — archaeology.** `include_identity=True` reproduces the pre-fix `np.abs(coeffs).sum()`
  values, so recorded numbers stay reachable.
- **G3 — the recorded boundary survives the fix.** λ_{H²} > λ_H still holds on every
  `SPEC_certified_noise` case after identity exclusion; that spec's G4 finding is not overturned.
  (H₄ moves 63/10 → 50.7/7.1; the ratio *strengthens* 6.3 → 7.09.)
- **G4 — growth within a homogeneous family.** For Hₙ/STO-3G, n = 2, 4, 6, the ratio is strictly
  increasing. Killed if flat or decreasing.
- **G5 — the boundary that kills the naive law.** Across heterogeneous CAS spaces the ratio is
  **not** monotone in qubit count: LiH CAS(2,3) (6q) < H₂ CAS(2,2) (4q). Killed if LiH's ratio
  exceeds H₂'s — which would mean size alone does predict the undercount after all.
- **G6 — the undercount is material, and its direction is one-way.** With real SCDF-shifted λ_DF,
  pricing ⟨H²⟩ moves at least one verdict, and **every** move is near-term → FT. Killed if no
  verdict moves (the undercount would be immaterial to the advisor), or if any move runs the other
  way.

## 6. Implementation plan (test-first)

1. `tests/test_lambda_h2_bridge_spec.py` encoding G1–G6 (initially failing).
2. `certified_noise.hamiltonian_one_norms`: add `include_identity`, default `False`.
3. `make gates`.

## 7. Out of scope

- Changing `near_term_shots` / `cost_advisor` to charge ⟨H²⟩. This spec **measures** the undercount
  and its direction; picking the allocation model is a follow-up that needs the Temple variance
  decomposition, not a summed 1-norm.
- Grouped-Pauli / commuting-set measurement (would shrink λ_{H²}).
- Re-deriving any published crossover.

## 8. Caveats and risks

- **R1 — the inflation λ_H + λ_{H²} is a proxy, not a derivation.** Labelled as such in G6; the
  gate asserts direction and existence of movement, never a crossover value.
- **R2 — λ_{H²} here is the ungrouped 1-norm** from `(H @ H).simplify()`, and H² has far more terms
  than H (O(N⁸) growth). Gates stay at ≤12 qubits.
- **R3 — G3 depends on `SPEC_certified_noise`'s case list**; if that spec's cases change, G3 must be
  re-derived rather than silently narrowed.
- The spec-prose numbers in `SPEC_certified_noise` (H₄ "63 vs 10") are identity-inclusive and become
  50.7/7.1 after the fix. Recorded here rather than edited away.

## 9. Deliverables

- `certified_noise.py` — `hamiltonian_one_norms(mh, *, include_identity=False)`.
- `tests/test_lambda_h2_bridge_spec.py` — G1–G6.
