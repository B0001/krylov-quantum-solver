# SPEC: the chained overlap bound — the stub, implemented, and it works in self mode

**Status:** IMPLEMENTED once gates green. Closes the backlog hypothesis *"The `krylov_refine` stub is
not a marginal tightening"*. Unblocked by [`SPEC_reachability_tolerance.md`](SPEC_reachability_tolerance.md)
and [`SPEC_symmetry_reachability.md`](SPEC_symmetry_reachability.md), which made the certified target
well defined.

---

## 1. Goal

`hybrid_quantum_solver/certified_overlap/krylov_refine.py` has been a `NotImplementedError` stub
since 2026-07-17, with a docstring promising Lanczos-chained refinement. Implement it.

The direct SPEC-21 certificate bounds |⟨u|ψ₀⟩| using **u's own** residual, so it goes vacuous when
r_u ≥ δ_u — which for a Hartree–Fock guiding state happens exactly on the multireference systems the
certificate is wanted for. Chaining through the Krylov ground Ritz vector v uses **v's** residual
instead:

    θ(u,ψ₀) ≤ θ(u,v) + θ(v,ψ₀) ≤ θ(u,v) + arcsin(r_v/δ_v)
    ⇒ |⟨u|ψ₀⟩| ≥ cos( θ(u,v) + arcsin(r_v/δ_v) )

valid while the angle sum stays below π/2. It costs **no extra measurements — but the justification
has two halves, not one.** ⟨u|v⟩ is a linear combination of Krylov overlap-matrix entries the solver
has already formed; r_v additionally needs ⟨v|H²|v⟩, which `certified_gaps.gap_bracket` already
computes as its variance term. So the claim holds *given that the gap bracket has been run* — which
it has, since it supplies `e1_floor`. Quoting only the ⟨u|v⟩ half would overstate it.

**The claim: this is valid, strictly tightens the direct bound, rescues every vacuous case, and —
the part that decides whether it is usable in production — survives with the repo's own self-mode
E₁ floor instead of an oracle.**

## 2. Background and honest framing

Measured on a symmetry-clean set (**5 systems × M ∈ {6,8,12} = 15 cells**), all references built at
`conv_tol=1e-13` so the SCF-residue artifact of `SPEC_reachability_tolerance` is absent. (An earlier
draft quoted 7 systems / 21 cells. That set included square H₄ at a = 1.10 and 1.35 — the two
geometries this spec itself excludes as having no well-defined target, see R3. The counts below are
the post-exclusion re-run; the *mechanism* numbers were unaffected by the correction.)

- **Validity: 0 violations in 15 cells**, oracle and self mode alike. Largest excess over the exact
  reference **+2.22e-16** — pure float rounding on a saturated case.
- **Rescue: 6/15 cells have a VACUOUS direct bound** (all of linear H₆, all of square H₄). The
  chained bound rescues **6/6 in self mode as well as oracle** — identical coverage.
- **Self mode is nearly free.** γ_self/γ_oracle: min 0.9221, median 1.0000, max 1.0000. The cost
  concentrates at M = 6–8 on square-H₄ near-degeneracies and vanishes by M = 12.

**The mechanism worth recording.** Where the *direct* bound survives, a loose self-mode floor costs
it heavily — linear H₄ at M=6: **0.7765 (oracle) → 0.4791 (self), a 38% loss**. The *chained* bound
at the same point loses **0.08%** (0.9655 → 0.9647). The floor enters the chained bound only through
`arcsin(r_v/δ_v)` with a tiny r_v, not through the much larger HF residual, so **chaining absorbs a
loose E₁ floor far better than the direct bound does.** That is the practical argument for it, and
it is stronger than the tightening headline.

**What we cannot claim.**
- **It is NOT monotone in M.** On square H₄ the bound gets *worse* from M=6 to M=8 before jumping at
  M=12. Any claim of monotone improvement with depth is false; G5 pins this.
- **It SATURATES.** At a machine-converged Ritz vector the bound *equals* the exact overlap, and
  rounding puts it 1–3 ulp either side. A bare `assert chained <= exact` fails on correct behaviour;
  gates must carry `SATURATION_SLACK = 1e-14`. This is a real property, not a fudge.
- **It does NOT enforce the M ≥ 6 premise gate.** `certify_hf_overlap` hard-raises below M = 6 in
  self mode; `refine_via_lanczos` takes no `m` and enforces nothing — the caller owns the floor's
  provenance, and a floor derived below M = 6 will be consumed without complaint. That is a
  deliberate difference from its sibling, stated here because an earlier draft of this spec claimed
  the gate was inherited. It is not. The certificate remains conditional on whatever premise backs
  the floor, which is checkable but not self-verifiable.
- Reachable-sector scope and exact-statevector evaluation, inherited from the direct bound.

**What it does NOT do:** it does not moot the d=2 block certificate. The two bound *different*
quantities — |⟨u|ψ₀⟩| versus ‖P_S u‖, with ‖P_S u‖ ≥ |⟨u|ψ₀⟩| always — so "chained d=1 beats block
d=2" is a category error. That comparison is dropped, not resolved.

## 3. Approach

Compose validated primitives: the Ritz vector from `QuantumKrylovSolver.eigenstates`, the E₁ floor
from `certified_gaps.gap_bracket` (oracle or self mode, caller's choice and caller's provenance).
No new numerics beyond the angle arithmetic.

**Reference:** `hf_overlap_certificate.exact_reachable_overlap` (dense) — the killable check. Any
γ_chain above it, beyond saturation slack, kills the bound.

## 4. Public interface

```
certified_overlap.krylov_refine.refine_via_lanczos(H, u, v, e1_floor) -> float | None
certified_overlap.krylov_refine.SATURATION_SLACK = 1e-14
```

`None` means VACUOUS — no statement — never a fabricated zero-information positive.

## 5. Acceptance criteria (validation gates)

`tests/test_chained_overlap_spec.py`.

- **G1 — validity (DEFINITION OF DONE).** γ_chain ≤ exact + `SATURATION_SLACK` on every system × M,
  in **both** oracle and self mode. One real violation kills it.
- **G2 — it tightens.** Wherever the direct bound is non-vacuous, γ_chain > γ_direct strictly.
- **G3 — it rescues.** On linear H₆ and square H₄, where the direct bound is vacuous, γ_chain is
  non-vacuous — in self mode, not just oracle.
- **G4 — self mode is nearly free, and absorbs the floor better than the direct bound.**
  γ_self/γ_oracle ≥ 0.9 everywhere; and on linear H₄ at M=6 the direct bound's oracle→self loss
  exceeds 30% while the chained bound's is under 1%.
- **G5 — NOT monotone in M.** A geometry exists where γ_chain decreases from M=6 to M=8. Killed if
  the bound turns out monotone — the caveat would be unnecessary and should be removed.
- **G6 — vacuous is `None`, and inputs are checked.** Unnormalized `u` or `v` raises; a floor below
  λ_v returns `None`.
- **G7 — R2b is measured, not asserted.** On the gated set, ‖P_unreach v‖ < 1e-12 (so the neglected
  O(leak²) term stays ≥ 10 orders below `SATURATION_SLACK`), *and* unreachable levels really do sit
  below E₁_reach so the caveat is non-vacuous. Positive control: the excluded a = 1.10 / 1.35
  geometries must FAIL that leakage threshold — if they passed, the exclusion would be unnecessary
  and R2b overstated.

## 6. Implementation plan (test-first)

1. `tests/test_chained_overlap_spec.py` encoding G1–G7 (initially failing — stub raises).
2. Implement `refine_via_lanczos`.
3. `make gates`.

## 7. Out of scope

- Wiring the chained bound into `certify_hf_overlap` as the default. It changes recorded γ values in
  `SPEC_hf_overlap_certificate`; a follow-up should do that deliberately with the numbers re-run.
- The d=1 vs d=2 comparison (category error, see §2).
- Shot noise on r_v / ⟨u|v⟩ — the direct bound's noise sibling does not exist either.

## 8. Caveats and risks

- **R1 — saturation makes a naive gate wrong.** Mitigated by `SATURATION_SLACK`, stated in §2 rather
  than hidden in a tolerance. It is *not* a general safety net: it covers rounding at saturation
  only. A denormalized input or a large sector-leakage term would exceed the exact overlap by ~1e-5,
  nine orders past the slack — those are handled by the input guard and R2b, not by widening this.
- **R2 — self-mode provenance.** The floor is the caller's responsibility; the function does not
  know whether it was handed an oracle or a heuristic, and does not pretend to.
- **R2b — the sector restriction is inherited in NAME only.** `hf_overlap_certificate` justifies
  reading "spectrum" as "reachable spectrum" because *u*'s unreachable components vanish **by
  construction**. The chained bound applies Davis–Kahan to **v**, whose unreachable components vanish
  only *numerically*. Measured on the gated set: leakage ‖P_unreach v‖ runs 0 (stretched H₂) to
  **3.6e-14** (linear H₆, M=6), with square H₄ at 1.9e-14 — and 4–12 unreachable levels genuinely sit
  *below* E₁_reach, so the caveat is not vacuous. The induced error is O(leak²) = 3.6e-28 … 1.3e-27,
  ~14 orders below `SATURATION_SLACK`, so nothing here is violated.
  **But it scales as ε², and the margin is not structural — it is a property of this set.**
  `REACHABLE_TOL_CERTIFIED = 1e-10` admits amplitudes up to 1e-5, and the excluded a = 1.10 / 1.35
  square-H₄ geometries show what that costs: leakage **1.7e-6 / 1.6e-6**, eight orders above the clean
  set, giving a neglected term of ~3e-12 — **~200× above `SATURATION_SLACK`**, i.e. large enough to
  matter to G1 had those geometries been gated. R3's "no well-defined target" and this are therefore
  the same underlying issue seen twice, and the exclusion is load-bearing for R2b, not only for R3.
  G7 pins both halves; the excluded geometries serve as its positive control.
- **R3 — one clean set.** 5 systems; square H₄ is included only at a = 1.05, because a = 1.10 / 1.35
  have no well-defined target (`SPEC_symmetry_reachability`) — and, per R2b, are also where the
  neglected leakage term is largest. Two independent reasons to exclude them, one gated set.

## 9. Deliverables

- `hybrid_quantum_solver/certified_overlap/krylov_refine.py` — the implementation.
- `tests/test_chained_overlap_spec.py` — G1–G7.
