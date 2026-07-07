# SPEC: Cost advisor — a per-molecule near-term-vs-FT verdict with a parametrized crossover

**Status:** IMPLEMENTED — gates G1–G4 green.

---

## 1. Goal

Turn the precision-cost *scaling law* (`precision_cost`) into an *engineering artifact*: for a target
accuracy ε on a given Hamiltonian, return a **verdict** — near-term certified Krylov vs FT-QPE, which
is cheaper — with the crossover ε* and the verdict's robustness parametrized by the one honest
unknown, the per-query-to-per-shot cost ratio ρ. Falsifiable: a verdict that does not flip at the
computed ε*, or a robustness claim that a decade of ρ overturns, breaks it.

## 2. Background and honest framing

The bridge gives two cost curves (N(ε) ∝ (λ_meas/ε)², Q(ε) ∝ λ_DF/ε). A recommendation needs a
common cost model, whose only free parameter is ρ = cost_per_query / cost_per_shot (a shot is one
measurement circuit; a query is one qubitization/T-heavy block encoding — ρ ≫ 1 in practice, but its
value is hardware-specific and genuinely unknown). This spec is the repo's `validate_and_cost` /
`cross_check` instinct applied to the certified-vs-FT choice.

- **What we can claim:** a deterministic per-molecule verdict at (ε, ρ), the crossover ε*(ρ), and a
  cross-validated *robust* recommendation over a ρ-range (FT/near-term wins regardless of the
  unknown, or the ρ* where it flips).
- **What we cannot claim:** the true ρ (hardware-specific, not computed); absolute wall-clock or
  dollar cost. The verdict is a resource-count decision under a stated cost model; the exponent gap
  is the robust part, the location moves with ρ.

## 3. Approach

Reuse `precision_cost.{near_term_shots, ft_queries, crossover_epsilon}` (and the SCDF-shifted λ_DF).
Common cost: near-term = N(ε)·1, FT = Q(ε)·ρ. Verdict = cheaper of the two; ε*(ρ) = the precision at
which they cross (below ε* → FT). Robustness: sweep ρ over decades and report whether the verdict at
ε is invariant. **Reference:** the crossover algebra is exact (ε* = 2z²λ_meas²/(πρλ_DF)); gates check
the verdict flips exactly there and the robustness threshold ρ* = N(ε)/Q(ε).

## 4. Public interface

```
cost_advisor.Verdict(eps, rho, lam_meas, lam_df, shots, queries, eps_star, cheaper, rho_flip)
cost_advisor.advise(lam_meas, lam_df, eps, rho=1.0, z=2.0) -> Verdict
cost_advisor.advise_from_integrals(h1, eri, e_core, nelec, norb, eps, rho=1.0, shift=True) -> Verdict
cost_advisor.robust_over_rho(lam_meas, lam_df, eps, rho_lo, rho_hi, z=2.0) -> (verdict:str, robust:bool)
cost_advisor (CLI)                                            # verdict table for a few molecules
```

## 5. Acceptance criteria (validation gates)

`tests/test_cost_advisor_spec.py` (test-first).

- **G1 — the verdict flips at ε*.** At fixed ρ, `advise` returns FT for ε just below ε* and near-term
  just above; `cheaper` is consistent with the direct cost comparison at every ε tested.
- **G2 — ε* is parametrized by ρ correctly (DEFINITION OF DONE).** ε*(ρ) ∝ 1/ρ (a decade more
  per-query cost shrinks ε* tenfold), and the verdict-flip ρ* at a fixed ε equals N(ε)/Q(ε).
- **G3 — cross-validated robustness.** At chemical accuracy (1.6 mHa) the N₂ flip-ρ is ~2×10⁵, so FT
  is robustly cheaper across ρ ∈ [1, 10⁵]; `robust_over_rho` returns FT/robust below the flip and
  "mixed" for a range straddling it.
- **G4 — integrals path + honesty.** `advise_from_integrals` (with the SCDF shift) agrees with
  `advise` on the same λ's; the verdict object exposes ρ so the unknown is never hidden.

## 6. Implementation plan (test-first)

1. `tests/test_cost_advisor_spec.py` encoding G1–G4 (initially failing — no module).
2. `cost_advisor.py` composing `precision_cost` into a verdict + robustness sweep.
3. `make gates`.

## 7. Out of scope

- Computing the true ρ (hardware-specific) or absolute cost.
- Depth/coherence constraints, error-correction overhead.
- Systems beyond statevector validation of the near-term side (the ladder spec handles scale).

## 8. Caveats and risks

- **R1 — ρ is the honest unknown.** Mitigated by making it an explicit input and reporting
  ρ-robustness, never a single hidden number.

## 9. Deliverables

- `cost_advisor.py`, `tests/test_cost_advisor_spec.py`, this spec, BACKLOG entry.
