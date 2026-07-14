# SPEC: ADAPT-VQE — does gradient-greedy operator selection actually buy a more compact ansatz?

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

---

## 1. Goal

`adapt_vqe.py` implements ADAPT-VQE (exact statevector): grow the ansatz one operator at a time,
each round picking the pool operator with the largest energy gradient, then re-optimize all
parameters. Its `__main__` block already asserts the result is variational (`>= CASCI`), but that
assertion has never been a CI gate, and the method's actual selling point — that GREEDY gradient
selection reaches a target accuracy with a MORE COMPACT ansatz than an arbitrary fixed operator
order — has never been checked at all. This spec gates both: the existing informal variational
check, and the new, sharper, comparative claim. False if a fixed/random operator order matches or
beats greedy on operator count to chemical accuracy on a system with real electron correlation.

## 2. Background and honest framing

- `adapt_vqe.py` already reuses validated primitives (`qubitization_blueprint`'s exact JW operators)
  and is exact-statevector — this spec adds no new physics, only a falsifiable gate around
  behavior the module already exhibits, closing a gap the repo's own culture flags ("a capability
  without a test that could have killed it is not done").
- **What you can claim if the gates pass:** ADAPT-VQE's ansatz stays variational at every growth
  step (never below CASCI) and converges to chemical accuracy within budget on three systems; and
  on a system with real multi-orbital correlation (H4 CAS(4,4)), greedy operator selection reaches
  chemical accuracy in markedly fewer operators than a fixed order does within a matched budget —
  the mechanism ADAPT-VQE is named for is real here, not just asserted.
- **What you cannot claim:** an advantage on every system — a small, weakly-correlated active space
  (LiH CAS(2,2)) shows NO advantage at all (both greedy and every tested random order need exactly
  1 operator), the recorded boundary; this is exact statevector (no hardware gradient-measurement
  cost, no shot noise); the pool is generalized singles+doubles only (no orbital-optimization,
  no UCC constraints); "compact" means operator count, not circuit depth or gate count.
- **Reference:** CASCI (dense diagonalization in the active space, via PySCF `mcscf.CASCI`) — the
  same reference `adapt_vqe.py`'s own `__main__` already uses.

## 3. Approach

Reuse `adapt_vqe.py`'s existing `adapt_vqe`/`build_pool`/`hf_state`/`adapt_ground_state` unchanged.
Add one small comparison baseline, `fixed_order_vqe`, which grows the SAME pool in a CALLER-GIVEN
fixed order (not by gradient) — otherwise byte-identical growth/re-optimize loop (BFGS, `gtol=1e-7`,
warm-started parameters) — so the only variable between greedy and baseline is the selection rule.
For each system: run greedy once (deterministic — no randomness in gradient argmax); run the fixed
baseline over several seeded random permutations of the pool, and record the operator count each
needs to first reach chemical accuracy (1.6 mHa) relative to CASCI, within a shared operator budget.

## 4. Public interface

```
adapt_vqe.fixed_order_vqe(H, hf, pool, order, max_ops) -> history
    # history: list of (n_ops, electronic_energy), same shape as adapt_vqe()'s history entries
```

## 5. Acceptance criteria (validation gates)

- **G1 — variational floor at every growth step (pins the existing informal assertion).** On H2
  CAS(2,2), LiH CAS(2,2), H4 CAS(4,4): every intermediate electronic energy in `adapt_vqe`'s
  history is `>= CASCI_electronic - 1e-6` Ha, not just the final one.
- **G2 — convergence to chemical accuracy within budget.** Greedy ADAPT-VQE reaches CASCI to within
  1.6 mHa on all three systems within `max_ops=30`. *Measured: reached in 1 / 1 / 9 operators
  (H2 / LiH / H4).*
- **G3 — THE FINDING (definition of done): greedy beats a fixed order on the correlated system.**
  On H4 CAS(4,4), greedy reaches chemical accuracy in 9 operators; a fixed fixed random order
  (5 seeds, budget = 9+5=14 ops) reaches it in NONE of the 5 seeds — a decisive, not marginal,
  advantage. *Measured: greedy=9; random (5 seeds) = [None, None, None, None, None] within 14 ops.*
- **G4 — boundary, recorded not smoothed over: NO advantage on a trivially simple active space.**
  On LiH CAS(2,2), greedy and every one of 5 random-order seeds reach chemical accuracy in EXACTLY
  1 operator — adaptivity buys nothing when the pool's first operator (in any order) already
  resolves the correlation. *Measured: greedy=1; random (5 seeds) = [1, 1, 1, 1, 1].*

> Definition of done: **G3**. If a future correlated system shows fixed-order matching or beating
> greedy, that contradicts ADAPT-VQE's core premise here and must be recorded, not silently dropped.

## 6. Implementation plan (test-first)

1. Write `tests/test_adapt_vqe_compactness_spec.py` encoding G1-G4 (RED — `fixed_order_vqe` doesn't
   exist yet).
2. Add `fixed_order_vqe` to `adapt_vqe.py` (~15 lines, mirrors `adapt_vqe`'s existing growth loop).
3. `make gates` / targeted pytest to green; ruff clean.

## 7. Out of scope

- Circuit depth / T-gate cost of the compact ansatz (operator count only, per §2).
- Orbital optimization, UCC-constrained pools, or larger active spaces (CAS(4,4) is already ~400
  pool operators; a bigger system is a follow-up, not this spec).
- Hardware gradient measurement (shot noise on the ADAPT gradient step) — a separate spec.

## 8. Caveats and risks

- **R1 — the specific operator counts are system- and pool-dependent** (H2/H4 pool sizes 21/406);
  a different active space could show a smaller or larger greedy-vs-random gap. The falsifiable
  claim (G3) is that a decisive gap EXISTS on at least one correlated system tested, not a universal
  compactness ratio.
- **R2 — BFGS is a local optimizer**, not guaranteed global; a fixed-order run could in principle
  get stuck and understate its own best-case performance. Mitigation: same optimizer settings,
  same warm-start convention, for both greedy and baseline — an apples-to-apples comparison, not a
  claim about either optimizer's absolute quality.
- Honest limitation: 3 systems, one pool type, exact statevector.

## 9. Deliverables

- `adapt_vqe.py` — add `fixed_order_vqe` (baseline comparison).
- `tests/test_adapt_vqe_compactness_spec.py` — gates G1-G4.
- Results summary (with R1/R2 caveats) in the PR description / BACKLOG entry.
