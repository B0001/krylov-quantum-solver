# SPEC: Interval-dominance screening — sound pruning, and the economics it actually buys

**Status:** v2.0 — IMPLEMENTED. G1 and G4 pass; **G2 FALSIFIED** (the ≥75% saving needs an unsound
bound); **G3 was vacuous** and is replaced by the ratio that governs pruning; G1b and G5 added.
v1.0 is preserved verbatim in §8 with what resolved each claim.
**Depends on:** `temple_bounds.py`, `hybrid_quantum_solver/quantum_krylov_solver.py`
**Gates:** `tests/test_screening_loop_spec.py` (G1, G1b, G2, G3, G4, G5 — 6 tests, 16 s) ·
`screening_loop.py` · `data/screening_loop_performance.csv`

---

## 1. Goal

Prune a candidate library using **certified** energy brackets: if candidate *j*'s certified lower
bound exceeds the best certified upper bound anywhere in the library, *j* provably cannot be the
global minimum and can be dropped before paying for convergence.

**Answers.** The dominance arithmetic is sound and the loop never loses the winner. But the saving
v1.0 advertised is not a property of the method, and the number it named is not reachable soundly:

- With a **rigorous** bracket the saving is **62.5–66.7%**, not ≥75%.
- **79.2%** is reachable — but only in Temple's **self** mode, which this spec measures to be
  **not a bound at all** (up to 0.695 mHa *above* the true energy). The headline was bought with
  an invalid certificate.
- How much any library saves is set by **one ratio**: its energy spread over the bracket width at
  the cheap dimension.

## 2. What v1.0 assumed, and what is actually true

**The certificate a screening loop actually has is the unsound one.** Temple's lower bound needs a
separator ε ≤ E₁. Oracle mode takes the exact E₁ — but knowing E₁ presupposes the candidate is
already solved, which is precisely what screening exists to avoid. Self mode estimates
ε = θ₁ − σ₁ from the same Krylov data, and **that estimate is not rigorous**: across the WIDE
library at M = 2…16 it exceeds the exact FCI energy on 5 (candidate, M) pairs, worst by
**0.695 mHa**, while oracle mode violates on **none** (G1b). A lower bound that can sit above the
true energy is licence to prune a candidate that is still viable.

This is the same failure `SPEC_excited_state_certification` established one rung up, and screening
**inherits** it: there, a self-certified separator cannot be trusted because a Krylov subspace
cannot certify an eigenvalue *count*; here the same unverifiable premise decides whether a molecule
is thrown away.

**Why no false prune is observed anyway, and why that is not a guarantee.** On a **homogeneous**
library — one molecule, swept geometry — every candidate's bound is inflated together, so the
winner's margin `L_winner − U_best` stays negative (measured ≤ −0.43 mHa) and the error is
common-mode. That is a measured property of same-molecule libraries, **not** a theorem. A
heterogeneous library whose members converge at different rates has no such protection, and this
spec does not claim safety there.

**STO-3G H₂ cannot test this spec.** Its Krylov space saturates at M = 2 (bracket width 0.00 mHa,
7/8 pruned immediately), so "cheap versus converged" does not exist and the loop measures nothing.
The library is H₄ chains, the smallest system here still unconverged at M = 2. This is the same
category error `SPEC_excited_state_certification`'s G2 records for the same molecule.

## 3. Approach (as built)

1. **Library.** Linear H₄ chains. **WIDE** = a ∈ [0.8, 2.2] Å (spread 285 mHa); **TIGHT** =
   a ∈ [1.00, 1.07] Å (spread 18.9 mHa). Two libraries, because one library cannot distinguish a
   property of the method from a property of the library.
2. **Loop.** Start at M = 2; refine every survivor; deactivate each active candidate with
   `L_j > U_best`; step M → M + 2 until one survivor or `max_m`. A **vacuous** (−inf) Temple bound
   never prunes — the conservative direction, gated in G5.
3. **Cost model, stated.** Cost = **count of Krylov basis vectors**, one controlled
   time-evolution family each, assuming the basis is extended incrementally so a candidate stopped
   at M costs M. Brute force = N · max_m. This is a count, **not** a wall-clock or gate count, and
   no QPU is involved.

## 4. Public interface (`screening_loop.py`)

```python
@dataclass
class Candidate:                  # name, mh, eps, active, current_m, lower, upper, cost, pruned_at
    def refine(self, m: int) -> None

interval_dominance_sweep(library, max_m=12, m_step=2, start_m=2) -> (Candidate, dict)
brute_force_sweep(library, target_m=12)                          -> (Candidate, dict)
sweep_metrics(library, max_m, pool_curve)                        -> dict
h4_library(spacings, oracle=False)                               -> list[Candidate]
```

**Deviations from v1.0, all deliberate:**

- `Candidate` carries **`eps`** (the Temple separator) and **`pruned_at`**. Without `eps` the
  sound and unsound modes cannot be compared, which is the whole finding; without `pruned_at` a
  false prune cannot be detected after the fact.
- `sweep_metrics` returns **`spread_over_width`** beside the saving. A saving quoted without the
  library that produced it is not a result (§6.1), so the interface refuses to report one alone.
- `cumulative_circuits` is renamed **`cost`** and documented as a basis-vector count rather than
  circuits, because the circuit count per basis vector is not modelled here.

## 5. Acceptance gates (as gated, with the measured numbers)

- **G1 — Pruning never loses the true winner (unchanged).** WIDE and TIGHT × oracle and self:
  **4/4** agreements with brute force, and the winner is never among the pruned. *Passes.*
- **G1b — The cheap certificate is not a certificate (NEW, the inherited kill).** Self mode:
  **5 violations**, worst **0.695 mHa** above exact. Oracle mode: **0**. Winner margin on the
  homogeneous library stays ≤ **−0.43 mHa**, so no false prune results *here*. *Passes.*
- **G2 — THE KILL (replaces v1.0's ≥75%).** Basis-vector cost, 8 candidates, max_m = 12 (brute =
  96). Sound/oracle: WIDE **36 = 62.5%**, TIGHT **32 = 66.7%** — both real, both under the bar.
  Unsound/self on WIDE: **79.2%**, the only schedule that clears it. Gated as a ceiling on the
  *sound* saving plus an explicit assertion that the unsound mode beats it. *Passes.*
- **G3 — Pruning power is the spread/width ratio (replaces a vacuous gate).** v1.0's "strictly
  non-increasing staircase" is true **by construction** — an active pool cannot grow — and gates
  nothing. At M = 2: WIDE spread/width = **4.45** prunes **6/8** at once; TIGHT spread/width =
  **1.00** prunes **0/8** and cannot move until M = 4. Same loop, same budget. The monotonicity
  claim is retained as a one-line side assertion rather than a headline. *Passes.*
- **G4 — Degenerate candidates resolve without a false prune.** Found while gating: twins 1e-5 Å
  apart (0.02 mHa) are **not** a tie — by M = 10 the oracle bracket is 0.0012 mHa wide, dominance
  is proven honestly, and the loop stops early. The gate therefore uses **exactly** degenerate
  geometries, where no bracket can separate the pair: both reach max_m, neither prunes, the other
  three do. *Passes.*
- **G5 — A vacuous bound never prunes (NEW).** Temple goes −inf for 2–3 of 8 candidates at
  M = 6–8 (ε ≤ θ₀). −inf never exceeds a finite upper bound, so those candidates survive. This is
  why pruning power is **non-monotonic in M** even though the pool itself never grows. *Passes.*

## 6. Findings

1. **The saving is a property of the library, not the method.** spread/width at the cheap
   dimension decides everything: 4.45 → 6/8 pruned at M = 2; 1.00 → 0/8. Any percentage must be
   quoted with that ratio, and `sweep_metrics` returns it so it cannot be dropped.
2. **The economics and the soundness are in direct tension.** The ≥75% target is met only by the
   bound that is not a bound. Rigorous screening costs ~62–67% savings instead — still worth
   having, and honestly smaller.
3. **Screening inherits the self-certification kill.** An oracle separator presupposes a solved
   candidate. This is not a gap to be closed by better estimation: per
   `SPEC_excited_state_certification`, a Krylov subspace cannot certify the eigenvalue count the
   separator needs.
4. **Common-mode error is why homogeneous libraries survive an unsound bound.** Every member is
   inflated together. This makes same-molecule sweeps far safer than the arithmetic alone implies,
   and it is precisely why a heterogeneous library is the dangerous case.
5. **Near-degeneracy is resolved by converging, not by tie-breaking.** Only an exact tie exhausts
   the budget; anything separable is separated as soon as the bracket narrows past the gap.

## 7. Honest scope and caveats

- **Exact statevector.** No shot noise, no device noise, no Trotter error. Shot noise would widen
  every bracket and push pruning to larger M — the direction that makes the saving smaller.
- **Cost is a basis-vector count**, not gates, shots or wall-clock, and assumes incremental basis
  reuse. No QPU is involved anywhere in this spec.
- **Two libraries, one molecule family, 8 candidates.** The spread/width *law* should generalise;
  the specific percentages will not.
- **Soundness is gated only for homogeneous libraries.** No claim is made for heterogeneous ones,
  where common-mode protection does not apply.
- **Ground states only** — no binding energies, no free energies, nothing about actual drug
  screening beyond the shared arithmetic. v1.0's "materials and drug discovery" framing is
  withdrawn as scope inflation.
- The dominance test inherits `temple_bounds`'s own scope, including the HF-reachable sector
  restriction.

## 8. v1.0, and what resolved each claim

> **G2 (v1.0):** *"For a library of 8 molecules, the adaptive loop must achieve a total
> circuit/shot budget reduction of ≥ 75% compared to brute force, demonstrating massive economic
> acceleration."*
> **G3 (v1.0):** *"Plotting the size of the active pool vs Krylov dimension M must show a strictly
> non-increasing staircase curve."*
> **§1 (v1.0):** *"…save over 80% of the total quantum compute budget (gates, circuits, and
> shots)."*

Measured: sound saving **62.5%** (WIDE) and **66.7%** (TIGHT); **79.2%** only in the mode G1b shows
is unsound; no gates, no shots and no QPU were modelled, so "80% of the quantum compute budget" is
withdrawn in favour of a stated basis-vector count. v1.0's G3 is true by construction and gates
nothing; it is demoted to a side assertion and replaced by the spread/width ratio. G1 and G4
survive, G4 with its premise corrected. The spec was revised rather than the tolerances loosened —
`specs/README.md` step 5.
