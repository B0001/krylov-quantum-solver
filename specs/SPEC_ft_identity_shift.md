# SPEC: The identity term inflates the FT (qubitization) λ too — the third strike

**Status:** IMPLEMENTED (gates G1–G5 green in `tests/test_ft_identity_shift_spec.py`).

---

## 1. Goal

`SPEC_shift_both_sides` and `SPEC_lambda_meas_identity` caught the identity Pauli term
inflating the **near-term** measurement 1-norm λ_meas, twice, independently, in this repo's own
history. The **FT/qubitization** side was never touched: `build_walk_operator` loads every term
from `pauli_decompose` into PREPARE, identity included, so `lambda_1norm` carries dead constant
mass — the identity term has zero variance and needs no ancilla amplitude at all, so it is free
to drop and add back classically. Claim: the identity fraction is material (not a rounding
effect) on real active-space Hamiltonians, the drop is free and exact (spectrum-preserving), and
it earns a genuine walk-step/T-gate *budget* reduction at fixed target precision. Falsifiable: if
the identity fraction is below ~25%, or if recentering does not preserve the spectrum exactly, or
if the reduced λ does not translate into a smaller walk-step budget at a fixed target ε, the claim
breaks.

## 2. Background and honest framing

- **What we can claim if the gates pass:** the FT walk-operator's PREPARE 1-norm λ can be
  losslessly reduced by exactly `|c_I|` (the qubit Hamiltonian's identity coefficient) for every
  active-space system tested, with the ground state and full spectrum preserved exactly via a
  classical additive correction; this reduction shrinks the walk-step budget `Q(ε) ~ O(λ/ε)` by
  the same ratio `λ/λ'` at any fixed target precision ε.
- **What we cannot claim:** that this changes exponents (it is a constant, not a scaling
  improvement — same caveat as `SPEC_shift_both_sides` R2); that the *point-estimate* error at a
  single literal fixed phase-bit count `t` scales predictably with any simple closed form
  (§3 KILLED FINDING — it does not, and this spec records why); or that spectral centering is
  novel — it is standard qubitization practice, reproduced here as an audit of this repo's own
  path (matches the bead's own caveat).

## 3. THE FINDINGS

### 3a. The identity fraction is real and large (reconfirms the scout probe)

Recomputed directly from `qubitization_blueprint.pauli_decompose` on CASCI active spaces
(sto-3g): `|c_I|/λ` = **30.1% (H₂ CAS(2,2)), 44.2% (LiH CAS(2,2)), 48.5% (N₂ CAS(3,4)), 45.6%
(H₂O CAS(3,4))**. These match the backlog scout probe (30.1/44.2/49.3/45.6%) to within the
expected geometry/rounding noise — large fractions, not an artifact. (G1)

### 3b. The drop is free and exact

For every system, `λ' = λ − |c_I|` to floating-point precision, and re-encoding the identity-free
terms into a fresh walk operator `W'` reproduces the **exact same** Hamiltonian spectrum once
`c_I` is added back classically: `λ'·cos(θ'_k) + c_I == E_k` to < 1e-8 Ha for every eigenvalue,
on both `H₂` and `LiH` (G2). This is the FT-side analogue of `shift_both_sides.shot_lambda`.

### 3c. KILLED — the naive "fixed-t point-estimate ratio" does not track any simple prediction

The bead's fix concept proposed `λ_eff = sqrt(λ² − E₀²)` (the *local* slope `λ·sin θ₀` of
`E = λ cos θ` at the true ground eigenvalue) as "the honest driver," predicting that the realized
QPE point-estimate error at a fixed phase-bit depth `t` should shrink by the ratio
`λ'_eff/λ_eff` when recentered. **Measured and killed:** sweeping `t = 4…22` on all four systems,
the ratio of recentered to raw point-estimate error (via `qpe_walk_readout.run_qpe`'s argmax
decode) ranges from **0.003 to 5.8** at *adjacent* `t` values, never settling near either the
`λ_eff`-ratio prediction (0.64–0.80) or the plain `λ`-ratio (0.51–0.70) — even at `t` as deep as
22 bits, where the point estimate is already at the 1e-6 Ha level. This is the same staircase/
dyadic-rounding phenomenon the (separate, not-yet-implemented) BACKLOG entry on the `π·sin θ₀`
envelope bound warns about: the argmax decode's error is `λ·sin θ₀·δ(t)` where `δ(t)` is the
*rounding offset* of the true phase to the nearest dyadic grid point — a quantity that does not
vary smoothly in `t`, so a literal single-`t` comparison is not a meaningful test of any
1-norm-based prediction, "predicted" or otherwise. **Revision:** the falsifier is restated at the
*envelope* level (G3, below), where it is well-defined and passes.

### 3d. The envelope/budget claim survives

Over a wide `t`-sweep (`t = 4…20`), the worst-case constant `max_t[err(t)·2^t/λ]` stays in a
narrow band (2.2–3.0) for the **raw** encoding and (2.7–3.0) for the **recentered** one, on every
system — i.e. the *same* empirical bound style already validated in
`SPEC_qpe_readout_laws` G2 (`err(t) ≤ 3·λ/2^t`) extends to the identity-free walk operator with
the same constant, just with the smaller λ'. Algebraically, that means the number of phase bits
(walk steps `Q ~ 2^t`) needed to *guarantee* a target precision ε via that bound drops by exactly
`λ/λ'` — a real, exploitable, and non-noisy T-gate/walk-step budget reduction, confirmed
end-to-end by actually running `run_qpe` at the computed `t*` for both variants and checking both
land under ε (G3, the definition of done).

### 3e. `df_lambda` vs identity-excluded naive λ — NOT-COMPARABLE-AND-VACUOUS

The scout also found `df_factorization.df_lambda` **losing** to the identity-*excluded* naive
Pauli λ on all four systems (df_lambda 2.60/1.54/7.58/8.67 vs excluded-λ 1.89/0.89/4.28/5.40),
which superficially looks like a flip of `SPEC_scdf_lambda` G1(b) (`df_lambda ≤` naive Pauli λ).
**It is not a flip — it is a category error, confirmed two ways (G5):**

1. `SPEC_scdf_lambda` G1(b) compares `df_lambda` against the **identity-included** naive λ
   (`lambda_ladder.lambda_and_terms`), and that comparison still holds on all four systems
   (`df_lambda ≤ λ_incl` — unchanged, unaffected by anything in this spec).
2. `df_lambda(leaves, h1, norb)` is a **pure function of `(h1, eri)`** — the double-factorization
   tensors — and never touches a Pauli-identity coefficient at all; there is no notion of
   "`df_lambda` with the identity removed" because DF's rotated-number-operator block encoding
   never separately counted one. Confirmed directly: `df_lambda` is bit-for-bit identical whether
   or not the corresponding qubit Hamiltonian's Pauli-identity term is tracked or dropped, because
   `(h1, eri)` do not change under that operation. Comparing it against a quantity
   (identity-*excluded* naive λ) that is only defined *because* one specific circuit construction
   (brute-force Pauli LCU) happens to expose a removable slack the other construction never
   counted the same way is exactly the "two λ's may not be comparable" caveat the bead flagged —
   confirmed, not merely suspected. **G1(b) is not flipped; the identity-excluded comparison is
   vacuous.**

## 4. Public interface

```
qubitization_blueprint.split_identity(terms) -> (non_identity_terms, c_identity)
qpe_walk_readout.run_qpe(h1, eri, norb, e_core, trial_vec, t_bits, recenter=False) -> tuple
    # recenter=True: PREPARE skips the identity load; lambda is the reduced, identity-free
    # 1-norm; the returned energy is unaffected (c_identity added back classically).
```

## 5. Acceptance criteria (validation gates)

`tests/test_ft_identity_shift_spec.py`.

- **G1 — identity fraction reconfirmed.** `|c_I|/λ > 25%` on all four systems (H₂, LiH,
  N₂ CAS(3,4), H₂O CAS(3,4)); matches the scout probe order of magnitude.
- **G2 — the drop is free and exact.** `λ' == λ − |c_I|` to `< 1e-9`; re-encoding the
  identity-free terms into a fresh walk operator and adding `c_I` back reproduces every original
  Hamiltonian eigenvalue to `< 1e-8` Ha (H₂, LiH — small enough for the exact-matrix build).
- **G3 — the budget claim (DEFINITION OF DONE).** For a target ε, the phase-bit count
  `t* = ceil(log2(3λ/ε))` from the validated bound style of `SPEC_qpe_readout_laws` G2 gives a
  walk-step ratio `2^(t*_raw) / 2^(t*_recentered) == λ/λ'` to machine precision (exact algebra);
  running `run_qpe` end-to-end at each variant's own `t*` on all four systems confirms the
  realized error lands under ε for **both** variants (validating the bound is not just algebra
  but an honored guarantee), and the recentered variant needs a strictly smaller (or equal)
  walk-step budget on every system. Dies if the recentered run ever needs *more* steps or misses
  its own bound.
- **G4 — the point-estimate ratio is killed at the literal fixed-`t` level (the recorded
  revision).** Sweeping `t = 4…20`, the ratio of recentered/raw point-estimate error is NOT
  confined near either the `λ`-ratio or the `sqrt(λ²−E₀²)`-ratio prediction — it must range over
  at least a factor of 5 across the sweep on every system, evidencing the staircase/rounding
  noise that makes the naive fixed-`t` framing unfalsifiable-by-simple-formula. This gate goes
  red (i.e. the finding would be overturned) if some future numerical fix makes the ratio
  suddenly well-behaved — that would be worth knowing.
- **G5 — the df_lambda comparability question (explicitly resolved).** (a) `df_lambda ≤`
  identity-included naive λ still holds on all four systems (SPEC_scdf_lambda G1(b) unbroken);
  (b) `df_lambda` is bit-identical regardless of whether the qubit-operator identity term is
  tracked or dropped (it is a function of `(h1, eri)` alone); (c) `df_lambda >`
  identity-excluded naive λ on all four systems (reconfirms the scout), but per (a)+(b) this is
  not a valid apples-to-apples comparison — recorded as NOT-COMPARABLE-AND-VACUOUS, not a flip.

## 6. Out of scope

- The (separate, not-yet-implemented) BACKLOG entry deriving a tight `π·sin θ₀` envelope bound
  for the Trotter/Krylov staircase error — this spec only *uses* the qualitative phenomenon
  (§3c) as the reason the naive fixed-`t` framing fails; deriving/gating the tight constant is
  that other entry's job.
- Applying the identity-drop to `qubitization_blueprint.build_walk_operator`'s ancilla-count
  accounting in `ft_resource_estimator.py` / ancilla-index bookkeeping beyond `λ` itself (the
  T-gate cost model's identity-loading cost, if any, is not itemized here — λ is the deliverable,
  matching `SPEC_scdf_lambda`'s own scope note).
- Re-deriving the SCDF/BLISS number-operator shift on the FT side (already done,
  `SPEC_scdf_lambda`) — this spec is about the *separate*, purely-additive identity-Pauli term,
  which is present whether or not the number-operator shift has also been applied.

## 7. Caveats and risks

- **R1 — reproduction, not novelty.** Spectral centering/identity-dropping in a block encoding is
  standard qubitization practice; this is an audit of this repo's own path, not a new technique
  (matches the bead's own caveat).
- **R2 — constant, not exponent.** Same as `SPEC_shift_both_sides` R2: this moves the walk-step
  *budget* by a constant factor `λ/λ'`, never the `O(λ/ε)` scaling law itself.
- **R3 — small-CAS-only exact verification.** G2's exact spectral round-trip is checked on H₂/LiH
  (≤4 qubits) where the brute-force Pauli/exact-matrix construction is tractable; larger CAS
  systems (N₂, H₂O CAS(3,4), 6 qubits) are checked on G1/G3/G5 (which do not require the exact
  `2^n × 2^n` matrix build at every step) but not on G2's exact-eigenvalue round-trip.

## 8. Deliverables

- `qubitization_blueprint.py` — `split_identity`.
- `qpe_walk_readout.py` — `run_qpe(..., recenter=False)`.
- `tests/test_ft_identity_shift_spec.py` — gates G1–G5.
- This spec; `specs/BACKLOG.md` close-out.
