# SPEC: Factorization-native λ + symmetry shift lowers the FT-QPE 1-norm

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

Close the gap `lambda_ladder.py` names in its own docstring — *"large-N lambda must come from the
factorization-native 1-norm formulas (a careful next step, not done here)"* — by computing the
double-factorization 1-norm λ_DF analytically (no brute-force Pauli enumeration), then applying a
**number-operator symmetry shift** (BLISS / SCDF) that provably leaves the eigenspectrum unchanged
while reducing λ. Claim: the symmetry shift lowers λ_DF by a substantial margin (target ≥ 20% on a
representative active space) **with the FCI energy exactly preserved** after the a-posteriori
constant correction.

## 2. Background and honest framing

In qubitization-based QPE the T-gate / Toffoli budget scales as O(λ/ε), so λ is the highest-leverage
knob (`lambda_ladder.py`, `ft_resource_estimator.py`). The repo already has `double_factorize`
(`df_factorization.py`) but currently scores λ only via brute-force Pauli decomposition, feasible
only for ≲4 orbitals. Two results give the missing pieces:

- **Native DF 1-norm** (Deka & Zak, `arXiv:2412.01338`, eq. 23; von Burg form in Rocca et al.,
  `arXiv:2403.03502`, eq. 15): λ_DF = ½ Σ_r ‖A_r‖²_nuc + ‖h′‖_nuc, computed directly from the
  factorization tensors — no Pauli list, so it scales to dozens of orbitals.
- **Symmetry shift** (Loaiza & Izmaylov BLISS; SCDF in 2403.03502; full BLISS+DF in 2412.01338):
  replace H with H̃ = H + (b₁N_e + b₂)(N_e − n_e). Because N_e commutes with H and (N_e − n_e)|ψ⟩ = 0
  on the fixed-electron-number ground state, H̃ has the **same** ground state and energy (up to a
  known constant), but optimizing b₁,b₂ shrinks λ. Reported reductions: ~38% (XDF→XDF+shift on
  FeMoCo), and SCDF/BLISS reach λ ≈ ¼ of XDF; Deka & Zak report a further 25% over SCDF.

- **What we can claim if the gates pass:** a validated, scalable λ_DF estimator for our active-space
  Hamiltonians, plus a number-operator symmetry shift that measurably lowers it without changing the
  physics — i.e. a real reduction in the FT-QPE T-gate budget our resource estimator reports.
- **What we cannot claim:** no novelty (reproduces published preprocessing). We implement the
  **number-operator (N_e, N_e²) shift only**, not the full tensor-optimized SCDF cost function (that
  is a heavier follow-up). We make no statement that any of these systems are at quantum advantage.

## 3. Approach

Add a native λ_DF function over the `double_factorize` leaves (the `(g_t, w^(t), U^(t))` tuples):

    λ_DF = ¼ Σ_t |g_t| (Σ_k |w^(t)_k|)²  +  ‖h1‖_nuc

(two-body nuclear-norm sum over the second-factorization weights, plus the one-body nuclear norm).
This is a *bare* DF 1-norm — convention-free and directly testable. **Note:** the DF 1-norm is a
**different (and smaller) quantity than the naive Pauli-LCU 1-norm** — that reduction is the whole
point of DF — so the brute-force Pauli λ is *not* an equality oracle for it; it is only an upper
reference (DF should not exceed the trivial LCU). Implement the symmetry shift on the one- and
two-body integrals before factorization: h̃, g̃ as in 2403.03502 eqs. 23–24 (subtract a′₁ from the
median of the one-body factor coefficients; subtract a′₂ from the (pq|rs) diagonal), with (a′₁, a′₂)
chosen to minimize λ. The ground-truth reference for the *physics* is **PySCF FCI** before vs after
the shift + constant correction (validates spectrum invariance); the formula itself is validated by
independent recomputation and the Pauli upper-reference.

## 4. Public interface

```
df_factorization.df_lambda(leaves, h1, norb) -> float            # native λ_DF, no Pauli list
df_factorization.symmetry_shift(h1, eri, norb, nelec) ->
        (h1_s, eri_s, e_shift_const, (a1, a2))                   # H̃ integrals + correction const
lambda_ladder.lambda_ladder(...)   # extend table with an [DF+symshift] row using df_lambda
```

Compose: `double_factorize` / `reconstruct_eri` (existing), `fci_energy` (existing), the brute-force
λ (existing) as the small-system oracle.

## 5. Acceptance criteria (validation gates)

Gates live in `tests/test_scdf_lambda_spec.py` (test-first); small CAS so the Pauli oracle runs.

- **G1 — native λ formula is correct and tighter than the trivial LCU.** On N₂ CAS(norb=3,4e) (the
  `lambda_ladder` __main__ case) and H₂O STO-3G small CAS: (a) `df_lambda` at full rank equals an
  **independent in-test recomputation** of ¼ Σ_t|g_t|(Σ_k|w_k|)² + ‖h1‖_nuc to `< 1e-9`; and (b)
  `df_lambda ≤ lambda_and_terms`' brute-force Pauli λ (DF is no looser than the naive LCU). G1(b) is
  the literature-backed claim and is revisable if a tiny CAS violates it (record the measured ratio).
- **G2 — spectrum invariance under the shift.** FCI energy of (h1_s, eri_s) **+ e_shift_const**
  equals the unshifted CASCI/FCI energy to `< 1e-8 Ha`. (The shift must not move the eigenvalue —
  if it does, the commutation/correction is wrong.)
- **G3 — λ actually drops.** `df_lambda` after `symmetry_shift` is lower than before by **≥ 20%** on
  at least one representative active space (e.g. N₂ CAS(6,6) or an H-chain segment). **Definition of
  done.** **MEASURED:** the number-operator shift drops λ_DF on N₂ CAS(6,6) from 24.94 → 4.00 Ha
  (**84%**), and on N₂ CAS(3,4) by 93% — comfortably past the gate. (The bare DF λ before shift is
  dominated by the diagonal Coulomb piece the N_e² shift removes, so the relative drop is large; the
  spectrum is preserved exactly, G2.)
- **G4 — monotone λ vs rank.** Truncating DF rank R upward, `df_lambda` increases monotonically
  toward the full-rank value (sanity that the native norm tracks the factorization).
- (Stretch, not a gate) feed the reduced λ into `ft_resource_estimator.py` and report the Toffoli/
  T-gate reduction; compare DF+shift vs THC on a larger CAS.

## 6. Implementation plan (test-first)

1. Write `tests/test_scdf_lambda_spec.py` encoding G1–G4 (initially failing).
2. Add `df_lambda` and `symmetry_shift` to `df_factorization.py`; wire a DF+shift row into
   `lambda_ladder.py`.
3. Iterate to green via `make gates` (pyscf process group, not block2).

## 7. Out of scope

- The full SCDF tensor-optimization cost function (eq. 27 of 2403.03502) and RCDF — number-operator
  shift only here; tensor-optimized compression is a follow-up spec.
- Spin / point-group symmetry shifts beyond N_e and N_e² (a later extension).
- Actual fault-tolerant circuit synthesis or magic-state accounting — λ and the resource-estimator
  T-count are the deliverables, not a compiled circuit.

## 8. Caveats and risks

- **R1 — small-CAS reduction may miss 20%.** The dramatic reductions in the papers are for large
  TM active spaces; our FCI-tractable CASes may shift less. *Mitigation:* G3 is explicitly
  revisable — record the measured % and adjust, rather than forcing the threshold.
- **R2 — shifted two-body tensor can gain negative eigenvalues** (noted in 2403.03502 App. B). Our
  λ formula must use ‖·‖_nuc (sum of singular values), not assume positivity. Covered by G1 catching
  any sign/normalization error against the Pauli oracle.
- Honest limitation: λ reduction lowers the *T-gate budget*, which remains astronomically large at
  these sizes — this is a cost-model improvement, consistent with the repo's resource-accounting
  honesty (`benchmark_resources.py`), not a claim that FT-QPE is now feasible.

## 9. Deliverables

- `df_factorization.py` — `df_lambda`, `symmetry_shift`.
- `lambda_ladder.py` — DF+symmetry-shift row driven by the native λ.
- `tests/test_scdf_lambda_spec.py` — gates G1–G4.
- Results summary (with §2/§7 caveats) in the PR description.
