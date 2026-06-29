# Spec backlog — falsifiable hypotheses

One line each: **a claim + the cheap check that would prove or kill it.** When you pick one up,
copy `SPEC_TEMPLATE.md` to `SPEC_<slug>.md` and turn the check into gates. Keep the bar:
*if you cannot write the check, it is not a candidate.* Reproductions are fine and useful — just
label them as such; novelty claims need a reference that could embarrass them.

Status key: `[ ]` open · `[~]` specced · `[x]` done (link the spec) · `[-]` killed (record why).

## Open

- [~] **Hₙ to larger n, done right** — *Claim:* the TDL fit tightens at large n **with adequate D**
  (≈ 400/800/1600) **or** a bulk per-site estimator `(E(n) − E(n−Δ))/Δ`. *Gate:* leave-one-out
  shift < 0.1 mHa/atom. (Supersedes the cheap version below.) → [`SPEC_hchain_largen2.md`](SPEC_hchain_largen2.md)
- [ ] **Be₂ toward experiment** — *Claim:* core-valence correlation + a cc-pVxZ→CBS extrapolation
  moves the FCI/DMRG well depth from ~305 cm⁻¹ toward the experimental 929.7. *Gate (reproduction):*
  `|D_e − 930| < 100 cm⁻¹` at CBS+CV; `R_e` within 0.1 Å of 2.45. (Honest: reproduces a settled
  result.)
- [ ] **1D Hubbard vs Bethe ansatz** — *Claim:* the solver/DMRG reproduce the exact 1D Hubbard
  ground-state energy. *Gate:* `|E − E_Bethe| < 1 mHa/site` at half-filling, L=10, several U/t.
  (A clean analytic reference — strong falsifiability.)
- [ ] **GPU backend on real hardware** — *Claim:* `device="gpu"` reproduces CPU energies and reaches
  larger qubit counts. *Gate (on an NVIDIA node):* `|E_gpu − E_cpu| < 1e-6 Ha` at 16q; completes a
  28q statevector run. (Code exists, CPU-validated; needs a GPU node — see `tests/test_gpu_backend.py`.)
- [ ] **DMRG-referenced transition-metal active space** — *Claim:* for a low-spin NbN (or a smaller
  TM dimer) CAS large enough that FCI is intractable, DMRG gives a converged correlation energy.
  *Gate:* bond-dim extrapolation stderr < 1 mHa **and** agreement between two independent DMRG
  sweep schedules < 1 mHa. (Honest: a reference number, not a materials claim — finite cluster.)
- [~] **SKQD reproduces the exact-Krylov floor** — *Claim:* a sample-based Krylov path (diagonalize
  in determinants sampled from e^(−ikΔtH)|HF⟩) converges to FCI from above on H₄ / N₂ CAS(6,6).
  *Gate:* `E_skqd ≥ E_fci − 1e-6` always, and `|E_skqd − E_fci| < 1.6 mHa` at depth ≥ 6, high shots.
  (Reproduction of SKQD/SqDRIFT — `arXiv:2501.09702`, `2508.02578`; bridges the Krylov and SQD code
  paths. No advantage at this scale.) → [`SPEC_skqd.md`](SPEC_skqd.md)
- [~] **Factorization-native λ + symmetry shift** — *Claim:* a native double-factorization 1-norm
  (no brute-force Pauli, so it scales past ~4 orbitals) plus a number-operator symmetry shift lowers
  λ_DF by ≥ 20% with the FCI energy exactly preserved. *Gate:* native λ matches the brute-force Pauli
  oracle < 1e-6 on a small CAS; FCI invariant < 1e-8 Ha under the shift; λ drops ≥ 20%. (Reproduction
  of BLISS/SCDF — `arXiv:2403.03502`, `2412.01338`; closes the gap named in `lambda_ladder.py`.)
  → [`SPEC_scdf_lambda.md`](SPEC_scdf_lambda.md)

## Done

- [x] **single-ramp DMRG extrapolation** — → [`SPEC_singleramp.md`](SPEC_singleramp.md). One ramping
  run via block2 `get_dmrg_results()` agrees with the per-D protocol < 0.1 mHa, lands at FCI, and
  uses half the sweeps (gates G1–G3 in `tests/test_singleramp_spec.py`). `protocol="ramp"`.
- [x] **Hₙ thermodynamic limit** — bond-dim + n→∞ extrapolation. → [`SPEC_hchain_tdl.md`](SPEC_hchain_tdl.md)
  (e_∞ = −0.539967 ± 0.000107 Ha/atom; gates G1–G5 in `tests/test_hchain_tdl_spec.py`).

## Killed

- [-] **Hₙ to larger n, *cheaply*** (ramp + D=100/200/400) — → [`SPEC_hchain_largen.md`](SPEC_hchain_largen.md).
  *Killed:* D=400 truncates too hard as chain entanglement grows — stderr balloons to ~5 mHa, the
  discarded-weight extrapolation falls back to `invD` by n=30, and leave-one-out = 1.07 mHa/atom
  (gate < 0.1). The ramp protocol is fine; cheap bond dims are not. Superseded by "done right" above.
