# Spec backlog — falsifiable hypotheses

One line each: **a claim + the cheap check that would prove or kill it.** When you pick one up,
copy `SPEC_TEMPLATE.md` to `SPEC_<slug>.md` and turn the check into gates. Keep the bar:
*if you cannot write the check, it is not a candidate.* Reproductions are fine and useful — just
label them as such; novelty claims need a reference that could embarrass them.

Status key: `[ ]` open · `[~]` specced · `[x]` done (link the spec) · `[-]` killed (record why).

## Open

- [ ] **single-ramp DMRG extrapolation** — *Claim:* reading block2 `get_dmrg_results()` from ONE
  ramping run gives the same `E(D→∞)` as the per-D protocol at ≳3× less wall-time. *Gate:*
  `|E_singleramp − E_perD| < 0.1 mHa` on H₁₂ **and** wall-time < 0.4×. (Unblocks large-n Hₙ.)
- [ ] **Hₙ to larger n** — *Claim:* extending the TDL fit to n=40,60 tightens `e_∞`. *Gate:*
  leave-one-out shift < 0.1 mHa/atom with n up to 60 (needs the single-ramp speedup first).
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

## Done

- [x] **Hₙ thermodynamic limit** — bond-dim + n→∞ extrapolation. → [`SPEC_hchain_tdl.md`](SPEC_hchain_tdl.md)
  (e_∞ = −0.539967 ± 0.000107 Ha/atom; gates G1–G5 in `tests/test_hchain_tdl_spec.py`).

## Killed
<!-- record dead ideas + the one-line reason, so they stay dead -->
