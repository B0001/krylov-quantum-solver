# Spec backlog — falsifiable hypotheses

One line each: **a claim + the cheap check that would prove or kill it.** When you pick one up,
copy `SPEC_TEMPLATE.md` to `SPEC_<slug>.md` and turn the check into gates. Keep the bar:
*if you cannot write the check, it is not a candidate.* Reproductions are fine and useful — just
label them as such; novelty claims need a reference that could embarrass them.

Status key: `[ ]` open · `[~]` specced · `[x]` done (link the spec) · `[-]` killed (record why).

## Open

- [x] **Nb₃X₈ vs magnetometry — a parameter-free prediction against the lab** *(the profitable turn:
  internal cross-check → falsifiable prediction against experiment)* — the ab-initio interlayer
  singlet–triplet gap J (from the same downfolded bilayer as [`SPEC_odmd_spin.md`](SPEC_odmd_spin.md),
  cRPA params, **no fit**) predicts the *scale and ordering* of the **measured** magnetic-singlet
  transitions of Nb₃Cl₈ (~90 K, Sheckelton 2017 / `arXiv:1701.05528`) and Nb₃Br₈ (~382 K, Haraguchi
  2017): the exact dimer χ(T) peaks at k_BT_max ≈ 0.625 J → 479 K (Cl) < 862 K (Br), reproducing the
  observed Cl<Br ordering within an order of magnitude (G2). **THE FINDING (G3, DoD):** the isolated
  dimer **overpredicts Tc by 5.3× (Cl) / 2.3× (Br)** — an overcoupling that *weakens monotonically
  down the series*, the isolated-cluster→cooperative-lattice renormalization; and it exposes **two
  distinct couplings** — −J/4 = −192 K overshoots the measured Curie–Weiss θ_W = −13.1 K by 15×, so
  the interlayer J sets Tc while a separate weak *in-plane* exchange (absent from the bilayer dimer)
  sets θ_W. Numbers the cluster papers never reported. **Boundary (G4):** Nb₃I₈ has no
  interlayer-singlet transition (moment-retaining ground state) — excluded; the predictor sets
  *scales*, not a first-order cooperative transition. 100% primitive reuse (`susceptibility`,
  `dimer_exchange_analytic`) + a cited experimental table. Gates G1–G4 in
  `tests/test_nb3x8_magnetometry_spec.py`; `nb3x8_magnetometry.py`.
  → [`SPEC_nb3x8_magnetometry.md`](SPEC_nb3x8_magnetometry.md) (comparison vs measured references,
  not a fit; isolated bilayer dimer, density-density only).
- [x] **Be₂ toward experiment** — *Claim:* core-valence correlation + a cc-pVXZ→CBS extrapolation
  moves the well depth from ~305 cm⁻¹ toward the experimental 929.7. **THE FINDING: the original
  gate does not hold.** CASSCF orbital optimization tried first and rejected — numerically unstable
  for Be₂'s near-degeneracy (unconverged at several basis/R points; well depth swung
  non-monotonically 792→251 cm⁻¹ across one basis step). Falling back to the validated fixed-orbital
  CASCI(4,8) (`study_be2.py`) + NEVPT2 (core unfrozen, so core-valence correlation is automatic),
  CBS-extrapolated cc-pVTZ/QZ: **Re=2.58 Å, De=469 cm⁻¹** — genuinely bound at the right general
  location (unlike the baseline, whose "~305 cm⁻¹" turns out to be a spurious minimum at R≈4.5 Å,
  not a real near-Rₑ well) and closer to experiment than the baseline (461 cm⁻¹ off vs 625 cm⁻¹
  off — gated), but **not** within the original 100 cm⁻¹/0.1 Å tolerance — attributed to 2nd-order
  MRPT2 + an 8-orbital active space missing what CCSD(T)-F12/large-CAS-MRCI+Q needed in the
  literature. Gates G1–G4 in `tests/test_be2_cbs_spec.py`; `be2_cbs.py`.
  → [`SPEC_be2_cbs.md`](SPEC_be2_cbs.md) (reproduction attempt that falls honestly short, not a
  929.7 cm⁻¹ match; isolated dimer, no vibrational correction).
- [x] **1D Hubbard vs Bethe ansatz** — the validated stack (`hubbard_chain_integrals` → FCI/Krylov)
  reproduces the exact Lieb–Wu Bethe-ansatz energy: closed-shell-BC finite-L FCI extrapolates to the
  TDL integral (U=2→1.6, U=4→5.5, U=8→1.4 mHa/site), free-fermion (−4/π) + dimer limits at machine
  precision. *Gate revised from 1 mHa/site to < 8 mHa/site:* L≤12 FCI extrapolation is finite-size-
  limited (intermediate coupling slowest); sub-mHa needs DMRG at larger L. Finding: HF-referenced
  Krylov is slow on the strongly-correlated Mott chain (L≥6, large U). Gates G1–G4 in
  `tests/test_hubbard_bethe_spec.py`; `model_hamiltonians.py`.
  → [`SPEC_hubbard_bethe.md`](SPEC_hubbard_bethe.md) (repro of Lieb & Wu 1968, exact analytic ref).
- [x] **Nb₃X₈ / Hubbard model loader** — *Claim:* a tight-binding hopping matrix + Hubbard/cRPA
  interaction maps onto the universal `(h1, eri, e_core, nelec, norb)` interface and the
  number-conserving solver reproduces the analytic 2-site Hubbard dimer `(U−√(U²+16t²))/2` and PySCF
  FCI for a rank-4 cRPA tensor. *Gate:* `|E_solver − E_analytic| < 1e-6 Ha` across U/t = 0…20;
  rank-4 mapping `< 1e-9` vs FCI. → [`SPEC_nb3x8_hubbard.md`](SPEC_nb3x8_hubbard.md). **G1–G5 green
  (2026-06-29, re-confirmed 2026-07-12); `hybrid_quantum_solver/model_hamiltonians.py`.** Validates
  the Nb₃X₈ Model-Database
  drop-in path (handoff). Honest: not yet checked against the DB's DMFT/Hubbard-I gaps (DB files not
  bundled); Nb₃I₈ is the weakly-correlated anchor, not a strong-correlation headline.
- [x] **DMRG-referenced transition-metal active space** — NbN CAS(14,14) (half-filling sector
  comb(14,7)² ≈ 1.18×10⁷ determinants, beyond the 5×10⁶ FCI cutoff): two independent sweep
  schedules (perD 400/800/1200 vs ramp 300/600/1200, distinct seeds/scratch) agree to
  **3×10⁻⁸ Ha** — five orders below the 1 mHa gate — at **E = −110.046028 Ha** (high-spin
  nelec=(10,4)). **Findings:** (i) the CAS is a *soft* target — the high-spin sector converges to
  sub-nHa by D=400 (discarded weight ~1e-9), so "FCI-intractable by count" ≠ "strongly
  correlated"; (ii) a **near-degeneracy the SCF spin scan doesn't surface** — the low-spin
  nelec=(7,7) sector sits just **3.5 mHa above** the ground, so the reference is only meaningful
  once the sector is named. CI gates use cheap dims (≤300, ~2 min): |E_A′−E_B′| = 1.2 µHa, both
  in the discarded-weight regime. Gates G1–G3 in `tests/test_nbn_dmrg_reference_spec.py`;
  `nbn_dmrg_reference.py`. → [`SPEC_nbn_dmrg_reference.md`](SPEC_nbn_dmrg_reference.md) (reference
  number, not a materials claim — finite cluster, ECP, fixed geometry; a hard multireference TM
  benchmark needs the low-spin sector at real bond dimension — a follow-up).
## Done

- [x] **Nb₃X₈ strain / pressure response — Grüneisen parameters in the hopping** *(prediction toward
  strain-tuning experiments)* — with |t| as the compression knob (uniaxial strain ↑ the dimerization
  overlap), the leading response γ = dln X/dln|t| of the exact dimer. **Three findings:** (i) the
  **spin gap stiffens everywhere** (γ_J > 0) and runs monotonically from the atomic-limit **2**
  (Nb₃F₈, J∝t²) toward **1** (strong hopping) — 2.00→1.89→1.78→1.52, the halide series tracing the
  correlation crossover (closed-form dJ/dt == finite-diff); (ii) the **charge gap is non-monotonic**
  (a minimum at |t\*|=271→152→122→76 meV) and the family **straddles it** — γ_gap < 0 for F/Cl (below
  their minima, compression softens the Mott gap), > 0 for Br/I (above) — so spin & charge respond
  *oppositely only for the light halides*; (iii) **the sharp prediction:** Nb₃Cl₈ sits almost exactly
  at its charge-gap minimum (|t|=136 ≈ |t\*|=152), so its strain response is **spin-charge decoupled**
  — strong γ_J≈1.9, near-zero |γ_gap|<0.05 (>30× split): straining Nb₃Cl₈ moves its singlet-triplet
  gap (and, ∝J, its χ-max and Schottky peak) while barely touching its Mott gap. χ/Schottky Grüneisen
  == γ_J to <0.02% (Cl/Br), the iodide deviating 3–6% — the recurring E_s/J charge boundary. 100%
  reuse (`exact_charge_gap`, `dimer_exchange_analytic`, `chi_max_temperature`,
  `schottky_peak_temperature`). Gates G1–G4 in `tests/test_nb3x8_strain_spec.py`; `nb3x8_strain.py`.
  → [`SPEC_nb3x8_strain.md`](SPEC_nb3x8_strain.md) (isolated dimer, |t| the sole strain proxy, linear
  response; experimental hook: strain-tunable Nb₃Cl₈, `arXiv:2601.14524`).
- [x] **Nb₃X₈ magnetic heat capacity & entropy — Schottky anomaly + the R ln 2 plateau** *(thermo
  study, reproduction)* — completes the thermodynamic triad (χ(T) → C(T), S(T)) from the same exact
  N=2 trace. **Two findings:** (i) the magnetic Schottky peak sits at **T ≈ 0.352 J for the whole
  family** (pinned by the analytic two-level singlet/triplet result), a J-scale fingerprint whose
  ratio to the χ(T) peak is **universal & material-independent** — C-peak/χ-peak = 0.564/0.564/0.580
  ≈ 0.3515/0.625 (both exact two-level features); (ii) the localized-moment entropy plateau **R ln 4
  /dimer (= R ln 2/cluster)** is clean *only* when the charge scale E_s ≫ J — its flatness and its
  deviation from ln 4 worsen strictly Cl→Br→I (0.061→0.207→0.253; −1.3%→+1.6%→+17%) as E_s/J =
  16.9→8.1→3.1 shrinks, so the iodide has **no clean plateau** — the *same* charge-scale boundary
  that bounds the Bleaney–Bowers regime, now read off the entropy. Numbers `arXiv:2501.10320` did
  not report. 100% reuse (`n2_spectrum`, `chi_max_temperature`). Gates G1–G4 in
  `tests/test_nb3x8_thermo_spec.py`; `nb3x8_thermo.py`.
  → [`SPEC_nb3x8_thermo.md`](SPEC_nb3x8_thermo.md) (reproduction of the two-level Schottky/entropy
  laws with ab-initio params; isolated dimer, no lattice/phonon C, no 90 K transition).
- [x] **single-ramp DMRG extrapolation** — → [`SPEC_singleramp.md`](SPEC_singleramp.md). One ramping
  run via block2 `get_dmrg_results()` agrees with the per-D protocol < 0.1 mHa, lands at FCI, and
  uses half the sweeps (gates G1–G3 in `tests/test_singleramp_spec.py`). `protocol="ramp"`.
- [x] **Hₙ thermodynamic limit** — bond-dim + n→∞ extrapolation. → [`SPEC_hchain_tdl.md`](SPEC_hchain_tdl.md)
  (e_∞ = −0.539967 ± 0.000107 Ha/atom; gates G1–G5 in `tests/test_hchain_tdl_spec.py`).
- [x] **SKQD reproduces the exact-Krylov floor** — sample-based Krylov in determinants sampled from
  e^(−ikΔtH)|HF⟩ converges to FCI from above on H₄ / N₂ CAS(6,6) (`E_skqd ≥ E_fci − 1e-6`,
  `|E_skqd − E_fci| < 1.6 mHa` at depth ≥ 6). Gates G1–G4 in `tests/test_skqd_spec.py`;
  `hybrid_quantum_solver/skqd.py`. → [`SPEC_skqd.md`](SPEC_skqd.md) (repro of `arXiv:2501.09702`,
  `2508.02578`; no advantage at this scale).
- [x] **Factorization-native λ + symmetry shift** — native DF 1-norm (no brute-force Pauli) + a
  number-operator shift drops λ_DF on N₂ CAS(6,6) from 24.94 → 4.00 Ha (84%, gate ≥ 20%) with FCI
  invariant < 1e-8 Ha. Gates G1–G4 in `tests/test_scdf_lambda_spec.py`; `df_factorization.py`.
  → [`SPEC_scdf_lambda.md`](SPEC_scdf_lambda.md) (repro of BLISS/SCDF — `arXiv:2403.03502`,
  `2412.01338`; closes the gap named in `lambda_ladder.py`).
- [x] **QKSD excited states** — the same real-time Krylov subspace carries the low-lying *excited*
  spectrum in its Ritz values. Every Ritz value is variationally above its exact target (Cauchy
  interlacing); on H₄ the lowest 3 of 12 HF-reachable states converge to < 0.05 mHa by depth M=24
  (excited states need deeper M than the ground state — the recorded finding), and the first
  excitation gap matches FCI < 1.6 mHa. Gates G1–G4 in `tests/test_qksd_excited_spec.py`;
  `solve_excited` in `quantum_krylov_solver.py`. → [`SPEC_qksd_excited.md`](SPEC_qksd_excited.md)
  (repro of `arXiv:2109.06868`; exact statevector, no advantage at this scale).
- [x] **QKSD molecular properties** — the Krylov eigenstates (not just energies) give dipoles,
  transition dipoles, and oscillator strengths via ⟨Ψ_m|μ̂|Ψ_n⟩, matching dense-diagonalization
  FCI: HeH⁺ permanent dipole + bright transition (|μ|≈0.85 a.u.), H₂ recovered dipole-zero/dark by
  symmetry; property matrices Hermitian, eigenstates normalized. Gates G1–G4 in
  `tests/test_qksd_properties_spec.py`; `eigenstates`/`ritz_pairs` + `qksd_properties.py` +
  `build_dipole_operators`. → [`SPEC_qksd_properties.md`](SPEC_qksd_properties.md) (repro of
  `arXiv:2501.05286`; property values, not the hardware RDM/QSP measurement scheme).
- [x] **Tensor hypercontraction (THC) factorization + λ** — `(pq|rs) ≈ Σ χχ ζ χχ`; linear-LS THC
  reconstructs exactly at rank norb(norb+1)/2 (28<196 for H₂O vs DF-THC), FCI preserved, and
  `thc_lambda` == `df_lambda` on the DF-derived structured THC (validates the 1-norm against vetted
  code). **Finding:** with unoptimized collocation λ_THC ≈ 62× λ_DF — the THC λ advantage needs
  ISDF/optimized points (out of scope). Gates G1–G4 in `tests/test_thc_lambda_spec.py`;
  `thc_factorization.py`. → [`SPEC_thc_lambda.md`](SPEC_thc_lambda.md) (repro of `arXiv:2011.03494`;
  the asymptotic λ win is not reachable at this scale with naive collocation — the recorded boundary).
- [x] **Excited-state QKSD under shot noise** — `solve_excited` under the finite-sampling noise model
  degrades gracefully: H₂ gap error shrinks 0.062→0.012 Ha as shots go 4096→262144, stays bounded
  (no blow-up). **Finding:** excited states are ≈ 24–36× more noise-fragile than the ground state,
  and weakly-overlapped excited states fall below the noise-aware overlap floor (rank collapse) at
  low shots. Gates G1–G4 in `tests/test_qksd_noise_spec.py` (no new code — reuses the noise
  machinery). → [`SPEC_qksd_noise.md`](SPEC_qksd_noise.md) (repro of QKSD sampling-error analysis,
  `arXiv` Lee-Lee-Huh / Kirby 2024; idealized i.i.d. shot noise).
- [x] **Finite-T magnetic susceptibility χ(T) of the Nb₃X₈ dimers** *(thermodynamic study,
  reproduction)* — exact N=2 Boltzmann trace over the same downfolded cluster gives χ(T),
  Curie–Weiss θ, and room-T moments for all 10 sets from the cRPA parameters (numbers
  `arXiv:2501.10320` didn't tabulate). **Finding:** the pure-spin Bleaney–Bowers law holds up to
  a temperature set by the **charge scale E_s (first ionic singlet), not J** — the 5%-deviation
  onset is **T ≈ 0.40·E_s universally** (0.407–0.444 across the family), so the iodides
  (E_s/J ≈ 3, spin/charge least separated) break BB at the lowest reduced temperature (T₅/J ≈ 1.3
  vs > 3 for Br, > 10³ for F, strictly ordered I<Br<Cl<F). θ_CW = −J/4 (gated); room-T (300 K):
  fluorides/HT near free-pair Curie (χ≈1.9e-2 meV⁻¹, μ_eff≈1.22), iodides deep singlet
  (χ≈6e-6, μ_eff≈0.02). Gates G1–G4 in `tests/test_nb3x8_susceptibility_spec.py`;
  `nb3x8_susceptibility.py`. → [`SPEC_nb3x8_susceptibility.md`](SPEC_nb3x8_susceptibility.md)
  (**reproduction** of Bleaney–Bowers 1952 / the exact Hubbard-dimer thermodynamics
  `arXiv:1502.05038` with ab-initio params; prior Nb₃X₈ χ(T) exists — Sheckelton 2017, Haraguchi
  2017, Grytsiuk/Rösner `arXiv:2305.04854`; isolated single dimer, no 90 K transition).
- [x] **The visibility law, made predictive (a calibrated shot-cost law)** — the rule recorded
  qualitatively in three specs becomes an experiment-planning tool: for the unnormalized
  correlator, a line of weight w costs **shots* ∝ 1/(w²K)** (σ* = w√(dm)/(c(√d+√m))). Gated:
  log-log slope **−2.000** over four orders of magnitude in w (eight orders in shots — the
  Nb₃I₈ optical line costs ~2 shots/element, the near-dark Nb₃F₈ line ~2·10⁸); **one
  calibration on Nb₃Br₈ predicts every other line to ~10%**, including all three components of
  the multi-line H₄ removal signal; the edge component scales as 1/K (3.83 vs law 4).
  **Findings/boundaries:** line *attribution* needs a tolerance below the line spacing or the
  strong neighbor masquerades (45× too-early apparent detection — gated); tight attribution
  adds a pole-accuracy cost at shallow K (61% at K=8, 4% at K=32 — G3 split accordingly, the
  SDD-honest revision); false-positive rate 0.5%. Gates G1–G4 in
  `tests/test_visibility_law_spec.py`; `visibility_law.py`.
  → [`SPEC_visibility_law.md`](SPEC_visibility_law.md) (known-σ Gaussian noise; prices
  detection, not precision; c calibrated, not derived; fold device damping into w first).
- [x] **Spin spectroscopy — the interlayer exchange J of the Nb₃X₈ dimers** — the response
  trilogy's third channel: the staggered-magnetization kick S₁ᶻ−S₂ᶻ (cannot move charge) exposes
  exactly the state P leaves dark — the m=0 triplet — as a **single line at ω = J**, pinned by
  the closed form √(((U₀−Us)/2)²+4t²) − (U₀−Us)/2. The family's interlayer exchange constants
  (unreported in `arXiv:2501.10320`): **J = 0.051 / 66.2 / 119.1 / 245.9 meV (F/Cl/Br/I)**.
  **Findings:** the Heisenberg superexchange 4t²/(U₀−Us) fails progressively — 0% (F) → 6.3% →
  14.1% → **46.5% (I): the iodide is beyond the Heisenberg regime**; channel complementarity is
  exact (spin and optical lines ≥ 528 meV apart — symmetry partitions the spectrum); the
  local-moment fraction ‖Sᶻ|ψ₀⟩‖² falls 1.000 → 0.759 F→I (charge fluctuations eat a quarter of
  the iodide's moment). 100% machinery reuse (`absorption_lines`; Sᶻ|ψ₀⟩ is an exact eigenstate,
  re-exercising the degenerate-reference fix on a second operator). Gates G1–G4 in
  `tests/test_odmd_spin_spec.py`; `odmd_spin.py`. → [`SPEC_odmd_spin.md`](SPEC_odmd_spin.md)
  (isolated-dimer interlayer J only — no kagome in-plane exchange; Nb₃F₈'s 0.051 meV is below
  the model's own neglected terms — quote as ≈0).
- [x] **Optical absorption & exciton binding via ODMD** *(+ the eigenstate-kick hang, found and
  fixed)* — the two-particle side of ODMD spectroscopy: kicking ψ₀ with a same-sector operator
  (μ̂ for HeH⁺ — machine-exact lines, the bright transition cross-pins `SPEC_qksd_properties`'
  0.85² = 0.7224; P = n₁−n₂ for the dimers) gives the absorption spectrum with selection rules
  as missing lines. On the inversion-symmetric Nb₃X₈ dimers **the whole optical spectrum is one
  bright line** (the odd singlet at exactly U₀ — gated as a selection rule), so the optical gap
  is analytic and the **exciton binding Δc − Δopt collapses from 0.986·Us (F, atomic limit) to
  0.263·Us (I)** — the exciton unbinds with hopping; oscillator weight spans 1.1e-4 → 0.96
  (a 4-orders polarizability ladder). Numbers `arXiv:2501.10320` did not report. **Found
  defect:** `odmd_spectral.reference_signal` hung (~10¹⁷ expm substeps) when a kick lands on an
  exact eigenstate — width 0 → τ = π·10¹²; fixed with a degenerate-reference short-circuit,
  gated (G2), pinned spectral gates re-run green. Gates G1–G4 in
  `tests/test_odmd_optical_spec.py`; `odmd_optical.py`.
  → [`SPEC_odmd_optical.md`](SPEC_odmd_optical.md) (isolated-dimer optics, P as dipole stand-in;
  cluster gap difference, not a solid-state exciton).
- [x] **ODMD spectroscopy — Green's-function poles and weights from survival amplitudes** — a
  new observable class: the Lehmann representation is ODMD-shaped, so the survival amplitude of
  a_i|ref⟩ / a†_i|ref⟩ yields the photoemission / inverse-photoemission spectrum A(ω): DMD poles
  = ionization/affinity lines (machine-exact, ≤ 1e-13 Ha on strong lines), Vandermonde
  amplitudes = spectral weights (degeneracy-aggregated, < 1e-4 vs exact Lehmann). The Nb₃I₈
  dimer's A(ω) shows both Hubbard bands and its band gap min(ω⁺)−max(ω⁻) reproduces the capstone
  **842.44 meV to the last digit** (the cross-spec pin). **Findings:** HF-referenced weights
  differ from true Lehmann by ~2% on H₄ (real, bounded, gated); **intensities are
  damping-immune along with energies** (< 1e-6 at f=0.7 — uniform damping enters |λ| only,
  extending the device-ODMD phase immunity to the whole spectrum); the weakest satellite
  (w≈0.002) needed K=32 — the visibility law. Gates G1–G4 in `tests/test_odmd_spectral_spec.py`;
  `odmd_spectral.py`. → [`SPEC_odmd_spectral.md`](SPEC_odmd_spectral.md) (exact statevector
  signals; degenerate lines merge with summed weight; repro-adjacent — real-time a_i|ψ⟩
  propagation is established (Kosugi–Matsushita PRA 101, 012330), ODMD supplies the extraction).
- [x] **Coverage-gated error bars for ODMD (union bootstrap from a single signal)** — a real
  experiment gets one noisy signal and no ground truth; the union of a parametric bootstrap
  (refit → rebuild → re-noise at the *known* σ) and BOP-DMD-style bagging (random Hankel-column
  subsets, cf. arXiv:2107.10878) yields a 90% CI with measured coverage **0.895–1.000 on every
  system × budget** (H₂/H₄/N₂ CAS(6,6) × 10⁴/10⁵ shots, 200 trials each) at a conservatism cost
  of 1.9–2.9× the median error. **Findings:** each arm alone is broken in a complementary regime
  — the parametric bootstrap is anti-conservative up to 18× (coverage 0.05!) because the DMD fit
  absorbs realized noise and never sees threshold rank-switching (mean-shift bias correction
  doesn't fix it — probed 0.03), while bagging under-spreads few-mode signals (H₂ 0.785); and
  **no resampling of one signal can see model-misspecification bias** — at K=8 (8.7 mHa
  truncation bias) coverage collapses to 0.000, gated as the boundary: pair every interval with
  a depth-convergence check. Gates G1–G4 in `tests/test_odmd_uq_spec.py`; `odmd_uq.py`.
  → [`SPEC_odmd_uq.md`](SPEC_odmd_uq.md) (known-σ Gaussian noise, α=0.1, K=24; conservative by
  design — an anti-conservative error bar is worse than none).
- [x] **Nb₃X₈ cluster charge gaps through the simulated-hardware pipeline** *(capstone
  study/composition, not a method rung)* — the materials thread and the device-validated ODMD
  stack compose end-to-end: Δ = E(3)+E(1)−2E(2), one depolarizing-immune ground-state ODMD run
  per particle sector on genuinely-Trotterized Hadamard-test circuits under an Aer device noise
  model. Ladder (Nb₃I₈ LT-bulk, exact 842.44 meV): statevector ODMD exact to 1e-10; circuit
  eigenphases at reps=1 miss by **−101 meV (12%) — sector Trotter biases do NOT cancel in the
  gap** (order-2 ratio 4.65; near-commuting Nb₃F₈: 0.1 meV — the contrast); circuit-exact
  Richardson(2,4) lands 0.27 meV; the noisy device (cx=1e-4, 32768 shots) with Richardson(1,2)
  measures the gap to a **median error of 10.4 meV = 1.2%** (9.6× below raw reps=1).
  **Crossover finding:** at cx=3e-4 the noise floor exceeds the reps=2 bias and Richardson stops
  paying (16.3 vs 9.6 meV) — `SPEC_trotter_odmd` R1 demonstrated on a material. Gates G1–G4 in
  `tests/test_nb3x8_device_gap_spec.py`; `nb3x8_device_gap.py` (pure composition).
  → [`SPEC_nb3x8_device_gap.md`](SPEC_nb3x8_device_gap.md) (ISOLATED-CLUSTER gap — upper bound
  on the ~600–650 meV broadened solid, per the corrected `SPEC_nb3x8_gaps.md`; simulated device;
  nearly-free charged sectors — a pipeline demonstration, not a correlated-electron benchmark).
- [x] **Device-noise ODMD — eigenphases are depolarizing-immune (until the noise edge)** — a
  global depolarizing channel damps the signal s_k → f^k·s_k, multiplying every DMD eigenvalue
  by f but leaving its **phase** untouched: ODMD is *exactly* damping-invariant (< 1e-6 Ha even
  at f=0.7, where KQD on identically damped rows drifts 6.4 mHa noiselessly and fails by
  **~2600×/~150×** under 10⁵-shot noise at f=0.9/0.7). Local gate noise is *not* global: on the
  full Aer hardware stack (transpiled ancilla-controlled Trotter circuits, `hardware_krylov`
  Hadamard tests, depolarizing NoiseModel) the H₂ eigenphase survives **70% amplitude loss with
  0.05 mHa error** (cx=3e-4), and immunity ends exactly where the visibility law says — at
  cx=1e-3 the signal (2% retained) falls under the shot floor → 5 mHa. **Mechanism finding:**
  `odmd_energy`'s unit-modulus window misidentifies damped modes (3.75× worse); the device
  estimator = noise-edge cutoff + wide window + amplitude floor (pure composition, no new DMD
  code). Gates G1–G4 in `tests/test_device_odmd_spec.py`; `device_odmd.py` +
  `HardwareKrylovSolver.measure_signal`. → [`SPEC_device_odmd.md`](SPEC_device_odmd.md)
  (simulated device — no coherent errors/crosstalk/drift; energies are eigenphases of the
  *Trotterized* unitary, quote vs the circuit eigenphase or Richardson-remove per
  `SPEC_trotter_odmd`).
- [x] **Circuit-real ODMD + Richardson Trotter-bias removal** *(and a found+fixed silent-exactness
  bug)* — probing this spec exposed that `Operator()`/`Statevector.evolve()` evaluate an opaque
  `PauliEvolutionGate` via its **exact matrix**, ignoring the SuzukiTrotter synthesis: the repo's
  "Trotterized" statevector path (`TrotterKrylovSolver`) had been doing **exact evolution all
  along** (`order`/`reps` were no-ops; the old "within Trotter error" gate passed vacuously —
  QCIVET's failure mode in a new coat). Fixed by materializing the synthesized circuit in
  `build_trotter_step`; G1 is the regression gate (pre-fix deviation ~1e-16 vs gated > 0.05).
  With circuits now real: ODMD returns the ground eigenphase of U_trot to < 5e-11 Ha (DMD adds
  no approximation of its own); the eigenphase bias vs FCI obeys the order-2 law (ratios
  4.65/4.14 on H₂, 4.12/4.03 on H₄; 19.3→1.0 / 12.0→0.72 mHa over reps=1→4); two-point
  Richardson removes it to 0.046/0.0065 mHa (22×/112× below the fine bias) and beats the plain
  estimate ~9× under 10⁶-shot noise. **Boundary:** the price is circuit *depth* (reps=4 → 4×),
  and large-δt pairs leave higher-order residue (H₂ reps 1+2: 0.9 mHa) — extrapolate the fine
  pair, only when bias > noise. Krylov-subspace paths self-correct Trotter basis error (H is
  measured exactly), so the raw-eigenphase bias is the worst case. Gates G1–G4 in
  `tests/test_trotter_odmd_spec.py`; `trotter_odmd.py` + the `trotter_krylov.py` fix.
  → [`SPEC_trotter_odmd.md`](SPEC_trotter_odmd.md) (repro of Trotter effective-Hamiltonian theory
  `arXiv:1912.08854` + extrapolated eigenphases cf. `arXiv:2212.14144`; statevector circuits, no
  device noise).
- [x] **Certified two-sided energy brackets (Temple/Weinstein on Krylov Ritz states)** — every
  number in this repo was a variational *upper* bound; now each Krylov solve carries a rigorous
  **lower** bound too, from ONE extra expectation ⟨Ψ₀|H²|Ψ₀⟩: Temple's inequality on the QKSD
  ground eigenstate gives a certified bracket [E_Temple, E_Ritz] containing the exact
  reachable-sector energy at **every** system × depth tested (H₂/H₄/LiH/N₂ CAS(6,6) × M=2…24,
  zero escapes), closing to µHa width (N₂: 3.8 µHa at M=16). **Findings:** certification costs
  only ~2.7× the uncertified Ritz error at the same depth (gate < 5×); the oracle-free mode
  ε = θ₁ − σ₁ is valid at M ≥ 6 on all systems but its premise ε ≤ E₁ *fails at M ≤ 4* (the
  Krylov space must resolve the excited state before it can self-certify — the recorded
  boundary); Temple beats Weinstein ~430× on N₂. Gates G1–G4 in
  `tests/test_temple_bracket_spec.py`; `temple_bounds.py`.
  → [`SPEC_temple_bracket.md`](SPEC_temple_bracket.md) (repro of Temple 1928 / Pollak–Martinazzo
  applied to QKSD data; sector-restricted, exact statevector, ⟨H²⟩ hardware cost not modeled).
- [x] **Certified brackets on a spectral GAP (no FCI oracle)** — the Temple ground-energy bracket,
  lifted to the object spectroscopy actually measures: the reachable gap Δ = E₁−E₀ gets a two-sided
  interval from Krylov data alone, **Δ_hi = θ₁−τ₀** (Cauchy interlacing θ₁≥E₁ + Temple τ₀≤E₀) and
  **Δ_lo = (θ₁−σ₁)−θ₀** (Weinstein self-ε), at one extra ⟨H²⟩ on the first-excited Ritz state. The
  exact reachable gap is inside at **every M ≥ 6 across H₄ / LiH / N₂ CAS(6,6) (zero escapes)** and
  the interval **closes with depth** (H₄ 342→0.7, N₂ CAS(6,6) 391→37 mHa over M=6…24). **Finding /
  boundary:** certification inherits the temple_bracket M≤4 boundary *exactly* — at M=4 the Weinstein
  premise ε₁≤E₁ fails for the multireference cases (H₄, N₂) and the **lower** certificate escapes;
  the asymmetry is the mechanism — the upper certificate (interlacing+Temple) holds at every depth,
  the lower is premise-sensitive because a real-time Krylov space has no lower bound on E₂ to anchor
  a rigorous E₁ floor. 100% primitive reuse (`eigenstates`, `temple_bounds`). Gates G1–G4 in
  `tests/test_certified_gaps_spec.py`; `certified_gaps.py`.
  → [`SPEC_certified_gaps.md`](SPEC_certified_gaps.md) (extends Temple/Kato–Weinstein to gaps;
  sector-restricted, exact statevector, premise checkable-not-self-verifiable).
- [x] **Oracle-free trustworthiness certificate for the certified gap** — removes the open
  limitation of [`SPEC_certified_gaps.md`](SPEC_certified_gaps.md): its lower certificate rested on
  a premise (ε₁≤E₁) verifiable only against an FCI oracle. Replace the oracle with **cross-depth
  consistency** — a bracket is *corroborated* iff it overlaps the deep-anchor (intersection of the
  deepest brackets); a premise failure inflates/shifts the bracket so it no longer overlaps, and is
  caught with no FCI. **The certificate is ADAPTIVE, not a blanket "distrust shallow M":** on H₄/N₂
  it rejects M=4 (where certified_gaps escapes) but on LiH it *accepts* M=4 (premise already holds
  there) — it distinguishes real failures from merely-shallow-but-valid brackets. The self-checked
  interval (intersection of corroborated brackets) is validated to **contain the exact reachable
  gap** on all three (H₄/LiH/N₂) while a naive all-depth intersection is empty (the M=4 outlier).
  **Boundary:** necessary-not-sufficient — the deep anchor is trusted, and a consistently-biased
  sequence is invisible (the `SPEC_odmd_uq` model-misspecification blind spot); pairs with a
  depth-convergence check. 100% reuse (`gap_bracket_ladder`). Gates G1–G4 in
  `tests/test_gap_selfcheck_spec.py`; `gap_selfcheck.py`.
  → [`SPEC_gap_selfcheck.md`](SPEC_gap_selfcheck.md) (oracle-free self-verification; exact
  statevector; closes the certified_gaps premise-verifiability gap).
- [x] **Certified error bars on a molecular PROPERTY (dipole)** — extends the certified arc from
  energies/gaps to the observables spectroscopy reports. Rigorous interval [μ ± half_width] on the
  ground-state dipole from Krylov data, via Davis–Kahan sinθ ≤ σ₀/Δ_lo — where **Δ_lo is exactly the
  certified gap lower bound** of [`SPEC_certified_gaps.md`](SPEC_certified_gaps.md) — and the SHARP
  fluctuation bound half_width = 2σ_A·s + W_A·s² (using the dipole fluctuation σ_A, not ‖μ‖: LiH
  σ_A≈1.1 vs ‖μ_z‖≈6.9, ~6× tighter). **Zero FCI-dipole escapes** on HeH⁺/LiH; the interval closes
  to **−1.818 ± 0.065 a.u.** at M=24 (exact −1.817). **The finding — a property certificate INHERITS
  the gap certificate:** half_width is finite iff s<1 iff σ₀<Δ_lo, so it is vacuous exactly where
  Δ_lo is weak (LiH M=8–16) and sharp where healthy (M≥20) — pair with `gap_selfcheck` to know when
  Δ_lo is trustworthy. **Bug caught:** the reference must be the HF-*reachable* ground, not the
  global lowest eigenvector (a different particle-number sector for the charged HeH⁺ — G1 pins it).
  100% reuse (`gap_bracket`, `_mean_and_variance`, `build_dipole_operators`). Gates G1–G4 in
  `tests/test_certified_dipole_spec.py`; `certified_dipole.py`.
  → [`SPEC_certified_dipole.md`](SPEC_certified_dipole.md) (Davis–Kahan on QKSD + certified gap;
  sector-restricted, exact statevector, inherits the certified_gaps premise).
- [x] **Certified error bars on a RELATIVE energy (reaction / dissociation)** — the certified arc
  reaches chemistry's currency: relative energies. Δ = E(B)−E(A) gets a rigorous interval
  [τ_B−ρ_A, ρ_B−τ_A] by composing the Temple/Ritz brackets at two geometries — no FCI. On the H4
  symmetric stretch (0.9→2.3 Å) the exact Δ = **8.2255 eV** is inside at every depth (**zero
  escapes**), two-sided and < 0.01 eV wide at M=20. **The finding — the error bar localizes at the
  strongly-correlated endpoint:** the equilibrium bracket closes ~25× faster than the stretched
  (0.001 vs 0.025 eV at M=6), and at intermediate M the stretched Temple lower bound goes vacuous →
  Δ carries a **one-sided (upper-only)** certificate there before closing two-sided — the harder
  geometry sets the temple premise regime, exactly as [`SPEC_certified_gaps.md`](SPEC_certified_gaps.md)
  / `gap_selfcheck` chart. Reaches experiment (reaction/dissociation energies are measured; the
  certificate contains the *in-basis* FCI, basis error separate). 100% reuse (`krylov_bracket`).
  Gates G1–G4 in `tests/test_certified_thermochem_spec.py`; `certified_thermochem.py`.
  → [`SPEC_certified_thermochem.md`](SPEC_certified_thermochem.md) (composition of two temple
  brackets; sector-restricted, exact statevector, in-basis).
- [x] **The certified energy bracket under shot noise — a probabilistic certificate** *(the hardware
  capstone of the certified arc)* — re-grounds the whole arc (all rungs rest on ⟨H⟩, ⟨H²⟩) under
  i.i.d. shot noise (stderrs set by the 1-norms λ_H, λ_{H²}). **Two surprising findings:** (i)
  **sampling breaks the certificate** — at converged depth the raw bracket covers E₀ only ~0.40 of
  the time and the *variational upper bound holds ~0.49* (a coin flip), because ρ₀→E₀ makes symmetric
  noise land below E₀ half the time (the tighter the Ritz state, the more fragile); (ii) **shots do
  NOT buy coverage** — raw coverage is *N-independent* (identical at N=10⁴/10⁶/10⁸): the variational
  knife-edge is structural, not a finite-sample effect. **The repair:** z·se inflation restores
  coverage ≥0.9 (≈0.98 at z=2, conservative like `odmd_uq`), and the inflated half-width scales as
  z·λ_H/√N (54→5.4→0.54 mHa as N×100) — **inflation buys coverage, shots buy tightness** (the
  shot-cost law, cf. the visibility law). **Boundary:** λ_{H²}≫λ_H (H4: 63 vs 10) so the Temple lower
  bound is the noise-expensive side; idealized i.i.d. Gaussian, oracle gap, can't see systematic
  bias. Reuses `_mean_and_variance`. Gates G1–G4 in `tests/test_certified_noise_spec.py`;
  `certified_noise.py`. → [`SPEC_certified_noise.md`](SPEC_certified_noise.md) (Monte-Carlo coverage;
  the probabilistic counterpart to the exact-statevector certified arc).
- [x] **The precision-cost crossover — near-term certified (1/ε²) vs FT-QPE (1/ε)** *(bridges the
  repo's two halves)* — certifying the energy to precision ε costs the near-term arc N=(z·λ_meas/ε)²
  shot-measurements (standard limit, exponent −2, from `certified_noise`) vs FT-QPE
  Q=π·λ_DF/(2ε) queries (Heisenberg, exponent −1), so the resource ratio ~1/ε and **FT wins the
  exponent**. **The finding — the FT win is the exponent, not the constant:** the *raw* qubitization
  λ_DF does NOT uniformly beat the measurement 1-norm — for **N₂ CAS(6,6) λ_DF=24.94 > λ_meas=22.84**
  — so double factorization alone doesn't shrink the FT constant below measurement; only the
  **symmetry shift** ([`SPEC_scdf_lambda.md`](SPEC_scdf_lambda.md)) drops λ_DF to 0.97/1.83/4.00,
  below λ_meas for every molecule by a margin that *grows* with size (2.8×→5.7×) — the shift is
  load-bearing, not a nicety. R@1.6mHa ~ 4·10³–3·10⁴; crossover ε* parametrized by the per-query
  cost. **Boundary:** shots-vs-queries needs a common cost model (per-query T-cost = chem-ft, not
  computed); the exponent gap is unit-independent. Reuses `df_lambda`/`symmetry_shift`. Gates G1–G4
  in `tests/test_precision_cost_spec.py`; `precision_cost.py`.
  → [`SPEC_precision_cost.md`](SPEC_precision_cost.md) (reproduction-adjacent scaling laws; the
  composition on real shifted λ's is new; absolute T-gate cost out of scope).
- [x] **Cost advisor — a per-molecule near-term-vs-FT verdict** *(the bridge, closed into a decision
  tool)* — turns [`SPEC_precision_cost.md`](SPEC_precision_cost.md) into an engineering artifact:
  `advise(λ_meas, λ_DF, ε, ρ)` returns which is cheaper (near-term certified Krylov vs FT-QPE) at
  target accuracy ε, under a common cost model whose only free parameter is the honest unknown
  ρ = cost_per_query/cost_per_shot. The verdict flips exactly at ε*(ρ) ∝ 1/ρ; the flip-ρ at fixed ε
  is N/Q; and `robust_over_rho` gives a cross-validated recommendation. **Quantified verdict:** for
  N₂ CAS(6,6) at chemical accuracy the flip-ρ is **~2×10⁵** — FT is cheaper unless a query costs more
  than ~2×10⁵ shots (ρ is never hidden — it's an explicit input). Reuses `precision_cost`. Gates
  G1–G4 in `tests/test_cost_advisor_spec.py`; `cost_advisor.py`.
  → [`SPEC_cost_advisor.md`](SPEC_cost_advisor.md) (resource-count decision under a stated cost
  model; true ρ hardware-specific, not computed).
- [x] **Excited-state ODMD via noise-edge thresholding** — the *same* survival-amplitude signal
  carries the low-lying spectrum in its higher DMD eigenphases (no extra measurements): noiseless
  E₁/gap < 1e-5 Ha (H₄ K=24, N₂ CAS(6,6) K=48). Under shot noise the ground-state spec's relative
  5σ·σ_max floor loses the excited mode in 100% of seeds; replacing it with the absolute
  random-matrix noise edge c·σ(√d+√m) recovers the H₄ gap to **5.8 mHa at 10⁵ shots/element —
  ~31×/21× below noisy QKSD `solve_excited` at M=16/24** (matched per-element σ, and QKSD is
  noiselessly converged there, so noise — not depth — is its limit). **Finding (the visibility
  law):** mode n is recoverable iff p_n·√(dm) clears the noise edge, so depth buys excited
  visibility as √K (K=16: 74% unresolved; K=48: 0%); dark/weakly-overlapped states stay invisible
  (the rodeo/SKQD physics, now quantitative). Gates G1–G4 in `tests/test_odmd_excited_spec.py`;
  `odmd.py` (`noise_edge`/`odmd_spectrum`/`sample_odmd_spectrum`).
  → [`SPEC_odmd_excited.md`](SPEC_odmd_excited.md) (single observable; the MP edge on *structured*
  Hankel noise is a calibrated heuristic, c=1.2).
- [x] **ODMD — ground state from the survival amplitude alone** — the complex overlap time series
  s_k = ⟨φ₀|e^(−ikτH)|φ₀⟩ (the *first row* of the S matrix QKSD already measures; **no Hamiltonian
  elements at all**) recovers FCI via SVD-truncated Hankel DMD: < 1e-5 Ha noiseless at K=20 on
  H₂/H₄/N₂ CAS(6,6), and at matched measured-element count + shots the median error beats KQD by
  ~57× (10⁴ shots) / ~10× (10⁵) on N₂ — the λ-scaled H-measurement noise simply never enters.
  **Findings:** ODMD is *non-variational* (dips ~2.6 mHa below FCI at K=8 — never quote it as a
  bound), and the SVD truncation IS the robustness (removing it inflates the noisy median ~900×).
  Gates G1–G4 in `tests/test_odmd_spec.py`; `odmd.py`. → [`SPEC_odmd.md`](SPEC_odmd.md) (repro of
  `arXiv:2306.01858`; complex-quadrature signal, exact statevector, i.i.d. shot noise; advantage
  measured against this repo's LCU/Hadamard KQD noise model — MSD narrows it).
- [x] **Mirror subspace diagonalization (MSD)** — estimates H from central finite-differences of
  shifted-time overlaps instead of per-Pauli measurement, so sampling variance scales with the
  stencil 1-norm fd1 not the Hamiltonian 1-norm λ. With an energy-level shift + order-8 stencil,
  fd1=5.48 < λ=14.75 on N₂ CAS(6,6) and MSD's median error is ≈ 3.2× below KQD at 10⁵ shots.
  **Boundary:** H₂ (λ/W≈1.3) gives fd1>λ → no advantage; the win is a λ/W effect, modest at this
  scale (the paper's 10–10⁴× needs larger λ/W). Gates G1–G4 in `tests/test_msd_sampling_spec.py`;
  `msd.py`. → [`SPEC_msd_sampling.md`](SPEC_msd_sampling.md) (repro of `arXiv:2511.20998`; idealized
  shot noise, exact statevector).
- [x] **Hamiltonian-moment energies (PDS / CMX)** — ground-state energy from the moments ⟨H^n⟩ of
  the HF reference, no time evolution. PDS(K) is a variational upper bound (≥ FCI at every K) that
  converges (H₄ 67→10→2.0→0.43 mHa over K=1..4; PDS(1)=⟨H⟩); CMX(2) dips below FCI on H₂
  (−0.27 mHa) — non-variational, the recorded boundary. Gates G1–G4 in
  `tests/test_moment_pds_spec.py`; `moment_expansion.py`. → [`SPEC_moment_pds.md`](SPEC_moment_pds.md)
  (repro of PDS/CMX — `arXiv:2101.08526`, JCP 153 201102; exact statevector moments).
- [x] **Classical shadows** — estimate ⟨H⟩ from randomized single-qubit (random-Pauli) measurements.
  Unbiased on HF and FCI states (within 4·stderr), ~1/√shots convergence, single-shot variance
  bounded by the HKP shadow norm Σ|c_k|²3^{w_k} (1.95 ≤ 3.09 on H₂). **Finding:** the 3^{weight} factor
  makes high-weight terms sample-expensive (weight≥3 ≈ 22% of the norm). Gates G1–G4 in
  `tests/test_classical_shadows_spec.py`; `classical_shadows.py`.
  → [`SPEC_classical_shadows.md`](SPEC_classical_shadows.md) (repro of HKP — `arXiv:2002.08953`;
  random-Pauli shadows, exact statevector).
- [x] **Rodeo algorithm** — stochastic spectral filter: K cycles of random-time evolution + ancilla
  band-pass the spectrum to a target E. The expected survival probability peaks at the eigenvalues
  with height = reference overlap; the dominant low-energy peak recovers FCI (≈0.13 mHa, H₂/H₄), the
  peak sharpens with K (H₄ 4.63→0.13 mHa over K=3→12), off-resonance suppressed as (<1)^K. **Finding:**
  ground-peak height = |⟨HF|E_0⟩|², so a poor reference gives a weak peak. Gates G1–G4 in
  `tests/test_rodeo_spec.py`; `rodeo.py`. → [`SPEC_rodeo.md`](SPEC_rodeo.md) (repro of
  `arXiv:2009.04092`; expected-value simulation, exact statevector).
- [x] **Quantum imaginary-time evolution (QITE)** — reach the ground state by replacing the
  non-unitary step e^{-ΔτH} with a unitary e^{-iΔτÂ}, Â from the McLachlan system
  S a = b (S_IJ=Re⟨σ_Iσ_J⟩, b_I=Im⟨σ_I H⟩). Exact ITE monotone/variational →FCI; full-domain QITE
  reproduces it and reaches FCI on H₂ (update equations correct); step error → 0 with Δτ. **Finding:**
  a weight-≤2 domain stalls at Hartree–Fock (+20.5 mHa) — H₂'s correlation is a weight-4 operator, so
  QITE accuracy is set by the domain (Motta's locality). Gates G1–G4 in `tests/test_qite_spec.py`;
  `qite.py`. → [`SPEC_qite.md`](SPEC_qite.md) (repro of `arXiv:1901.07653`; exact statevector, full
  4ⁿ domain so H₂-scale).

- [x] **Exact Nb₃X₈ cluster gaps vs Hubbard-I** *(scientific study, not a method rung)* — the Nb₃X₈
  bilayer downfolds to an exactly-diagonalizable generalized Hubbard dimer; from the paper's cRPA
  parameters, exact ED gives charge gaps the paper never reported (Nb₃I₈ bulk 842 / bilayer 1961 …
  Nb₃F₈ bulk 2581 / bilayer 3979 meV). **Corrected conclusion:** on the isolated cluster Hubbard-I
  underestimates the iodides by ~29% (bulk)/12% (bilayer) and the single-ratio U₀/|t| law fails
  (Spearman −0.86) — but a thermodynamic-limit check (coordination scan + DMRG chain) shows this is an
  artifact of neglecting band broadening: restoring coordination collapses the exact Nb₃I₈ gap
  (842→650 at z=3; ~708 at the 1-D DMRG TDL; ~600–650 with 3-D coordination) back toward Hubbard-I
  (~599). So it **supports** the paper's cluster-DMFT/Hubbard-I gaps — a null/confirming result whose
  headline was corrected by the TDL step (the honest turn). Gates G1–G6 in
  `tests/test_nb3x8_gaps_spec.py`; `nb3x8_gaps.py`. → [`SPEC_nb3x8_gaps.md`](SPEC_nb3x8_gaps.md)
  (data from `arXiv:2501.10320`, Tables I & IV).
- [x] **ScreenLoop — screening with *zero false eliminations*** — *Claim:* eliminating candidates only
  when their certified bracket is strictly disjoint from the target region can never discard a true
  hit, a guarantee no point-estimate screener can make. *Gate (the check that could kill it):* on a
  50-candidate space with exhaustive ground truth, `screen` eliminates **0** of the 11 true hits
  (pilot precision and through refinement) while `point_estimate_screen` eliminates > 0 on the *same*
  space with the *same* oracle; the loop reaches the exhaustive-at-full-precision hit set for < half
  the cost. Soundness is *relative to* bracket validity — a lying oracle voids it (stated, not hidden).
  Gates 1–5 in `tests/test_screenloop_spec.py`; `screenloop/`.
  → [`SPEC_screenloop.md`](SPEC_screenloop.md).
- [x] **Shift BOTH sides — the bridge overstated FT's advantage** — the BLISS/SCDF number-operator
  shift that the bridge banked on the FT side (λ_DF) **also lowers the measurement 1-norm λ_meas**,
  so it cuts the near-term shot cost too; comparing raw-λ_meas against shifted-λ_DF was an
  inconsistency. **THE FINDING:** the draft's headline (λ_meas −54/−51/−73%, crossover 4–14×) was
  scored with the *wrong* 1-norm — `precision_cost.measurement_lambda` sums **every** Pauli
  coefficient *including identity*, but the identity is a zero-variance constant that costs **zero
  shots**, and a large part of what the shift does is dump weight into it (N₂: 8.55 → 0.10 Ha). On
  the honest identity-excluded 1-norm the reductions are **42.9/39.4/57.9%** and the fair crossover
  moves **2.7–5.7×**, not 4–14× (N₂ flip-ρ 8.14e4 → 1.44e4). The claim survives, smaller; G1's bar
  was revised 40%→35% because H₂O honestly reduces by 39.4%. Gates G1–G5 in
  `tests/test_shift_both_sides_spec.py`; `shift_both_sides.py`.
  → [`SPEC_shift_both_sides.md`](SPEC_shift_both_sides.md).
- [x] **λ_meas includes the identity term — the published crossovers inherit the bias** *(follow-up
  from `SPEC_shift_both_sides`)* — fixed in place: `precision_cost.measurement_lambda` now excludes
  the zero-variance identity term by default (`include_identity=True` retained for archaeology).
  Re-scored: λ_meas 2.70/9.93/22.84 → 1.89/5.40/14.30 (identity was 30/46/37% of the 1-norm); the
  shift margin 2.8/5.4/5.7× → 1.9/3.0/3.6× (precision_cost G3 bar revised 5.0 → 3.0, reason
  recorded); ε\* and R move by (λ_old/λ_new)². **THE FINDING: one qualitative verdict flips** —
  H₂ at ρ=10⁴ goes FT → near-term; the conjecture recorded here earlier ("no qualitative verdict
  flips") was **wrong**. Every structural finding of `SPEC_precision_cost` survives (raw λ_DF > λ_meas
  for N₂ holds *more strongly*; margin still grows with size). Gates G1–G4 in
  `tests/test_lambda_meas_identity_spec.py`.
  → [`SPEC_lambda_meas_identity.md`](SPEC_lambda_meas_identity.md).
- [x] **`nb3x8_device_gap` G2: a Trotter-leakage resolution floor, exposed as a flaky gate** —
  resolved by [`SPEC_trotter_resolution_floor.md`](SPEC_trotter_resolution_floor.md). Root cause,
  measured: **(1)** Trotter synthesis was **nondeterministic across processes** — the Pauli term
  order followed Python hash randomization (observed orderings ≡ coeff-asc/coeff-desc; F8 sector-2
  reps=1 deviation 1.045e-2 vs 5.964e-3), and SuzukiTrotter products depend on term order; **(2)**
  the F8 sector-2 reps=1 eigenphase sits **below the method's resolution floor** — genuine reference
  population 1.36e-5 < Trotter leakage dev² ≈ 3.6e-5–1.1e-4 — so its computed population was
  interference noise (1e-5…2.5e-13) and the extracted branch flipped by the 2π/τ wrap quantum
  (−1171.1 vs +2580.7 meV; flake rate 3/5 pytest, 2/2 make gates). **Fix:** canonical term ordering
  (largest-|coeff| first) in `build_trotter_step` — device_gap now 4/4 green on every run — plus the
  law as a gated invariant: *an eigenphase is resolvable only if the reference population exceeds
  ‖U−U_exact‖²* (`trotter_resolution_floor.is_resolvable`; the ordering-spread probe: 3756.7 meV at
  reps=1 = the wrap quantum, 5.5 meV at reps=2). **Casualties honestly recorded:** the spec's prose
  values were hash artifacts (I8 bias −100.9/12%/4.65 → deterministic −292.9/35%/4.29) and G4's
  noise crossover moved a decade (cx=3e-4 → 1e-3; Richardson still pays at 3e-4 and 6e-4) — the
  crossover *exists* as claimed, its recorded location was never reproducible. Latent `np.min`
  unphysical-branch issue fixed below (G5). Gates G1–G4 in
  `tests/test_trotter_resolution_floor_spec.py`; `trotter_resolution_floor.py`.
- [x] **`trotter_odmd` eigenphase selection can pick unphysical branches** — *(latent, from the
  resolution-floor work)* `build_trotter_odmd_problem` took `min(-angle(λ)/τ)` over the population
  cut alone; eigenphases outside the reachable band exist in the same spectrum (an e = −1657.45 meV
  image for F8 sector-2, |e| > width/2) and `np.min` would select one if it ever cleared the
  population cut. **Fixed:** extracted `select_ground_eigenphase` and restricted selection to the
  physical band `|e| ≤ width/2` (with a relative-tolerance boundary guard — a naive strict `<=`
  drops the closed-boundary branch to roundoff and silently flips I8 sector-1 reps=2's sign, caught
  while verifying "unchanged"). Pure hardening: every committed F8/I8 `circuit_gap` number
  reproduces exactly; the exclusion only bites a synthetic decoy branch (gated directly, since no
  real committed system exercises it — "untriggered under canonical ordering" confirmed). Gates
  G5 in `tests/test_trotter_resolution_floor_spec.py`; `trotter_odmd.py`.
- [x] **SenseForge — Nb₃X₈ strain/field sensor screening: the harness is real, the *ranking* is
  vacuous** — a resumable sweep harness (validated YAML + provenance hash, crash-resume to a
  byte-identical CSV, ADR-0003 header automation, closed-form-vs-exact cross-check on all four
  halides) over strain (`t(ε)=t₀(1+ε)`) and Zeeman (`g·μ_B·B·Sz`) perturbations of the validated
  dimer. **THE FINDING (G9):** the FoM ranking **does not discriminate on this model, on either
  axis.** Brackets are exact (zero-width) so `FoM ≡ |S|`, and |S| is *exactly constant* on the field
  axis (Zeeman response is exactly linear → **all 19 operating points tie**; the rank order is sort
  noise) and *strictly monotone* on the strain axis (→ **rank 1 is just the window edge**: widen
  ±2%→±5% and the "optimum" moves +1.75%→+4.75%). No interior optimum exists on either axis. This
  was **shipping as a real recommendation** — the pre-fix field design card published "+2 T" as the
  #1 operating point, provably no better than any other. Fixed: `ranking_verdict()` classifies every
  ranking (degenerate/monotone/interior) and every artifact now leads with the verdict. **Second
  defect:** `krylov_dim` was dead *and* stamped on every published design card as `krylov_dim=12` —
  false provenance for an exact-diagonalization result with no Krylov subspace; removed. Also
  recorded: `certified_gaps.gap_bracket` is *actively wrong* here (returns the bright ionic ~1117
  meV, not the dark spin gap J≈66 meV). The sensitivity itself is real and reproduces the gated
  Grüneisen sign. **Boundary:** a non-vacuous screener needs a system that is (a) too large for FCI
  yet Krylov-resolvable (so `/width` does work) *and* (b) has an interior |S| optimum — the Nb₃X₈
  dimer is neither. Gates G1–G9 in `tests/test_senseforge_spec.py`; `senseforge/`.
  → [`SPEC_senseforge.md`](SPEC_senseforge.md).

- [x] **The certified relative-energy bracket under shot noise — composition beats either endpoint
  alone** — extends `certified_thermochem`'s two-geometry composition into the noise regime
  `certified_noise` opened for a single bracket. **THE FINDING:** the composed relative-energy
  bracket inherits the coin-flip collapse (raw coverage ~0.70, shot-count-independent, same
  mechanism as the single-bracket ~0.40) but is **measurably better calibrated than either endpoint
  alone** — restoring 90% coverage costs **z\* ≈ 0.55–0.60**, roughly 3–4× *less* inflation than the
  single-bracket z=2 rule (which still works if reused unchanged, just over-conservative). Two
  independent one-sided coin-flips, composed into a two-sided difference, do not simply compound.
  100% reuse (`certified_noise`'s per-endpoint noise model + `certified_thermochem`'s composition);
  the one new seam is that cross-geometry energy offsets do NOT cancel (unlike a single-geometry
  gap) and must be added back before composing. Gates G1–G4 in
  `tests/test_certified_thermochem_noise_spec.py`; `certified_thermochem_noise.py`.
  → [`SPEC_certified_thermochem_noise.md`](SPEC_certified_thermochem_noise.md) (oracle E1 per
  geometry, i.i.d. Gaussian shot noise, H4 stretch only — the specific z* is a measurement on this
  system, not a derived constant; see R1).

- [x] **Gap self-check under shot noise — intersection concentrates noise where composition
  diluted it** — the noise counterpart to `certified_thermochem_noise`, but for the OTHER
  composition operator this repo uses: `gap_selfcheck`'s oracle-free trustworthiness certificate
  builds its self-checked interval by *intersecting* several independently-noisy corroborated
  brackets, rather than differencing two. **THE FINDING:** raw coverage is broken the same way
  (~0.14–0.20 H4, ~0.06–0.08 LiH, shot-count-independent — same coin-flip mechanism as
  `certified_noise`), but unlike the difference composition (which needed *less* than the
  single-bracket z=2 rule), intersection needs **MORE**: z=2 leaves coverage at 0.70 (H4) / 0.53
  (LiH) at shots=1e5, and the minimal z reaching 90% coverage is **z\* = 3.25 (H4) / 4.00 (LiH)** —
  roughly 1.6–2× the single-bracket rule, the opposite direction from
  [`SPEC_certified_thermochem_noise.md`](SPEC_certified_thermochem_noise.md). Padding does fix
  "inconclusive" well before "correct": the empty-interval rate (no bracket corroborates) drops
  from ~0.1–0.13 at z=0 to <0.005 at z=2, even while coverage is still broken there. 100% reuse
  (`gap_selfcheck.self_checked_gap` unmodified, `certified_noise.certified_half_width`); the one
  new design choice — a simple post-hoc pad on each depth's bracket rather than internal
  directional inflation — was forced by finding that the two-eigenstate gap formula has no single
  unambiguous "worst-case" direction for θ0 (it plays opposing conservative roles in gap_lower vs.
  the Temple term feeding gap_upper). Gates G1–G4 in `tests/test_gap_selfcheck_noise_spec.py`;
  `gap_selfcheck_noise.py`. → [`SPEC_gap_selfcheck_noise.md`](SPEC_gap_selfcheck_noise.md)
  (self-mode/oracle-free only; i.i.d. Gaussian shot noise; two systems — the z* numbers are a
  measurement, the *direction* is the falsifiable claim; no repair attempted, a follow-up).

- [x] **Certified dipole under shot noise — inflation cannibalizes the gap margin it needs to stay
  finite** — closes the noise trilogy: `certified_thermochem_noise` (difference composition, needs
  *less* z), `gap_selfcheck_noise` (intersection composition, needs *more* z), and now the property
  bracket, which composes THREE noisy quantities through the Davis–Kahan `s = sigma_0/Delta_lo < 1`
  gate. **THE FINDING:** on the healthy-margin system (HeH+) moderate inflation (z=1) still works
  (coverage 0.98–1.00 at shots ≥1e5), but at a tight shot budget (1e4) MORE inflation makes it
  *worse* — the finite-bracket rate falls from 0.728 (z=1) to 0.554 (z=3) — an **inflation ceiling**
  neither prior noise spec showed (both monotonic in z up to z=5–6): padding conservatively shrinks
  `Delta_lo`, the exact margin the certificate needs to stay finite, so widening the bracket to fix
  the coin-flip collapse cannibalizes the resource that keeps it from going vacuous. On the
  fragile-margin system (LiH at M=16) no tested z rescues it — finite-bracket rate stays <0.3
  throughout, a recorded boundary, not a fix. **Caught before it became a claim:** an earlier probe
  using internal directional perturbation (rather than post-hoc padding) produced a spurious
  monotonically-*decreasing* coverage curve from a sign error (θ0 is negative, inverting the
  naive "shift down is conservative" assumption) — the same class of ambiguity
  [`SPEC_gap_selfcheck_noise.md`](SPEC_gap_selfcheck_noise.md) flagged, now hit a second time and
  designed around the same way (simple post-hoc padding on raw, undirected noisy draws). 100% reuse
  (`certified_dipole.spectral_width`, `certified_gaps`' self-mode formula, `certified_noise`'s
  padding convention). Gates G1–G4 in `tests/test_certified_dipole_noise_spec.py`;
  `certified_dipole_noise.py`. → [`SPEC_certified_dipole_noise.md`](SPEC_certified_dipole_noise.md)
  (coverage target is the noiseless M=16 Krylov estimate, not dense-diagonalization FCI — isolates
  noise from depth-convergence; i.i.d. Gaussian noise; two systems, the ceiling's exact location is
  shot-budget/system-dependent — R1).

- [x] **ADAPT-VQE — gradient-greedy selection actually buys a more compact ansatz (on a system with
  real correlation)** — a different algorithmic family from everything else gated so far (variational
  ansatz growth, not Krylov/ODMD/DMRG/rodeo/QITE/shadows/moments): `adapt_vqe.py` already ran and
  asserted a variational floor informally in its `__main__`, but had never been CI-gated, and its
  core selling point (greedy > fixed order) had never been checked at all. **THE FINDING:** on H4
  CAS(4,4) (real multi-orbital correlation), greedy gradient selection reaches chemical accuracy in
  **9 operators**; NONE of 5 seeded random fixed orders reach it within a matched 14-op budget — a
  decisive win. **Boundary, recorded not smoothed over:** on LiH CAS(2,2) (trivially simple active
  space) adaptivity buys **nothing** — greedy and every one of 5 random orders reach chemical
  accuracy in exactly 1 operator. Variational floor pinned at every growth step (not just the final
  one) on all three systems (H2/LiH/H4), closing the gap between an informal `assert` and an actual
  CI gate. New code: `fixed_order_vqe` (~25 lines, the comparison baseline — byte-identical
  growth/optimize loop to `adapt_vqe`, differing only in selection rule), everything else reused
  unchanged (`qubitization_blueprint`'s JW operators, `build_pool`, `hf_state`). Gates G1–G4 in
  `tests/test_adapt_vqe_compactness_spec.py`; `adapt_vqe.py`.
  → [`SPEC_adapt_vqe_compactness.md`](SPEC_adapt_vqe_compactness.md) (exact statevector, generalized
  singles+doubles pool only, operator count not circuit depth/T-cost; the decisive gap is shown on
  one correlated system, not claimed universal — R1).

- [x] **Z2 qubit tapering preserves the FULL sector spectrum — independently verified, not just
  the ground energy** — `taper_qubits.py` had only ever checked its own ground-eigenvalue-matches-
  CASCI assertion informally, in `__main__` (the weakest possible falsifier: a bug scrambling every
  excited state while leaving the ground state alone would sail through it). **THE FINDING:** the
  full spectrum matches an INDEPENDENTLY-constructed sector projection (built directly from
  computational-basis parity, no Clifford rotation — deliberately not reusing
  `taper_hamiltonian`'s own code path) to machine precision on H2 CAS(2,2), LiH CAS(2,2), and — new
  scope — an OPEN-SHELL radical (H3, CAS(3,3), spin=1, nelec=(2,1)) that `taper_qubits.py`'s own
  `__main__` never exercises (it hardcodes closed-shell RHF); qubits removed exactly equals the
  independent Z-symmetry count on all three; and the check is proven non-vacuous — the same
  construction with a deliberately WRONG reference state gives a genuinely different spectrum.
  **Scalability boundary found while probing, recorded not fixed:** `pauli_decompose` (used inside
  `taper_hamiltonian` itself) costs **~1000 seconds at 8 qubits** vs ~0.02s at 4 and ~0.3s at 6 — a
  ~16x-per-qubit blowup — meaning `taper_qubits.py`'s own existing H4 CAS(4,4) `__main__` case
  already takes 15-20 minutes to run; this spec's gates stop at 6 qubits by design rather than
  pretend the module scales. No library code changes — the independent verification lives entirely
  in the test file, a genuine external check. Gates G1–G4 in `tests/test_taper_spectrum_spec.py`.
  → [`SPEC_taper_spectrum.md`](SPEC_taper_spectrum.md) (Z-type symmetries only; minimal-basis
  STO-3G; composing tapering with a downstream method like ADAPT-VQE is a follow-up, not attempted).

- [x] **The qubitization walk operator recovers EVERY Hamiltonian eigenvalue, exactly twice** —
  `qubitization_blueprint.py`'s own `verify_qubitization` only ever checked one direction of
  `eig(W) = e^{±i arccos(E_k/lambda)}`: that every recovered value is *valid* (near some true
  eigenvalue), never that every true eigenvalue is actually *recovered* — the direction that matters
  for QPE, since a construction bug that silently dropped the ground state from `W`'s spectrum would
  sail through the existing check. **THE FINDING:** the reverse holds to machine precision (H2
  2.0e-15, LiH 1.0e-15), and sharpened to an exact counting invariant — of `W`'s (non-deduplicated)
  eigenphases, exactly `2 * dim(H)` land within tolerance of a true eigenvalue on both systems (32 =
  2×16), i.e. every eigenvalue slot contributes exactly one `theta/-theta` pair, no eigenvalue
  missing, no extra collision; every individual multiplicity count is even, confirming the pairing
  holds without exception (H2/LiH counts up to 6, all even). Construction bookkeeping
  (`lambda == sum|coeffs|` independently recomputed, `dim(W) == ancilla_dim * dim(H)`) pinned too.
  **Bug caught while writing the gate, not shipped:** an early version of the counting invariant
  summed local match-counts from the (possibly-degenerate) `eigH` side, double-counting whenever two
  eigenvalues sat within tolerance of each other (60 instead of the correct 32) — fixed by counting
  from the `W`-recovered side instead, and left as a comment in the test explaining why the other
  direction is wrong. No library code changes — the check lives entirely in the test file, a genuine
  external verification of `qubitization_blueprint.py`'s existing `build_walk_operator`. Gates G1–G4
  in `tests/test_qubitization_spectrum_spec.py`.
  → [`SPEC_qubitization_spectrum.md`](SPEC_qubitization_spectrum.md) (dense-matrix verification
  oracle, not circuit-level; two CAS(2,2) systems; T-gate/query cost claims not gated here).

- [x] **Iterative QPE — a fat-tailed bit-flip cascade the median error hides** — unlike
  `qpe_walk_readout.py`'s noiseless exact-simulation oracle, `iterative_qpe.py`'s per-bit
  measurement is genuinely stochastic (a majority vote over `shots_per_bit` Bernoulli trials), and
  because IQPE reads bits least-significant-first with feedback, a single flipped LOW bit corrupts
  every subsequent (coarser) bit's measurement — a cascade the docstring's `"precision ~ 1/2^bits"`
  claim (checked only by eye on one seed) never characterized. **THE FINDING:** median error DOES
  improve monotonically with more bits at the module's own default `shots_per_bit=15` (a staircase,
  same character as `qpe_walk_readout`'s precision law), but at low `shots_per_bit` the WORST-CASE
  error over seeds is dramatically worse than the median — at `n_bits=8, shots_per_bit=3`:
  **max/median = 10.4×** (134.3 vs 12.9 mHa) — while at the SAME `n_bits` with the module's own
  default `shots_per_bit=15` the ratio is exactly 1 (fully deterministic, no cascade — the default
  already avoids the pathology a naive shot-cutting choice would reintroduce). The tail has a
  measured threshold, not an open-ended risk: at `n_bits=12`, `shots_per_bit=1` shows the cascade
  (ratio 10.3×) but `shots_per_bit>=3` is fully deterministic. **Boundary, recorded not smoothed
  over:** even at the module's own default shots, `n_bits=8` NEVER reaches chemical accuracy
  (>10 mHa on every one of 40 seeds) — the `__main__` demo's own smallest printed bit counts are
  honestly still in the gross-error regime. No library code changes — the checks live entirely in
  the test file, exercising `iqpe_ground_energy` unmodified. Gates G1–G4 in
  `tests/test_iqpe_cascade_risk_spec.py`.
  → [`SPEC_iqpe_cascade_risk.md`](SPEC_iqpe_cascade_risk.md) (exact-`p1` Bernoulli simulation, no
  device noise beyond that; one system, H2 CAS(2,2) — the shots_per_bit threshold is a measurement
  on this system's rotation-angle sequence, not a portable formula; no mitigation attempted).

- [x] **QPE readout — the precision law is a staircase, and the state-prep law is
  overlap-independent** — `qpe_walk_readout.py`'s docstring claims "resolution ~ lambda/2^t" and
  "success probability set by ground-state overlap," neither ever gated. **THE FINDING:** the
  literal "~lambda/2^t" smooth-scaling reading is FALSE — the realized argmax point-estimate error
  is a monotonically non-increasing STAIRCASE in `t` (adding phase bits never hurts, but the error
  often doesn't move for several consecutive `t` before a sudden drop), bounded above by a real,
  checked constant (measured max ratio `err/(lambda/2^t)` = 2.175, gated at 3x). The state-prep
  claim sharpens into a precise law: `p_success(window)/overlap` depends almost entirely on `t`, NOT
  on the overlap value — verified band width `< 0.05` across an overlap sweep spanning 0.02 to 0.99
  (a 50x range) at each of `t=8,10,12` (measured widths 0.029 / 0.002 / 0.0001). **Boundary, recorded
  not smoothed over:** that `t`-dependent prefactor is itself NOT monotonic in `t` — t=10's ratio
  (0.928) is measurably lower than t=8's (0.960), an unexplained artifact of window/dyadic-grid
  alignment, gated explicitly rather than glossed into a false "increases with t" story. No library
  code changes — the checks live entirely in the test file, exercising `qpe_walk_readout.py`'s
  existing `run_qpe` unmodified. Gates G1–G4 in `tests/test_qpe_readout_laws_spec.py`.
  → [`SPEC_qpe_readout_laws.md`](SPEC_qpe_readout_laws.md) (exact Fejer-kernel simulation oracle,
  not circuit-level; one system, H2 CAS(2,2) — the specific constants are measurements, not derived).

- [x] **krylov_subspace_solver — the documented bug fix gets a regression test, and the cross-check
  finally happens** — a second, independent real-time Krylov implementation (FCI-direct, not
  qubit-mapped) whose own docstring names its purpose as cross-checking the qubit-based path and
  documents a historical bug ("the collapse fix": an absolute overlap-eigenvalue cutoff that
  "dropped almost every vector and nulled the eigenproblem") — neither claim was ever gated.
  **THE FINDING:** the collapse fix is proven, not asserted — under a benign rescaling of the
  reference state (scale = 1, 1e-3, 1e-6, 1e-9), the RELATIVE cutoff keeps an IDENTICAL 7/12 basis
  vectors at every scale (scale-invariant, as a well-posed generalized eigenproblem requires), while
  an independently-reimplemented ABSOLUTE cutoff at the same numeric threshold collapses to **0/12
  kept vectors at scale <= 1e-6** — reproducing the exact historical failure mode. The cross-check
  the module was built for finally ran: its converged H4 CAS(4,4) energy agrees with
  `hybrid_quantum_solver.QuantumKrylovSolver`'s (an entirely independent, qubit-mapped statevector
  implementation) to **0.00037 mHa** — two genuinely separate code paths, sub-microhartree agreement.
  Variational floor pinned across closed- and open-shell systems (H2, H4, O2 triplet); the
  docstring's "condition numbers of 1e6+ are normal and harmless" claim measured (up to 9.5e6,
  error always < 3 mHa), not just asserted. No library code changes — G2's absolute-cutoff
  comparison is a deliberate reimplementation of the module's FORMER behavior in the test file, never
  imported from the module. Gates G1–G4 in `tests/test_krylov_subspace_solver_gates_spec.py`.
  → [`SPEC_krylov_subspace_solver_gates.md`](SPEC_krylov_subspace_solver_gates.md) (three systems,
  CAS(4,4) and smaller — the module's own documented `_dense_active_H` scope; the artificial
  rescaling in G2 demonstrates the mechanism on this system's real S matrix, not a claim that
  natural parameter ranges would themselves collapse without the fix — see R1).

- [x] **cross_check's own trust semantics — what "reference" means when CASCI can't run** — the
  "capstone for the near-term half of the stack" (four independent solvers — CASCI, Krylov, ADAPT,
  SQD — agreeing is what justifies trusting a number before spending QPU time on it) had an informal
  `__main__` assertion and an unexamined fallback: `reference = available.get("CASCI",
  next(iter(available.values())) ...)`. **THE FINDING:** with CASCI forced unreachable (cost-cap
  stress test, not a code change), the reference becomes EXACTLY the Krylov value — proving the
  fallback is a literal first-available pick in source insertion order, not a documented trust
  ranking; with CASCI AND Krylov both unreachable, it becomes ADAPT, never SQD — the accidental
  insertion order (CASCI, Krylov, ADAPT, SQD in the source) happens to keep the
  configuration-sampling method last in line, though nothing enforces this as a policy. **Boundary,
  recorded not smoothed over:** SQD has NO cost-cap knob and cannot be suppressed through the public
  API — with all three other caps forced to zero, SQD still ran unconditionally and became the sole
  reference; there is no way to drive `cross_check` into a fully-empty "no reference available"
  state through its public parameters alone, an asymmetry among the "four independent methods" worth
  knowing before assuming the caps make all four symmetric knobs. Baseline agreement (all four
  reachable) pinned too: H4 CAS(4,4) max deviation 0.0047 mHa, well inside the 5 mHa tolerance — and
  SQD's own deviation there (1.0e-9 mHa) was tighter than Krylov's or ADAPT's, a system-specific
  reminder that the accidental fallback order isn't a general accuracy ranking either. No library
  code changes — every gate drives `cross_check.py`'s existing `fci_max_dim`/`krylov_max_dim`/
  `qubit_dense_max_orb` cost caps, never touching internals. Gates G1–G4 in
  `tests/test_cross_check_trust_semantics_spec.py`.
  → [`SPEC_cross_check_trust_semantics.md`](SPEC_cross_check_trust_semantics.md) (one system, H4
  CAS(4,4), for the fallback-order gates; a mechanism check on the harness's own logic, not an
  accuracy claim about any individual solver — see R1).

- [x] **validate_and_cost — one threshold quietly governs two independent regime boundaries** — the
  three-stage capstone (taper → cross-check → FT cost) that composes every already-gated method in
  the near-term stack had never been gated itself. **THE FINDING:** `qubit_dense_max_orb`, one
  parameter forwarded through `validate_and_cost`, quietly controls TWO conceptually independent
  regime boundaries at once — this module's own taper-skip gate, and `cross_check`'s internal
  ADAPT-skip gate — with no documentation saying so. Proven on H4 CAS(4,4): at the default
  threshold both taper and ADAPT run; lowered below `norb`, both skip together, never
  independently. On a genuinely larger system (N2 CASCI(8,10), 8 orbitals) the pipeline degrades
  gracefully rather than crashing — CASCI/Krylov/SQD still ran and agreed (max deviation 0.0012
  mHa) despite taper and ADAPT both being out of regime. FT cost gracefully reports `None` in the
  standard `chem` env (no `openfermion`), confirmed rather than assumed. **Near-miss caught while
  probing, not shipped:** the first attempt to prove the threshold coupling RAISED
  `qubit_dense_max_orb` on the large N2 system instead of lowering it on the small one — this let
  `taper_hamiltonian` attempt `pauli_decompose` at 16 qubits (exponential in qubit count, per
  `SPEC_taper_spectrum.md`'s own R1) and OOM-killed the process. Redesigned to prove the identical
  finding in the safe direction (shrink an already-tractable system's cap) instead. No library code
  changes — every gate drives `validate_and_cost.py`'s existing public API. Gates G1–G4 in
  `tests/test_validate_and_cost_composition_spec.py`.
  → [`SPEC_validate_and_cost_composition.md`](SPEC_validate_and_cost_composition.md) (stage 3's FT
  cost numbers themselves are out of scope — needs the `chem-ft` env, outside `make gates`'s default
  flow; a mechanism check on the composition seam, not a fix for the coupling — see R1/R2).

## Killed

- [-] **Hₙ to larger n, *cheaply*** (ramp + D=100/200/400) — → [`SPEC_hchain_largen.md`](SPEC_hchain_largen.md).
  *Killed:* D=400 truncates too hard as chain entanglement grows — stderr balloons to ~5 mHa, the
  discarded-weight extrapolation falls back to `invD` by n=30, and leave-one-out = 1.07 mHa/atom
  (gate < 0.1). The ramp protocol is fine; cheap bond dims are not. Superseded by "done right" above.
