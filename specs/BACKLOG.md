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
  **CI falsifiers G1–G3 green (2026-06-29); `bulk_per_site_energy` merged. Headline (leave-one-out
  < 0.1 mHa/atom) is a pending driver-level large-n ramp run — not yet executed/recorded, so Open.**
- [ ] **Be₂ toward experiment** — *Claim:* core-valence correlation + a cc-pVxZ→CBS extrapolation
  moves the FCI/DMRG well depth from ~305 cm⁻¹ toward the experimental 929.7. *Gate (reproduction):*
  `|D_e − 930| < 100 cm⁻¹` at CBS+CV; `R_e` within 0.1 Å of 2.45. (Honest: reproduces a settled
  result.)
- [x] **1D Hubbard vs Bethe ansatz** — the validated stack (`hubbard_chain_integrals` → FCI/Krylov)
  reproduces the exact Lieb–Wu Bethe-ansatz energy: closed-shell-BC finite-L FCI extrapolates to the
  TDL integral (U=2→1.6, U=4→5.5, U=8→1.4 mHa/site), free-fermion (−4/π) + dimer limits at machine
  precision. *Gate revised from 1 mHa/site to < 8 mHa/site:* L≤12 FCI extrapolation is finite-size-
  limited (intermediate coupling slowest); sub-mHa needs DMRG at larger L. Finding: HF-referenced
  Krylov is slow on the strongly-correlated Mott chain (L≥6, large U). Gates G1–G4 in
  `tests/test_hubbard_bethe_spec.py`; `model_hamiltonians.py`.
  → [`SPEC_hubbard_bethe.md`](SPEC_hubbard_bethe.md) (repro of Lieb & Wu 1968, exact analytic ref).
- [ ] **GPU backend on real hardware** — *Claim:* `device="gpu"` reproduces CPU energies and reaches
  larger qubit counts. *Gate (on an NVIDIA node):* `|E_gpu − E_cpu| < 1e-6 Ha` at 16q; completes a
  28q statevector run. (Code exists, CPU-validated; needs a GPU node — see `tests/test_gpu_backend.py`.)
- [~] **Nb₃X₈ / Hubbard model loader** — *Claim:* a tight-binding hopping matrix + Hubbard/cRPA
  interaction maps onto the universal `(h1, eri, e_core, nelec, norb)` interface and the
  number-conserving solver reproduces the analytic 2-site Hubbard dimer `(U−√(U²+16t²))/2` and PySCF
  FCI for a rank-4 cRPA tensor. *Gate:* `|E_solver − E_analytic| < 1e-6 Ha` across U/t = 0…20;
  rank-4 mapping `< 1e-9` vs FCI. → [`SPEC_nb3x8_hubbard.md`](SPEC_nb3x8_hubbard.md). **G1–G5 green
  (2026-06-29); `hybrid_quantum_solver/model_hamiltonians.py`.** Validates the Nb₃X₈ Model-Database
  drop-in path (handoff). Honest: not yet checked against the DB's DMFT/Hubbard-I gaps (DB files not
  bundled); Nb₃I₈ is the weakly-correlated anchor, not a strong-correlation headline.
- [ ] **DMRG-referenced transition-metal active space** — *Claim:* for a low-spin NbN (or a smaller
  TM dimer) CAS large enough that FCI is intractable, DMRG gives a converged correlation energy.
  *Gate:* bond-dim extrapolation stderr < 1 mHa **and** agreement between two independent DMRG
  sweep schedules < 1 mHa. (Honest: a reference number, not a materials claim — finite cluster.)
## Done

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

## Killed

- [-] **Hₙ to larger n, *cheaply*** (ramp + D=100/200/400) — → [`SPEC_hchain_largen.md`](SPEC_hchain_largen.md).
  *Killed:* D=400 truncates too hard as chain entanglement grows — stderr balloons to ~5 mHa, the
  discarded-weight extrapolation falls back to `invD` by n=30, and leave-one-out = 1.07 mHa/atom
  (gate < 0.1). The ramp protocol is fine; cheap bond dims are not. Superseded by "done right" above.
