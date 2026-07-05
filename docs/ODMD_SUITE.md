# The ODMD suite — survival-amplitude methods, from measurement to material

Eleven spec-gated capabilities built on one primitive: the **survival amplitude**
`s_k = ⟨φ|e^(−ikτH)|φ⟩`, a 1-D time series the Krylov pipeline already measures as the first
row of its overlap matrix. Hankel-matrix dynamic mode decomposition (ODMD, `arXiv:2306.01858`)
turns that signal into eigenphases (energies) and Vandermonde amplitudes (spectral weights).
Everything below follows the repo's SDD loop — every claim has a spec, a falsifiable gate
against an exact reference (FCI / dense ED / analytic closed form), and a recorded boundary.
Run any of it: `make gates`, or each module's `__main__` demo.

## The ladder

| # | Spec | Claim (gated) | Headline number | Recorded boundary / finding |
|---|------|---------------|-----------------|------------------------------|
| 1 | [`SPEC_odmd`](../specs/SPEC_odmd.md) | FCI ground energy from the signal alone — no Hamiltonian matrix elements | beats KQD **57×** (10⁴ shots) / 10× (10⁵) at matched budget on N₂ CAS(6,6) | non-variational (dips below FCI at small K); the SVD truncation IS the robustness (~900× blow-up without) |
| 2 | [`SPEC_odmd_excited`](../specs/SPEC_odmd_excited.md) | the same signal carries the excited spectrum — if truncated at the random-matrix noise edge | H₄ gap to **5.8 mHa @ 10⁵ shots**, 21–31× below noisy QKSD | the **visibility law**: mode n recoverable iff pₙ√(dm) clears the noise edge — depth buys visibility as √K |
| 3 | [`SPEC_temple_bracket`](../specs/SPEC_temple_bracket.md) | every Krylov solve certifies itself from below (Temple + one ⟨H²⟩ matvec) | µHa-wide certified brackets; **zero containment escapes** over 4 systems × 8 depths | certification costs ~2.7× the raw error; the oracle-free mode is rigorous only at M ≥ 6 |
| 4 | [`SPEC_trotter_odmd`](../specs/SPEC_trotter_odmd.md) | circuit-real signals: ODMD = exact eigenphase of U_trot; Richardson kills the δt² bias | bias 12.0 → 0.0065 mHa (H₄, reps 2+4 extrapolation) | **found+fixed: the repo's Trotter path was silently exact** (opaque `PauliEvolutionGate`); extrapolate only when bias > noise |
| 5 | [`SPEC_device_odmd`](../specs/SPEC_device_odmd.md) | eigenphases are depolarizing-immune on real transpiled Hadamard-test circuits | **0.05 mHa through 70% amplitude loss** (Aer, cx=3e-4); KQD fails ~2600× on matched data | immunity ends at the shot floor, exactly where the visibility law predicts |
| 6 | [`SPEC_nb3x8_device_gap`](../specs/SPEC_nb3x8_device_gap.md) | a material's charge gap through the whole stack (3 sectors × device ODMD + Richardson) | Nb₃I₈ cluster gap **842.44 meV measured to 1.2%** on the noisy simulated device | sector Trotter biases do NOT cancel in gaps (12% at reps=1); bias-vs-noise crossover gated |
| 7 | [`SPEC_odmd_uq`](../specs/SPEC_odmd_uq.md) | error bars from ONE signal, validated by coverage | union CI: coverage **0.895–1.000** at nominal 90%, all systems × budgets | parametric bootstrap alone is anti-conservative 18×; no resampling sees model bias (K-too-small ⇒ coverage 0) |
| 8 | [`SPEC_odmd_spectral`](../specs/SPEC_odmd_spectral.md) | photoemission A(ω): poles = ionization lines, amplitudes = Lehmann weights | Nb₃I₈ Hubbard-band gap = **842.44 meV to the last digit** (cross-pin) | HF-referenced weights ≈ 2% off true Lehmann; intensities are damping-immune too |
| 9 | [`SPEC_odmd_optical`](../specs/SPEC_odmd_optical.md) | optical absorption + exciton binding (P-kick) | exciton binding collapses **0.986·Us (F) → 0.263·Us (I)** — the exciton unbinds with hopping | **found+fixed: eigenstate kicks hung `expm_multiply`** (~10¹⁷ substeps); one bright line per dimer (selection rule gated) |
| 10 | [`SPEC_odmd_spin`](../specs/SPEC_odmd_spin.md) | spin spectroscopy: the interlayer exchange J | **J = 0.051 / 66.2 / 119.1 / 245.9 meV** (F/Cl/Br/I) | Heisenberg 4t²/(U₀−Us) fails 0→46.5% across the family; local moment 1.000→0.759 |
| 11 | [`SPEC_visibility_law`](../specs/SPEC_visibility_law.md) | the visibility law is predictive: shots* ∝ 1/(w²K), one calibration transfers | slope **−2.000** over 4 orders in w; one Br calibration predicts every line to ~10% | attribution tolerance must sit below the line spacing (strong lines masquerade 45× otherwise); tight attribution costs 61% at K=8, 4% at K=32 |

## The Nb₃X₈ dimer scorecard (LT bulk, meV — cRPA parameters of `arXiv:2501.10320`)

Numbers the source paper did not report, all exact on the isolated cluster (its charge gaps
were the subject of [`SPEC_nb3x8_gaps`](../specs/SPEC_nb3x8_gaps.md), which also shows why
isolated-cluster values are **upper bounds** on the broadened solid — quote accordingly):

| | Nb₃F₈ | Nb₃Cl₈ | Nb₃Br₈ | Nb₃I₈ |
|---|---|---|---|---|
| charge gap Δc | 2580.8 | 1311.8 | 1086.0 | 842.4 |
| optical gap Δopt | 1876.0 | 1117.5 | 963.7 | 774.4 |
| exciton binding Δc − Δopt | 704.9 (0.99·Us) | 194.3 (0.49·Us) | 122.3 (0.36·Us) | 68.0 (0.26·Us) |
| interlayer exchange J | 0.051 (≈0)¹ | 66.2 | 119.1 | 245.9 |
| Heisenberg-J error | 0.0% | 6.3% | 14.1% | 46.5% |
| local moment ‖Sᶻψ₀‖² | 1.000 | 0.944 | 0.890 | 0.759 |
| oscillator weight ‖Pψ₀‖² | 1.1e-4 | 0.22 | 0.44 | 0.96 |

¹ below the model's own neglected non-density-density terms.

## Running it on real hardware

`run_hardware_odmd.py` takes the suite to an actual device (or a device-noise Aer simulation).
It measures the H₂ survival signal by ancilla Hadamard tests (reps 1+2), prints a transpiled
resource table, then runs the full stack — `odmd_energy` → Richardson (SPEC 4) → single-signal
error bar (SPEC 7) — against the circuit eigenphase and FCI. `--backend aer` works now;
`--backend <ibm_name>` submits via qiskit-ibm-runtime (needs a saved account — the script prints
the one-line setup if none is found); `--dry-run` stops at the resource table.

**Hardware-readiness finding (from `--dry-run`):** the ancilla-controlled construction is deep —
H₂ at K=8 transpiles to ~3.4k two-qubit gates (reps 1) / ~6.6k (reps 2), beyond today's NISQ
reach. In Aer the K=8 pipeline still recovers the circuit eigenphase to <0.01 mHa and Richardson
lands within ~1 mHa of FCI, so the *algorithm* is correct — **circuit depth, not the method, is
the bottleneck**. Levers: smaller K (depth ~ controlled-U^{K−1}), reps=1-only, or a shallower
signal scheme. The dry-run resource table is how you price a device run before spending queue time.

## Defects found (and permanently gated) along the way

1. **The silently-exact Trotter path** — `Operator()`/`Statevector.evolve()` ignore an opaque
   `PauliEvolutionGate`'s synthesis; `TrotterKrylovSolver`'s `reps`/`order` were no-ops and its
   "within Trotter error" test passed vacuously. Fixed by materializing the synthesis;
   regression gate `SPEC_trotter_odmd` G1.
2. **The anti-conservative bootstrap** — the DMD fit absorbs realized noise, so a naive
   parametric bootstrap under-covers by up to 18×. Fixed by the union with a bagging arm;
   coverage-gated in `SPEC_odmd_uq`.
3. **The eigenstate-kick hang** — a symmetry-forced single-line reference gives reachable width
   0 and a π·10¹²-long evolution. Fixed by the degenerate-reference short-circuit; gated in
   `SPEC_odmd_optical` G2.

## Scope, in one paragraph

Validation-scale systems (dense-diagonalizable; H₂/H₄/LiH/N₂ CAS(6,6) and the Nb₃X₈ dimers);
idealized i.i.d. shot noise except where an Aer device model is stated (and that is still a
simulated device — no coherent errors, crosstalk, or drift); ODMD estimates are never
variational bounds (use the Temple bracket when a bound is needed); spectra are
reference-weighted (HF vs exact-ψ₀ differences are gated, not hidden); all Nb₃X₈ numbers are
isolated-cluster values. The measured advantages are matched-budget comparisons against this
repo's own validated implementations, at this scale — not hardware claims.
