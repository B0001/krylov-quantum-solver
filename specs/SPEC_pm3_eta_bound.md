# SPEC: A sound Krylov-moment bound on the HF weight η missed by an overshooting self-mode floor (chem-pm3)

**Status:** PRE-REGISTERED 2026-10-02. The kill rule and gates below were fixed and committed
*before* any witness number was computed with the sound LP. Results go in §10.
**Depends on:** `SPEC_subspace_floor_resolvability.md` (the escape witnesses), `hf_overlap_subspace.py`,
`hybrid_quantum_solver.certified_overlap`, `reachability._dense_hf_projection` (the reference).

---

## 1. Goal

On the three guard-passing escape witnesses of `SPEC_subspace_floor_resolvability.md` (self-mode
floor β above the true E_d): **(G-a)** the leakage-corrected certificate
γ_corr = √max(0, 1 − r²/δ² − η) is a valid lower bound on ‖P_S u‖ for *any* β; **(G-b)** the Krylov
moments S₀ₖ = ⟨HF|e^{−ikΔtH}|HF⟩, k < M ≤ 24, alone bound η to **≤ 0.05**. G-a should hold by
construction; G-b is the claim that can die (the bead's author predicted it would).

## 2. Background and honest framing

- **G-a derivation.** S′ = {eigenstates with E < β} is a cluster β floors *by definition*, so the
  sin-θ bound gives 1 − r²/δ² ≤ ‖P_S′ u‖² unconditionally (δ = β − λ_u). S ⊂ S′, hence
  ‖P_S u‖² = ‖P_S′ u‖² − η with η the HF weight on levels in [E_d, β).
- **G-b.** The S₀ₖ are the trigonometric moments of the HF spectral measure μ (already the first
  row of the solver's overlap matrix: no extra measurements). For any real trig polynomial f of
  degree M − 1 with f ≥ 1_[E_d, β] on supp μ: η ≤ μ([E_d, β]) ≤ ∫f dμ = Σₖ fₖ S̄₀ₖ. Minimizing over
  f is a semi-infinite LP (a Markov–Krein-type extremal problem).
- **Prior state (the bug this spec fixes).** The scout LP enforced f ≥ 1_band only at grid points.
  The optimum then dips below the indicator *between* grid points, at the spectral atoms, so it
  reports too little: 0.1472 < η = 0.15 on the synthetic gate, and 0.14998 at a 10× finer grid.
  The scout numbers (and its "2 of 3 witnesses survive") are therefore void.
- **Can claim if the gates pass:** whether noise-free Krylov moments, given oracle band edges,
  carry enough information to bound η to 0.05 on these witnesses at M ≤ 24.
- **Cannot claim:** a certified sub-cluster anchor. The band edge E_d and the support
  [e_min, e_max] are oracle inputs (dense diagonalization), so a **survival only means "not ruled
  out"**; a **kill is robust** (any valid oracle-free band or support is wider, which can only
  raise the bound). Three H₆/STO-3G witnesses only; noise-free moments.

## 3. Approach

- **Reference:** exact dense diagonalization via `reachability._dense_hf_projection` (ZHEEVR
  fallback for the macOS ZHEEVD bug at 12 qubits): reachable levels (HF weight > 1e-8), η, E_d,
  ‖P_S u‖. Floor β = θ_d − σ_d from the Krylov Ritz pairs, exactly as the library's self mode.
- **Sound LP (cutting planes).** Solve on a grid; compute the *continuous* minimum of
  g = f − 1_[E_d, β] on [e_min, e_max] (domain and band endpoints plus every stationary point of f:
  the angles of all roots of z^{M−1} f′, Newton-polished); add the violated minimizers as
  constraints; re-solve until min g ≥ −tol. Return ∫f dμ + max(0, −min g)·S₀₀: f plus that gap is
  ≥ 1_band everywhere, so every returned value is a valid bound even if the loop stops early.
- **Box** |coefficients| ≤ cap (default 1000) keeps HiGHS stable; it only shrinks the feasible
  set, so it can raise, never invalidate, the bound.

## 4. Public interface

```
scripts/spec_pm3_subspace_eta_bound.py
  _build_moment_bound_lp(moments, e_min, e_max, dt, band_lo, band_hi,
                         n_grid=12000, coeff_cap=1000.0, ...) -> float | None
  run_witness(w) -> dict      # G-a, eta, bounds by M, verdict, box sensitivity
  main()                      # the 3-witness run behind §10
```

## 5. Acceptance criteria — pre-registered

- **G-a.** On each witness, γ_corr ≤ exact ‖P_S u‖ + 1e-9 (plus the algebraic identity on 20
  random synthetic spectra).
- **G-b soundness (a bug check, not the claim).** Synthetic 7-atom measure with closed-form
  moments: bound ≥ η − 1e-6 for M ∈ {4, 8, 12, 20} and cap ∈ {1, 10, …, 1e4}; bound non-increasing
  in M and in cap (1e-6); the stationary-point minimum matches a brute-force 10⁶-point minimum;
  and on the three real witnesses B(M) ≥ η at every M.
- **G-b kill rule (the claim).** B_w(M) at the default box (cap 1000) for M ∈ {8, 12, 16, 20, 24}.
  Witness w **SURVIVES iff min over those M of B_w(M) ≤ 0.05**, else **KILLED**. (The feasible sets
  are nested in M, so for an exact solve the min equals the M = 24 value; the min is the literal
  reading of the bead's "at M ≤ 24".) Overall: **full kill** iff all three are killed, otherwise
  the per-witness verdicts stand ("survives on k of 3").
- **Box rule.** B_w(24) is also reported at cap 1e4 and 1e5. A kill that a larger cap flips to
  ≤ 0.05 is recorded as *box-dependent*, not as a kill.
- **Descriptive, not gates.** For each killed witness: the smallest M ∈ {28, 32, …, 64} with
  B_w(M) ≤ 0.05 (the quantified resolution requirement), 1/(MΔt), and the band-edge gaps
  E_d − E_{d−1} and (next reachable level) − β. The bead's design target "1/(MΔt) below the
  overshoot" counts as **contradicted** if a witness is killed at an M with 1/(MΔt) < overshoot.

Done = G-a and soundness green, and the per-witness verdicts recorded in §10 and pinned by
`tests/test_pm3_eta_bound_spec.py`.

## 6. Implementation plan

1. Gates in `tests/test_pm3_eta_bound_spec.py` (the two synthetic soundness gates were red on main).
2. Make the LP sound (cutting planes + certified gap); route dense diagonalization through
   `_dense_hf_projection`; fold the two scout check scripts into the driver.
3. Run `main()`, record §10, correct the BACKLOG entry.

## 7. Out of scope

- Oracle-free band edges, i.e. the certified sub-cluster anchor itself.
- Shot noise on S₀ₖ: the bound moves by at most ‖x‖₁·max|ΔS₀ₖ| ≤ (2M − 1)·cap·max|ΔS₀ₖ|, so the
  box is also a noise-amplification budget. Not measured here.
- Changing the library's self mode (`hf_overlap_subspace.py` stays heuristic).

## 8. Caveats and risks

- **Oracle inputs** (§2): E_d and [e_min, e_max] are exact; Δt = π/width is the solver default.
- **Floating-point rigor:** stationary points come from `np.roots` + Newton polish; gated against
  brute force. Moments come from the solver's propagated basis, not closed form.
- **Box:** a binding cap means B(cap) ≥ B(∞); the box rule above keeps a cap-made kill out.
- **Reachability cut** 1e-8 defines "levels"; η counts only reachable weight (the LP bounds all of
  μ's band weight, which is ≥ η).

## 9. Deliverables

- `scripts/spec_pm3_subspace_eta_bound.py` — sound LP + 3-witness driver.
- `tests/test_pm3_eta_bound_spec.py` — gates.
- §10 results; the chem-pm3 entry in `specs/BACKLOG.md`.

## 10. Results

*(Filled after the pre-registered run.)*
