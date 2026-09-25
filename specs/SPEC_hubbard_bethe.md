# SPEC: The solver reproduces the exact 1D Hubbard Bethe-ansatz ground-state energy

**Status:** CLOSED — gates G1–G4 PASS (2026-06-30); `hubbard_chain_integrals` + `lieb_wu_energy`
merged. **Extended with DMRG at large L (bead chem-tjr, 2026-09-25, §10).** Using open chains with
L ≤ 100 and a bulk difference quotient, the TDL agreement with the Bethe-ansatz integral is
**U=2 → +0.024, U=4 → +0.018, U=8 → −0.004 mHa/site**. All three are below 1 mHa/site; the FCI
L ≤ 12 values were 1.6 / 5.5 / 1.4. The ring cross-check agrees within its bars at U=4 and U=8 but
**fails at U=2** (+0.41 mHa/site, 1/L² fit biased at L ≤ 32). One U=8 point failed its ladder check
twice (§10.4). Free-fermion and dimer limits are at machine precision. Finding: HF-referenced Krylov
converges on L=4, but the strongly correlated L≥6 Mott chain needs far deeper Krylov (recorded in G4).

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

The half-filled 1D Hubbard model has an exact analytic ground-state energy per site (Lieb & Wu 1968,
a Bethe-ansatz integral). Claim: the project's validated stack — `hubbard_chain_integrals` →
`build_hamiltonian_from_integrals` / PySCF FCI / `QuantumKrylovSolver` — reproduces it. Concretely,
finite-size FCI per-site energies (with a closed-shell boundary phase that removes the even/odd-L
shell zigzag) extrapolate to the Lieb–Wu thermodynamic-limit integral across coupling regimes, and
the exactly-solvable limits (U=0 free fermions, the 2-site dimer) match to machine precision. The
claim is false if the extrapolated energy disagrees with the integral beyond the finite-size
tolerance, or if an analytic limit is missed.

## 2. Background and honest framing

- **Prior art / reference.** Lieb & Wu, Phys. Rev. Lett. 20, 1445 (1968) — the exact 1D Hubbard
  Bethe-ansatz solution. A *clean analytic reference*, the strongest kind of falsifier (no DMRG, no
  fit to data). Uses the new `model_hamiltonians.py` lattice loader.
- **What we can claim if gates pass.** The number-conserving solver / FCI on the 1D Hubbard chain
  reproduces the exact Bethe-ansatz energy: machine-precision at the free-fermion and dimer limits,
  and a few-mHa/site finite-size agreement with the thermodynamic-limit integral across U.
- **What we cannot claim (stated up front).** (a) Reproduction of a settled exact result, not
  novelty. (b) **The FCI TDL agreement is finite-L-limited:** with L ≤ 12 and a 1/L² extrapolation,
  the residual is a few mHa/site, worst at U≈4 (5.5 mHa/site). *Update (§10): DMRG at L ≤ 100
  removes this limit.* Open-chain difference quotients reach |e − e_LW| ≤ 0.024 mHa/site for
  U ∈ {2, 4, 8}, read from the vendored table by `tests/test_hubbard_lieb_wu_spec.py`. The
  sub-mHa claim rests on route (a). The independent ring route reaches it too (≤ 0.41 mHa/site),
  but at U=2 it disagrees with route (a) beyond both self-assessed bars. So the *bars* are not
  validated at weak coupling, even though the energies are. (c) One-dimensional minimal lattice model.

## 3. Approach

`hubbard_chain_integrals(L, U, t)` builds the half-filled ring with the closed-shell boundary phase
(periodic if L/2 odd, antiperiodic if even) so the per-site energy converges smoothly. FCI in the
fixed (L/2, L/2) sector (`fixed_filling_energy`) gives e(L); a least-squares fit in 1/L² extrapolates
to e_∞. Reference: `lieb_wu_energy(U)` (the Bethe-ansatz TDL integral, by quadrature) for the TDL,
the analytic free-fermion sum for U=0, and `hubbard_dimer_energy` for L=2.

## 4. Public interface

```
model_hamiltonians.hubbard_chain_integrals(n_sites, U, t=1.0, closed_shell=True, ...) -> ModelIntegrals
model_hamiltonians.lieb_wu_energy(U, t=1.0) -> float       # exact TDL energy per site (Bethe ansatz)
```

(`fixed_filling_energy`, `hubbard_dimer_energy`, `QuantumKrylovSolver` already exist.)

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_hubbard_bethe_spec.py` (test-first). PySCF FCI (no block2); L ≤ 12 keeps the
gate fast. t = 1.

- **G1 — TDL agreement with the Bethe-ansatz integral (definition of done).** For U ∈ {2, 4, 8},
  the FCI per-site energies at L ∈ {6, 8, 10, 12} (closed-shell BC) extrapolated in 1/L² agree with
  `lieb_wu_energy(U)` to `< 8e-3` Ha/site. **MEASURED (FCI, L ≤ 12):** U=2 → 1.6, U=4 → 5.5,
  U=8 → 1.4 mHa/site; intermediate coupling converges slowest (the recorded finite-size finding).
  **MEASURED (DMRG, §10, vendored `specs/hubbard_lieb_wu_table.csv`, gate G6 of
  `test_hubbard_lieb_wu_spec.py`):** open chains, bulk quotient over L = 80, 100:
  U=2 → +0.024, U=4 → +0.018, U=8 → −0.004 mHa/site (< 1e-3 Ha/site; the gate asserts it).
- **G2 — free-fermion (U=0) limit, machine precision.** `lieb_wu_energy(0) = -4/π` to `< 1e-4`, and
  the U=0 chain FCI equals the analytic closed-shell free-fermion energy (twice the sum of the L/2
  lowest hopping eigenvalues) to `< 1e-8` Ha.
- **G3 — dimer limit.** The L=2 chain FCI equals `hubbard_dimer_energy(t, U)` to `< 1e-8` for several
  U (the U→∞ superexchange end is covered here).
- **G4 — solver reproduces FCI on the lattice.** `QuantumKrylovSolver` on the chain (fixed filling)
  matches FCI to `< 1e-6` Ha for L=4 at U ∈ {2, 4, 8} (Krylov depth 24). **Finding:** real-time
  Krylov from |HF⟩ converges on the small chain but the strongly-correlated half-filled Mott chain
  is hard — |HF⟩ has poor overlap with the true ground state, so L ≥ 6 at large U needs far deeper
  Krylov (L=6, U=8 is still ≈160 mHa off at depth 24). Gated where it converges; the limitation is
  recorded, not hidden.

> Definition of done: **G1** (with G2/G3 as the machine-precision analytic anchors). If the 1/L²
> extrapolation cannot reach the tolerance at L ≤ 12 (intermediate coupling), the residual is the
> finite-size finding — record it and note the DMRG-to-larger-L tightening, do not fake sub-mHa.

## 6. Implementation plan (test-first)

1. Write `tests/test_hubbard_bethe_spec.py` encoding G1–G4 (initially failing).
2. Add `hubbard_chain_integrals` + `lieb_wu_energy` to `model_hamiltonians.py`.
3. Iterate to green via `make gates` (own process; pyscf/qiskit, no block2).

## 7. Out of scope

- ~~Sub-mHa TDL agreement~~ — done with DMRG in §10 (bead chem-tjr).
- Away-from-half-filling, finite temperature, magnetic fields, the full finite-L Bethe-equation
  solver (the TDL integral is the reference used here).
- 2D Hubbard (no exact Bethe solution).

## 8. Caveats and risks

- **R1 — open-shell zigzag.** Plain periodic BC makes e(L) oscillate with L/2 parity, ruining the
  fit. *Mitigation:* the closed-shell boundary phase (built into `hubbard_chain_integrals`); G2's
  free-fermion check guards the BC.
- **R2 — intermediate-coupling finite-size error.** U≈4 converges slowest (≈5–6 mHa at L≤12).
  *Mitigation:* gate at `< 8e-3` Ha/site and record the per-U convergence; the residual shrinks with
  L_max.
- Honest limitation: minimal 1D lattice, reproduces a settled exact result.

## 9. Deliverables

- `model_hamiltonians.py` — `hubbard_chain_integrals`, `lieb_wu_energy`.
- `tests/test_hubbard_bethe_spec.py` — gates G1–G4.
- Results summary (per-U TDL agreement + the analytic-limit checks, with §2/§7 caveats) in the PR.

---

## 10. Extension: DMRG at large L (bead chem-tjr)

Caveat (b) above says the few-mHa TDL residual is finite-L-FCI-limited. This section tests that by
composing `hubbard_chain_integrals` with `dmrg_energy_extrapolated` (block2, SU(2)) at half filling
far past FCI. Driver `benchmark_hubbard_lieb_wu.py`; analysis `hubbard_tdl_analysis.py`; vendored
table `specs/hubbard_lieb_wu_table.csv`; gate `tests/test_hubbard_lieb_wu_spec.py`.

### 10.1 Budget probes and the stall they exposed (before pre-registration)

Hardware: `Linux 6.18.44-fc-v37 x86_64`, Intel Xeon @ 2.80GHz, nproc = 4, 15 GiB RAM, no GPU.

- **Stall, labelled `converged`.** Open L=60, U=4, default block2 random MPS, D = 100/200/400, 8
  sweeps per stage: E = −20.780 / −30.025 / −34.044 Ha (the D=100 stage sits ~13 Ha high), every
  discarded weight ≤ 4e-9, **`regime = "converged"`**, fit −38.666 Ha, 590 s. The chem-4e9 stall
  checks do not fire: the slope check runs on truncation ladders only, and the undershoot check
  allows 10× the last gap (4.0 Ha) while the fit sat 4.6 Ha below E(D_max). **Cause:** the
  default random MPS spreads the charge unevenly; at the Mott charge gap DMRG moves it ~1 Ha per
  sweep. Open L=40, U=4, D=100 for 12 sweeps gives −21.648 Ha from the default start (still falling
  by ~0.6 Ha/sweep) and −22.583594 Ha from a uniform-filling start (`occs = 1`, converged in 4
  sweeps). This repeats with the generic qc MPO and with a hand-built Hubbard MPO, so the Hamiltonian
  is not the cause.
- **Fix (library).** `dmrg_energy_extrapolated(..., init_occs=...)` passes block2's `occs` through.
  `ExtrapResult.stage_dE` also records each perD stage's last-sweep |ΔE|, a stall detector that
  does not depend on `regime`. Both are off by default, so existing callers are unchanged.
- **Probes with `occs = 1`, U=4** (seen before pre-registration, and not in the table):
  open L=40, D=100/200/400: E(400) = −22.5835938881, 30 s; open L=100: E(400) = −57.0053055086,
  discarded weight 2.5e-12, stage ΔE ≤ 1e-7, 173 s. Ring L=10: E(D=200) equals PySCF FCI to
  2e-11. Ring L=16 (folded order), D=400: discarded weight 3.6e-8, 30 s.

### 10.2 Pre-registered protocol and analysis (fixed before any production number)

The analysis code, `hubbard_tdl_analysis.py`, is committed with this section and must not change
after production numbers exist.

- **Points.** Every ladder uses `protocol="perD"`, 8 sweeps per stage, and `init_occs = 1`.
  Route (a) uses **open** chains, L ∈ {20, 40, 60, 80, 100}. Route (b) uses **closed-shell rings**
  (`hubbard_chain_integrals` default boundary phase, sites in folded order 0, L−1, 1, L−2, …),
  L ∈ {16, 20, 24, 28, 32}. Both routes run U ∈ {2, 4, 8}, t = 1. Bond-dimension ladders are in
  §10.3; they are chosen from probe **wall times and discarded weights only**.
- **Point check** (all must hold; `point_check`): regime ∈ {converged, truncation};
  E(D_max) − E_extrap ≤ the ladder's own spread (max − min of E(D)); E_extrap ≤ E(D_max) + 1 µHa;
  every stage's last-sweep |ΔE| ≤ 1 µHa. A point's uncertainty is
  σ = max(stderr, |E_extrap − E(D_max)|). A point that fails is re-run once with a larger ladder,
  and the re-run is recorded. The analysis does not change.
- **Route (a):** e_a = `bulk_per_site_energy` over L = 80 and 100.
  bar_a = |e_a − q(60, 80)| + (σ₁₀₀ + σ₈₀)/20.
- **Route (b):** e_b = intercept of a least-squares fit e(L) = e_b + a/L² over all ring L.
  bar_b = |e_b − (the same fit without the smallest L)| + intercept stderr + max_L σ(L)/L.
- **Acceptance (bead text, unchanged):** (1) |e_a − `lieb_wu_energy`(U)| < 1e-3 Ha/site for
  U ∈ {2, 4, 8}; route (a) is the headline and U=4 is the headline coupling.
  (2) |e_a − e_b| < bar_a + bar_b. (3) G2 and G3 above still pass. (4) G1's numbers and caveat (b)
  are updated from the vendored table, read by the spec gate.

### 10.3 Bond-dimension ladders (from probe cost and weights only)

The ring L=32, U=4 probe used folded order, `occs = 1`, D = 300/600/1200 and 4 threads. It
finished in 623 s, with discarded weights 1.7e-6 / 1.3e-7 / 7.0e-9 and every stage's ΔE ≤ 5e-7.
That places it in the truncation regime, and its σ of about 2e-6 Ha is about 1e-7 Ha/site.
U=2 has a smaller charge gap and more entanglement, so both ladders get headroom above the probes:
**open D = 200/400/800** and **ring D = 400/800/1600**, with 4 threads and the points run one at a
time. The production order is U=4 (open, then ring), U=2, then U=8.

### 10.4 Results (2026-09-25)

All numbers come from `specs/hubbard_lieb_wu_table.csv` through the pre-registered
`hubbard_tdl_analysis.py`; `tests/test_hubbard_lieb_wu_spec.py` (G5–G8) re-derives and pins them.
Hardware is as in §10.1. The 30 table points plus 1 rejected point took **13 500 s of wall time
(≈ 3.75 h)**, run one at a time on 4 threads. The slowest points were the U=8 L=100 re-run (2456 s)
and ring U=2 L=32 (1111 s); open L=100 took 484–570 s.

| U | e_LW (Ha/site) | (a) open, e_a ± bar_a | resid (a) | (b) ring, e_b ± bar_b | resid (b) | (2): \|e_a−e_b\| vs bar_a+bar_b |
|---|---|---|---|---|---|---|
| 2 | −0.84437434 | −0.84435007 ± 1.6e-5 | **+0.024 mHa** | −0.84396542 ± 8.5e-5 | +0.409 mHa | 3.85e-4 vs 1.01e-4 → **FAIL** |
| 4 | −0.57372937 | −0.57371181 ± 1.1e-5 | **+0.018 mHa** | −0.57357087 ± 1.4e-4 | +0.159 mHa | 1.41e-4 vs 1.48e-4 → pass |
| 8 | −0.32753053 | −0.32753404 ± 3.4e-5 | **−0.004 mHa** | −0.32752684 ± 2.4e-6 | +0.004 mHa | 7.2e-6 vs 3.7e-5 → pass |

- **Acceptance (1): PASS for every U.** The U=4 headline residual is **+0.018 mHa/site**, compared
  with 5.5 mHa/site from FCI at L ≤ 12, a ×300 improvement and 55× inside the 1 mHa target.
- **Acceptance (2): pass at U=4 (narrowly) and U=8; FAIL at U=2.** This is recorded, not tuned
  away. DMRG does not cause it: every ring σ is ≤ 1.2e-5 Ha, or ≤ 4e-7 Ha/site. A *post-hoc*
  diagnostic (not part of the pre-registered analysis) supports finite-size curvature. The
  quantity L²·(e_ring(L) − e_LW) should be flat when the 1/L² form holds. At U=4 it is flattening
  (−0.736, −0.683, −0.660, −0.651, −0.647 for L = 16…32). At U=2 it still drifts linearly
  (−1.824, −1.748, −1.675, −1.604, −1.536). Why: at U=2 the charge gap is small (Δ_c ≈ 0.17 t),
  so the charge correlation length is ≳ 10 sites, and rings with L ≤ 32 sit in the crossover
  where 1/L² is not yet the leading form. The drop-smallest-L term of bar_b cannot see a smooth
  drift. **Limiting factor for route (b): L** (the ring length relative to ξ_c), not D.
- **The bars are not conservative (finding).** Measured against the exact answer, which the
  bars did not use:
  - Route (a) residuals are 1.5× bar_a at U=2 and 1.6× at U=4. The successive-quotient difference
    captures only part of the O(1/(L₁L₂)) remainder.
  - Route (b) residuals are 1.2× (U=4), 4.8× (U=2) and 1.5× (U=8) bar_b.

  The energies are sub-mHa; the self-assessed uncertainties are too small by ≈ 1.5–5×. The
  analysis is left unchanged, as pre-registered. A future spec should widen the bars by
  construction, e.g. with a 1/L⁴ or exponential term in the fit envelope.
- **One point fails its ladder check twice: open U=8 L=100.**
  - First run (D = 200/400/800, `specs/hubbard_lieb_wu_rejected.csv`): the D=200 stage's final
    sweep still moved 1.07e-2 Ha.
  - Pre-registered re-run (D = 400/800/1600): the D=400 stage stalled 0.55 mHa high (final-sweep
    ΔE 5.6e-2). `regime` read `uncontrolled` and the 1/D fit fell 0.28 mHa below E(D_max).
  - E(800) and E(1600) agree with each other, and with the first run's E(800), to 1e-13 Ha.
    Only the first, un-annealed stage stalls: even from the uniform start, 8 sweeps are not
    always enough at U=8, L=100.
  - The pre-registered analysis uses the re-run row as recorded. Divided by the 20-site step, its
    0.28 mHa error costs 0.014 mHa/site. Substituting E(D_max) (sensitivity only) gives
    +0.010 mHa/site instead of −0.004. **Acceptance (1) at U=8 holds either way.** Limiting
    factor: first-stage sweeps, not D or L. G5 pins this as the only failed point.
- **Anchors.** DMRG through the production path matches the U=0 free-fermion energy, open and
  folded ring, L=20 (G8). The ring has a real truncation error of 6e-5 Ha at D=400, and the exact
  energy lies inside the pre-registered σ. DMRG also matches the dimer to 1e-8. Probe: ring L=10,
  U=4 E(D=200) equals PySCF FCI to 2e-11. The existing G2–G4 still pass (acceptance 3).
- **Regime-classifier gap (follow-up, not fixed here).** `truncation_regime` labels the §10.1
  stalled ladder `converged` (G7 pins it). The chem-4e9 checks are skipped or too loose on the
  1/D branch. This spec rejects stalls with `stage_dE`. The classifier itself should learn
  from it: either take stage convergence as an input, or bound the undershoot on a converged
  ladder by `ENERGY_NOISE` alone.
