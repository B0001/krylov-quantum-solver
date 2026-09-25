# SPEC: The solver reproduces the exact 1D Hubbard Bethe-ansatz ground-state energy

**Status:** CLOSED — gates G1–G4 PASS (2026-06-30); `hubbard_chain_integrals` + `lieb_wu_energy`
merged. TDL agreement with the Bethe-ansatz integral: U=2→1.6, U=4→5.5, U=8→1.4 mHa/site;
free-fermion + dimer limits machine-precision. Finding: HF-referenced Krylov converges on L=4 but
the strongly-correlated L≥6 Mott chain needs far deeper Krylov (recorded in G4).

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
  novelty. (b) **The TDL agreement is finite-L-FCI-limited:** with L ≤ 12 (FCI-tractable) and a
  1/L² extrapolation the residual is a few mHa/site (worst at intermediate coupling U≈4, ≈5–6
  mHa/site), not sub-mHa — tightening to < 1 mHa needs larger L via DMRG (a noted extension, not
  done here). (c) One-dimensional minimal lattice model.

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
  `lieb_wu_energy(U)` to `< 8e-3` Ha/site. **MEASURED:** U=2 → 1.6, U=4 → 5.5, U=8 → 1.4 mHa/site;
  intermediate coupling converges slowest (the recorded finite-size finding).
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

- Sub-mHa TDL agreement (needs DMRG at larger L — the Hn-TDL machinery would extend this).
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
