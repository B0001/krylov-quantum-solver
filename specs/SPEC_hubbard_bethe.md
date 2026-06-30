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
