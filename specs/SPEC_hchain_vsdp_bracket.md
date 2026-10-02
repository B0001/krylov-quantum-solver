# SPEC: Certified two-sided H-chain bracket — rigorous 2-RDM SDP lower bound + variational upper bound

**Status:** DONE (chem-3z8). PQG bracket **killed** at n=12 (2.49 mHa/atom > 1). PQG+T2 passes where it
fits on this machine (n ≤ 8, ≤ 0.05 mHa/atom) but cannot be evaluated at n=12 here — see §8.
chem-7ko evaluated it with SU(2) × reflection adapted T2: **0.126 mHa/atom at n = 12**, so the claim
survives, given the singlet premise ([SPEC_hchain_t2sym](SPEC_hchain_t2sym.md)).

---

## 1. Goal

For open H_n chains (STO-6G, Löwdin site orbitals, R = 1.8 bohr, the `benchmark_hchain_tdl.integrals`
geometry) produce a **floating-point-rigorous** lower bound on E0 and a variational upper bound, and
tabulate the bracket width per atom vs n. Pre-registered kill: width > 1 mHa/atom at n = 12.

## 2. Background and honest framing

- Chaykin, Jansson, Keil et al. (JCTC 2020) made 2-RDM SDP lower bounds rigorous (VSDP) for small
  molecules. This is that method applied to H chains beyond easy FCI, not a new method.
- **Can claim:** each tabulated lower bound is a mathematical lower bound on the ground-state energy of
  the Hamiltonian defined by the float64 integrals, independent of solver accuracy (the solver only
  supplies a dual vector; any vector gives a valid, possibly weak, bound).
- **Cannot claim:** rigor with respect to the integrals themselves (PySCF float64 is taken as the
  definition, as in VSDP); anything about larger basis sets; T2-tightness beyond n = 8.

## 3. Approach

- Blocks are Gram matrices `M[a,b] = <O_a† O_b>` (T2: `<{O_a†, O_b}>`) split by Sz change:
  D1, Q1 (per spin); D2, Q2 (αα, ββ, αβ); G (ΔSz = 0, ±1); optional T2 (four Sz groups). Every entry
  is reduced to 1-/2-RDM parameters by a generic normal-ordering routine; no P/Q/G/T2 formula is
  hand-typed. Constraints: entry definitions, spin-resolved contractions
  `Σ_j <a†_iσ a†_jτ a_jτ a_kσ> = (N_τ − δ_στ) <a†_iσ a_kσ>`, `tr D1σ = Nσ`. Only Sz = 0 with fixed
  (Nα, Nβ) is assumed; no S² or point-group constraint is imposed.
- **Rigorous bound** for any dual y: with Z_j = C_j − A_jᵀy,
  `E0 ≥ bᵀy + Σ_j tr(X*_j) · min(0, λ_min(Z_j))`, where tr(X*_j) is the exact (Nα, Nβ)-only trace
  of the true RDM block and λ_min(Z_j) is bounded below by a floating-point Cholesky of Z_j − σI
  (Higham, *Accuracy and Stability*, Thm 10.3: R̃ᵀR̃ = M + ΔM, |ΔM| ≤ γ_{n+1}|R̃ᵀ||R̃|) plus explicit
  entrywise rounding bounds on Z_j and on the objective coefficients. No directed rounding is needed.
- **Solver:** SCS 3 (eps 1e-8) is the untrusted producer; Clarabel agrees to 1e-7 Ha at n = 8 but costs
  ~(block size)^6, so it is a cross-check only.
- **Upper bounds:** exact FCI (`dmrg_reference.fci_energy`) for n ≤ 12; for n = 14, 16 the final
  two-site-sweep energy of SU(2) DMRG at D = 1500 — variational, **not** the D→∞ extrapolation.

## 4. Public interface

```
vsdp_hchain.build(h1, eri, na, nb, t2=False) -> dict      # standard-form SDP data
vsdp_hchain.solve(prob, solver="CLARABEL", **kw) -> (value, y, status)
vsdp_hchain.rigorous_lower_bound(prob, y, e_const) -> (bound, per-block λ bounds)
vsdp_hchain.certify(prob, y, e_const) -> best of the two dual-sign conventions
uv run python vsdp_hchain.py --ns 2,...,16 [--t2]  -> data/hchain_vsdp[_t2].json
uv run python vsdp_hchain.py --dmrg-upper N        -> data/hchain_vsdp_upper.json (block2 isolated)
```

## 5. Acceptance criteria (validation gates) — `tests/test_hchain_vsdp_bracket_spec.py`

- **G1 — linear maps vs brute force.** At the exact H4 ground state (Fock-space diagonalization, no
  PySCF FCI), every block entry, contraction, objective and trace identity holds to 1e-10, with and
  without T2. Kills any sign/convention error in the normal ordering.
- **G2 — N = 2 is exact.** H2 certified bound within 1e-6 Ha below E0 (D ⪰ 0 is exact for N = 2).
- **G3 — valid and tight.** n = 4, 6: certified bound ≤ FCI and within 1e-5 Ha of the solver value.
- **G4 — rigor under adversarial duals.** Zero, random, perturbed and sign-flipped y all give ≤ E0.
- **G5 — T2 tightens.** PQG bound < PQG+T2 bound ≤ E0 at n = 4.
- **Kill test (the finding, not a gate):** width per atom at n = 12.

## 6. Results

PQG, certified (every row: certified bound ≤ solver value by ≤ 0.001 mHa):

| n | certified lower bound (Ha) | upper bound (Ha) | upper method | width (mHa/atom) |
|---|---|---|---|---|
| 2 | −1.11873391 | −1.11873390 | FCI | 0.0000 |
| 4 | −2.19237877 | −2.19038422 | FCI | 0.499 |
| 6 | −3.27426596 | −3.26674310 | FCI | 1.254 |
| 8 | −4.35944020 | −4.34507940 | FCI | 1.795 |
| 10 | −5.44629315 | −5.42438538 | FCI | 2.191 |
| 12 | −6.53412465 | −6.50422696 | FCI | **2.491** |
| 14 | −7.62256267 | −7.58439082 | DMRG D=1500 | 2.727 |
| 16 | −8.71142884 | −8.66476140 | DMRG D=1500 | 2.917 |

PQG+T2, certified:

| n | width (mHa/atom) | largest T2 block |
|---|---|---|
| 4 | 0.0006 | 88 |
| 6 | 0.022 | 306 |
| 8 | 0.050 | 736 |
| 10 | 0.088 | 500 (SU(2) × reflection adapted; Sz-only 1450) |
| 12 | **0.126** | 864 (adapted; Sz-only 2520) |

The n = 10 and 12 rows are chem-7ko's ([SPEC_hchain_t2sym](SPEC_hchain_t2sym.md) §6). They are valid
for a nondegenerate singlet ground state, which converged FCI verifies at both n.

n = 8 took 91 min / 2.4 GB with SCS; its certified bound sits 0.026 mHa below the solver value
(the T2 dual is less accurate than PQG's). Only the two small T2 Sz-blocks (a†_s a_q a_p with p, q
both of one spin and s of the other) at n = 6: 0.445 mHa/atom — they recover about
two-thirds of the PQG→T2 tightening, not enough to argue n = 12 from.

**Verdict.** The PQG bracket is killed: 2.49 mHa/atom at n = 12, and the width per atom grows
monotonically with n (no sign of saturating below 1). T2 is what makes the bracket tight. At n = 12 its
two large Sz blocks are 2520 × 2520 each, out of reach for SCS on a 16 GB laptop with Sz-only
blocking. So the pre-registered PQG+T2 claim is **not killed and not confirmed** at n = 12.

**Update (chem-7ko).** Adapting T2 to SU(2) and chain reflection makes n = 12 fit. The certified
PQG+T2 widths are 0.088 mHa/atom at n = 10 and **0.126 mHa/atom at n = 12**, well under the 1
mHa/atom kill, so the claim **survives at n = 12**. This relies on the ground state being a
nondegenerate singlet, which converged FCI verifies (⟨S²⟩ ≈ 1e-14, gap 0.109 Ha at n = 12). Beyond
n = 12 the premise is unverified.

## 7. Out of scope

- Spin adaptation (S² = 0, spin-flip) and chain-reflection block reduction; both need the ground state
  to be a nondegenerate singlet (true by FCI for n ≤ 12) and are what n = 12 T2 needs.
- Basis sets beyond STO-6G; other geometries.

## 8. Caveats and risks

- **R1 — BLAS error model.** The Cholesky bound assumes LAPACK `potrf` obeys the standard γ_{n+1}
  inner-product error model (true for blocked/FMA implementations; Strassen-type products are not
  used by Accelerate/OpenBLAS potrf). This is the same assumption VSDP and Rump's `isspd` make.
- **R2 — Upper bounds at n = 14, 16** are variational DMRG energies, i.e. upper bounds only up to
  Davidson tolerance (1e-9 Ha), negligible at mHa resolution. At n = 10 the same path gives FCI + 3e-10.
- **R3 — T2 at n = 12** needs ≈ 4× smaller blocks (spin + reflection adaptation) or a custom
  boundary-point solver (as in Mazziotti's and DePrince's v2RDM codes). **Resolved by chem-7ko.**
  Adaptation cut the largest block 2520 → 864, which is 2.9×, not 4×. That was enough: SCS solved
  n = 12 in 1.37e4 s and 3.24 GB ([SPEC_hchain_t2sym](SPEC_hchain_t2sym.md) §6.3).

## 9. Deliverables

- `vsdp_hchain.py` — SDP builder, rigorous dual bound, SCS/Clarabel driver, isolated DMRG upper bound.
- `tests/test_hchain_vsdp_bracket_spec.py` — G1–G5.
- `cvxpy` added as a dependency (brings SCS and Clarabel).
