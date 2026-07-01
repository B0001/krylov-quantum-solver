# SPEC: Exact Nb₃X₈ cluster gaps expose where Hubbard-I breaks down (weakly-correlated Nb₃I₈)

**Status:** CLOSED — gates G1–G4 PASS (2026-06-30); `nb3x8_gaps.py` merged. Exact cluster charge gaps
(meV): Nb₃I₈ 842, Nb₃Br₈ 1086, Nb₃Cl₈ 1312, Nb₃F₈ 2581. Hubbard-I error grows monotonically as
U₀/|t| falls (0.2% → 0.8% → 5.2% → 29%); Hubbard-I underestimates the weakly-correlated Nb₃I₈ gap by
29% (244 meV). Both → U₀ as t→0 (validated). Scope: isolated cluster, impurity-solver error (not the
solid's gap).

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

The Nb₃X₈ family downfolds, per bilayer, to a generalized Hubbard dimer (two trimer orbitals; on-site
`U₀`, inter-layer hopping `t`, inter-site density-density `U_s⊥`) — a four-spin-orbital cluster the
source paper solves with the **Hubbard-I** approximation. That cluster is *exactly diagonalizable*.
Claim: exact diagonalization gives charge gaps the paper never reported, and the Hubbard-I error on
the same cluster **grows monotonically as correlation weakens (U₀/|t| falls)** — negligible for the
strongly-correlated Nb₃F₈/Nb₃Cl₈ but ≈ 29 % (≈ 244 meV) for the weakly-correlated Nb₃I₈ — while both
methods agree (→ U₀) in the atomic limit `t → 0`. The claim is false if the atomic-limit agreement
fails, or if the Hubbard-I error does not grow toward weak coupling / is small for Nb₃I₈.

## 2. Background and honest framing

- **Prior art / reference.** Aretz et al., *From strong to weak correlations in breathing-mode kagome
  van der Waals materials: Nb₃(F,Cl,Br,I)₈*, arXiv:2501.10320 — ab-initio cRPA downfolding + cluster
  DMFT (Hubbard-I). Parameters: their Table I (LT bulk). Uses the `model_hamiltonians.py` loader and
  the validated FCI reference already in this repo.
- **Ground truth.** Exact diagonalization (FCI) of the same cluster Hamiltonian.
- **What we can claim if gates pass.** Exact charge gaps for the Nb₃X₈ bilayer clusters from the
  paper's own parameters, and a validated, physically-consistent quantification of where the
  Hubbard-I (atomic self-energy) approximation breaks down across the correlation-tuned family.
- **What we cannot claim (stated up front, this is essential).** (a) **This is the *isolated*
  cluster.** The paper's cluster-DMFT embeds it in a self-consistent bath; the exact-vs-Hubbard-I
  comparison here measures the *impurity-solver error on the cluster*, **not** the solid's true gap.
  It is **not** a claim that the paper's material gaps are wrong by 29 %. (b) Density-density
  interactions only (the paper reports non-density-density terms of only a few meV). (c) A minimal
  two-orbital model; no phonons, no long-range tail beyond the strong inter-layer term.

## 3. Approach

Build the two-orbital cluster (`h1` = hopping, `eri` = on-site `U₀` + inter-site `U_s⊥`
density-density). Exact charge gap `Δ = E(3) + E(1) − 2E(2)` via `fixed_filling_energy`. Hubbard-I
gap from the atomic self-energy embedded in the bonding/anti-bonding dispersion (poles solve
`ω − Σ_at(ω) = ±t`). Reference: the exact gap; validation anchor: the `t → 0` atomic limit where
both equal `U₀`.

## 4. Public interface

```
nb3x8_gaps.NB3X8_LT_BULK                                   # {compound: {U0, t, Us}} (meV, Table I)
nb3x8_gaps.dimer_cluster_integrals(U0, t, Us) -> ModelIntegrals
nb3x8_gaps.exact_charge_gap(U0, t, Us) -> float           # meV
nb3x8_gaps.hubbard_i_gap(U0, t, Us) -> float              # meV
```

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_nb3x8_gaps_spec.py` (test-first). PySCF FCI; NumPy; no block2.

- **G1 — atomic-limit validation.** For every compound's `(U₀, U_s⊥)` at `t → 0`, both
  `exact_charge_gap` and `hubbard_i_gap` equal `U₀` to `< 1e-3` meV. (Confirms the machinery: the
  atomic Mott gap is `U₀`.)
- **G2 — exact gaps (the new numbers).** The exact charge gaps match the computed values to `< 1` meV:
  Nb₃I₈ 842, Nb₃Br₈ 1086, Nb₃Cl₈ 1312, Nb₃F₈ 2581 meV. All positive (insulating).
- **G3 — the finding (definition of done).** The Hubbard-I error `|Δ_HubI − Δ_exact| / Δ_exact` grows
  monotonically as `U₀/|t|` decreases across F → Cl → Br → I, is `< 1 %` for Nb₃F₈ and Nb₃Cl₈, and is
  `> 20 %` (measured 29 %) for the weakly-correlated Nb₃I₈. Hubbard-I **underestimates** the Nb₃I₈ gap
  (`Δ_HubI < Δ_exact`).
- **G4 — physical consistency.** The Hubbard-I error correlates with correlation strength: ordering
  the compounds by `U₀/|t|` orders them by `|error|` (Spearman = 1), and the exact gap is monotone in
  `U₀` across the family.

> Definition of done: **G3** (with G1 as the validation anchor). If a compound breaks the monotonic
> error trend, that is the finding — record it and check whether a non-density-density term or the
> bath matters there.

## 6. Implementation plan (test-first)

1. Write `tests/test_nb3x8_gaps_spec.py` encoding G1–G4 (initially failing — module absent).
2. Add `nb3x8_gaps.py` (cluster integrals; exact gap; Hubbard-I dimer gap; the parameter table).
3. Iterate to green via `make gates` (own process; pyscf/NumPy, no block2).

## 7. Out of scope

- Cluster DMFT with a self-consistent bath (the paper's actual method) and the solid's true gap.
- Non-density-density (Hund's/pair-hopping) terms and the long-range Coulomb tail.
- Monolayer/HT structures, doping, magnetism — additional data points, not the core finding.
- Spectral functions / dynamical quantities beyond the charge gap.

## 8. Caveats and risks

- **R1 — isolated cluster ≠ solid.** The single biggest caveat; stated in §2 and the module docstring,
  and the claim is scoped to the impurity-solver error, not the material gap.
- **R2 — Hubbard-I convention.** The atomic-self-energy embedding has sign/Hartree subtleties.
  *Mitigation:* the `t → 0` gate (G1) pins both methods to `U₀`, catching a mis-derivation loudly.
- Honest limitation: a minimal two-orbital density-density model; a methodological gap-error study,
  not a materials-prediction claim.

## 9. Deliverables

- `nb3x8_gaps.py` — `dimer_cluster_integrals`, `exact_charge_gap`, `hubbard_i_gap`, `NB3X8_LT_BULK`.
- `tests/test_nb3x8_gaps_spec.py` — gates G1–G4.
- Results summary (the exact-gap table + the monotonic Hubbard-I error, Nb₃I₈ ≈ 29 %, with the §2/§7
  caveats front and centre) in the PR description.
