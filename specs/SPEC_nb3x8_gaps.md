# SPEC: Exact Nb₃X₈ cluster gaps expose where Hubbard-I breaks down (weakly-correlated Nb₃I₈)

**Status:** CLOSED — gates G1–G4 PASS (2026-07-01); `nb3x8_gaps.py` merged. **Revised by the full
10-cluster dataset (LT bulk + LT bilayer + HT bulk).** Robust finding: strongly-correlated clusters
(Nb₃F₈, HT-phase Cl/Br) are near-exact (<2%); the iodides are consistently worst — Hubbard-I
underestimates Nb₃I₈ by 29% (bulk) / 12% (bilayer). **Recorded negative result:** the clean
single-parameter "error ∝ U₀/|t|" law from the 4-point LT-bulk subset does **not** survive the full
set (Spearman −0.86, non-monotone) — the error is multi-parameter (t and U_s⊥). Both → U₀ as t→0
(validated). Scope: isolated cluster (impurity-solver error, not the solid's gap).

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

The Nb₃X₈ family downfolds, per bilayer, to a generalized Hubbard dimer (two trimer orbitals; on-site
`U₀`, inter-layer hopping `t`, inter-site density-density `U_s⊥`) — a four-spin-orbital cluster the
source paper solves with the **Hubbard-I** approximation. That cluster is *exactly diagonalizable*.
Claim (tested across all 10 dimer-cluster parameter sets — LT bulk, LT bilayer, HT bulk): exact
diagonalization gives charge gaps the paper never reported, and Hubbard-I is near-exact (< 2 %) for
the strongly-correlated clusters (Nb₃F₈, HT Cl/Br) but **underestimates the weakly-correlated iodides
substantially** — ≈ 29 % (bulk) / 12 % (bilayer) — while both methods agree (→ U₀) at `t → 0`.
**Recorded negative result:** the tidy single-parameter "error ∝ U₀/|t|" law from the 4-point LT-bulk
subset does *not* survive the full dataset (Spearman ≈ −0.86, non-monotone) — the error depends on `t`
and `U_s⊥` together, so only the material-level statement (iodides worst) is robust. The claim is
false if the atomic-limit agreement fails, or if the iodides are not the worst / Hubbard-I is
small-error for Nb₃I₈.

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
nb3x8_gaps.NB3X8_CLUSTERS                                  # all 10 sets {name: {U0, t, Us}} (meV)
nb3x8_gaps.NB3X8_LT_BULK                                   # the 4 LT-bulk compounds (Table I)
nb3x8_gaps.dimer_cluster_integrals(U0, t, Us) -> ModelIntegrals
nb3x8_gaps.exact_charge_gap(U0, t, Us) -> float           # meV
nb3x8_gaps.hubbard_i_gap(U0, t, Us) -> float              # meV
```

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_nb3x8_gaps_spec.py` (test-first). PySCF FCI; NumPy; SciPy; no block2. Run over
all 10 dimer-cluster parameter sets (`NB3X8_CLUSTERS`).

- **G1 — atomic-limit validation.** For every cluster's `(U₀, U_s⊥)` at `t → 0`, both
  `exact_charge_gap` and `hubbard_i_gap` equal `U₀` to `< 1e-3` meV. (Confirms the machinery: the
  atomic Mott gap is `U₀`.)
- **G2 — exact gaps (the new numbers).** The exact charge gaps match the computed values to `< 1` meV
  for all 10 clusters (e.g. Nb₃I₈ bulk 842, bilayer 1961; Nb₃F₈ bulk 2581, bilayer 3979 meV). All
  positive (insulating).
- **G3 — the robust material-level finding (definition of done).** Strongly-correlated clusters
  (Nb₃F₈ LT, HT-phase Cl/Br) have `|error| < 2 %`; the iodides are the worst — Hubbard-I
  **underestimates** Nb₃I₈ bulk by `> 20 %` (≈ 29 %) and bilayer by `> 10 %` (≈ 12 %) — and Nb₃I₈-bulk
  is the single largest error over the whole set.
- **G4 — the recorded negative result.** The tidy single-parameter law does **not** hold across the
  full dataset: `Spearman(U₀/|t|, |error|)` is strongly negative but **not** −1
  (`−1 < ρ < −0.7`, measured −0.86), and `|error|` is **non-monotone** in `U₀/|t|` — the error depends
  on `t` and `U_s⊥` together, not a single ratio.

> Definition of done: **G3**. If a compound breaks the material-level finding (an iodide is *not* the
> worst, or a strongly-correlated cluster is *not* near-exact), that is the finding — record it and
> check whether a non-density-density term or the bath matters there. (The single-ratio law already
> broke on the extended dataset — G4 — which is why the robust claim is material-level, not a scaling
> law.)

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

- **R1 — isolated cluster ≠ solid (bounded, G5).** The biggest caveat, now *quantified* rather than
  just stated: enlarging the correlated region to include the inter-cluster weak link
  (`four_site_exact_gap`, a bath-fit-free proxy for the DMFT bath) moves the Nb₃I₈ gap by only ~5%
  (the strong intra-dimer bond isolates the cluster for the iodides) — far below the ~29% Hubbard-I
  error — and the 4-site Hubbard-I error *grows* to ~34%, so the finding survives cluster enlargement.
  The bath effect is largest for Nb₃F₈ (~22%, where `t_s ≈ t_w` makes the dimer ill-defined), but
  Hubbard-I is exact there anyway. Still not full cluster-DMFT (the paper's method); it bounds how far
  the finding travels, it does not compute the solid's gap.
- **R2 — Hubbard-I convention.** The atomic-self-energy embedding has sign/Hartree subtleties.
  *Mitigation:* the `t → 0` gate (G1) pins both methods to `U₀`, catching a mis-derivation loudly.
- Honest limitation: a minimal two-orbital density-density model; a methodological gap-error study,
  not a materials-prediction claim.

## 9. Deliverables

- `nb3x8_gaps.py` — `dimer_cluster_integrals`, `exact_charge_gap`, `hubbard_i_gap`,
  `four_site_exact_gap` (bath bound), `NB3X8_CLUSTERS`, `NB3X8_LT_BULK`, `NB3X8_LT_BULK_5P`.
- `tests/test_nb3x8_gaps_spec.py` — gates G1–G5 (G5 = the bath bound).
- Results summary (the 10-cluster exact-gap table, the robust iodides-worst finding, the falsified
  single-ratio law, and the bath bound, with the §2/§7 caveats front and centre) in the PR
  description — packaged to be sendable to the corresponding author.
