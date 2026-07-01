# Exact charge gaps of the Nb₃X₈ bilayer cluster, and where Hubbard-I departs from them

*A short computational note on the downfolded models of Aretz, Grytsiuk, Liu, …, van Loon & Rösner,
"From strong to weak correlations in breathing-mode kagome van der Waals materials: Nb₃(F,Cl,Br,I)₈",
arXiv:2501.10320.*

## Summary

The paper's per-bilayer impurity problem — two inter-layer-dimerized trimer molecular orbitals with
on-site `U₀`, strong inter-layer hopping `t_s⊥`, and inter-site density-density `U_s⊥` — is a
four-spin-orbital cluster and can be diagonalized exactly. Using the paper's own cRPA parameters
(Table I, LT bulk & bilayer; Table IV, HT bulk), I report the **exact charge gaps** of this cluster
and compare them to the **Hubbard-I** approximation the paper's cluster-DMFT employs. The exact
gaps are new numbers; the comparison quantifies where Hubbard-I departs from the exact cluster
result, and a cluster-enlargement test bounds the isolated-cluster approximation.

**This is a statement about the impurity solver on the cluster, not about the solid's gap** (see
Caveats). It most likely confirms, in numbers, the paper's own reason for using CTHYB rather than
Hubbard-I away from the atomic limit.

## Exact vs Hubbard-I charge gap (meV), 10 dimer-cluster parameter sets

Charge gap `Δ = E(N+1) + E(N−1) − 2E(N)` at half-filling (N=2), by exact diagonalization; Hubbard-I
gap from the atomic self-energy embedded in the bonding/anti-bonding dimer dispersion. Both reduce to
`Δ = U₀` in the atomic limit `t_s⊥ → 0` (a validation of the machinery).

| set        | U₀/\|t_s⊥\| | Δ exact | Δ Hubbard-I | Hubbard-I error |
|------------|-----------:|--------:|------------:|----------------:|
| Nb₃I₈  LT-bulk  |   3.6 |   842 |   599 | **−29 %** |
| Nb₃Br₈ LT-bulk  |   7.0 |  1086 |  1029 |   −5.2 % |
| Nb₃I₈  LT-bil   |   8.8 |  1961 |  1723 | **−12 %** |
| Nb₃Cl₈ LT-bulk  |  10.7 |  1312 |  1322 |   +0.8 % |
| Nb₃Br₈ LT-bil   |  14.2 |  2281 |  2233 |   −2.1 % |
| Nb₃Cl₈ LT-bil   |  19.8 |  2550 |  2565 |   +0.6 % |
| Nb₃Br₈ HT-bulk  |  54.9 |  1092 |  1109 |   +1.5 % |
| Nb₃Cl₈ HT-bulk  |  81.9 |  1369 |  1384 |   +1.1 % |
| Nb₃F₈  LT-bulk  |   529 |  2581 |  2586 |   +0.2 % |
| Nb₃F₈  LT-bil   |   798 |  3979 |  3984 |   +0.1 % |

## Findings

1. **Hubbard-I is near-exact for the strongly-correlated clusters** (Nb₃F₈, HT-phase Nb₃Cl₈/Nb₃Br₈:
   |error| < 2 %) and **substantially underestimates the weakly-correlated iodides** — Nb₃I₈ by
   ~29 % (bulk) and ~12 % (bilayer), i.e. ~240 meV on the bulk gap. The iodides carry the largest
   inter-layer hopping (`|t_s⊥|` ≈ 218 meV) against the smallest `U₀` (787 meV), so hybridization
   competes with `U`, exactly where an atomic self-energy must fail. This is consistent with the
   paper's use of CTHYB rather than Hubbard-I for the harder/doped cases.

2. **The error is not a single-parameter function of `U₀/|t_s⊥|`.** Across the full set the rank
   correlation between `U₀/|t_s⊥|` and |error| is only −0.86 (not −1), and |error| is non-monotone
   (the Nb₃I₈ bilayer sits out of order). Both `t_s⊥` and `U_s⊥` matter; the robust statement is
   material-level (iodides worst), not a scaling law.

3. **The isolated-cluster approximation is bounded and favourable for the iodides.** Enlarging the
   correlated region to two dimers joined by the weak inter-bilayer link (`t_w⊥`, `U_w⊥` from
   Table I) shifts the *exact* Nb₃I₈ gap by only ~5 % — far below the 29 % Hubbard-I error — because
   the strong intra-dimer bond isolates the cluster (`|t_s⊥|` ≫ `|t_w⊥|`) for the iodides. On the
   enlarged cluster the Hubbard-I error *grows* (to ~34 %) rather than washing out. The cluster
   enlargement matters most for Nb₃F₈ (~22 % gap shift, where `t_s⊥ ≈ t_w⊥` and the "dimer" is barely
   defined) — but Hubbard-I is exact there anyway. So the finding is most reliable exactly where it
   is largest.

## Caveats (stated plainly)

- **Isolated cluster, not the solid.** This compares an exact solver to Hubbard-I *on the isolated
  cluster*; it does **not** reproduce the paper's self-consistent cluster-DMFT, and it is **not** a
  claim that the reported material gaps are wrong by 29 %. Finding 3 bounds how far the cluster result
  travels toward the solid; it does not compute the solid's gap.
- **Density-density interactions only**, consistent with the paper's statement that non-density-density
  (Hund's / pair-hopping) terms are of order a few meV.
- **Minimal two-orbital model**; no phonons, no long-range Coulomb tail beyond the strong inter-layer
  term, no monolayer (which has no dimer).
- The direction (Hubbard-I worse toward weak correlation) is expected; the contribution is the exact
  cluster gaps and the quantified, bath-bounded size of the departure.

## Method / reproducibility

Two-orbital cluster: `h1 = [[0, t_s⊥],[t_s⊥, 0]]`; density-density ERIs `U₀` on-site, `U_s⊥`
inter-site. Exact charge gap via full CI in the fixed particle-number sectors. Hubbard-I: atomic
self-energy `Σ(ω) = ω − 2x(x−U₀)/(2x−U₀)`, `x = ω − h` (`h` the inter-site Hartree shift), embedded in
the dimer dispersion; gap = spacing between the occupied and unoccupied poles. Bath bound: the same
construction on a four-site chain (two dimers + weak link). Parameters taken verbatim from Table I
(LT) and Table IV (HT) of arXiv:2501.10320.

Code and a test-gated derivation (including the `t → 0 → U₀` validation and the bath-bound check) are
in `nb3x8_gaps.py` / `tests/test_nb3x8_gaps_spec.py` of the accompanying repository; the numbers above
are regenerated by `python nb3x8_gaps.py`.
