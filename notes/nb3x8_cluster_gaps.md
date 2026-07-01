# Exact charge gaps of the Nb₃X₈ bilayer cluster, band broadening, and the reliability of Hubbard-I

*A short computational note on the downfolded models of Aretz, Grytsiuk, Liu, …, van Loon & Rösner,
"From strong to weak correlations in breathing-mode kagome van der Waals materials: Nb₃(F,Cl,Br,I)₈",
arXiv:2501.10320.*

## Summary

The paper's per-bilayer impurity problem — two inter-layer-dimerized trimer orbitals with on-site
`U₀`, strong inter-layer hopping `t_s⊥`, and inter-site density-density `U_s⊥` — is a four-spin-orbital
cluster and can be diagonalized exactly. I report the **exact charge gaps** of this cluster from the
paper's cRPA parameters, and I use exact diagonalization of *growing* clusters + DMRG to check how the
gap evolves toward the solid. **The main, self-correcting result:** on the *isolated* cluster
Hubbard-I underestimates the weakly-correlated iodide gaps by ~30 %, but that gap **is an artifact of
neglecting band broadening** — restoring inter-dimer coordination (which the paper's cluster-DMFT
keeps) collapses the exact gap back toward the Hubbard-I value. So the exercise ends up **supporting
the paper's Hubbard-I / cluster-DMFT gaps**, not challenging them. The exact isolated-cluster gaps may
still be useful as a reference.

## 1. Exact vs Hubbard-I on the isolated cluster (meV)

Charge gap `Δ = E(N+1) + E(N−1) − 2E(N)` at half-filling, by exact diagonalization; Hubbard-I gap from
the atomic self-energy in the bonding/anti-bonding dimer dispersion. Both reduce to `Δ = U₀` at
`t_s⊥ → 0` (validation).

| set             | U₀/\|t_s⊥\| | Δ exact | Δ Hubbard-I | error |
|-----------------|-----------:|--------:|------------:|------:|
| Nb₃I₈  LT-bulk  |   3.6 |  842 |  599 | −29 % |
| Nb₃Br₈ LT-bulk  |   7.0 | 1086 | 1029 | −5.2 % |
| Nb₃I₈  LT-bil   |   8.8 | 1961 | 1723 | −12 % |
| Nb₃Cl₈ LT-bulk  |  10.7 | 1312 | 1322 | +0.8 % |
| Nb₃Br₈ LT-bil   |  14.2 | 2281 | 2233 | −2.1 % |
| Nb₃Cl₈ LT-bil   |  19.8 | 2550 | 2565 | +0.6 % |
| Nb₃Br₈ HT-bulk  |  54.9 | 1092 | 1109 | +1.5 % |
| Nb₃Cl₈ HT-bulk  |  81.9 | 1369 | 1384 | +1.1 % |
| Nb₃F₈  LT-bulk  |   529 | 2581 | 2586 | +0.2 % |
| Nb₃F₈  LT-bil   |   798 | 3979 | 3984 | +0.1 % |

On the isolated cluster, Hubbard-I is near-exact (<2 %) for the strongly-correlated members and
underestimates the weakly-correlated iodides (up to ~29 %). The size of the error is not a clean
function of `U₀/|t_s⊥|` (rank correlation −0.86, non-monotone): both `t_s⊥` and `U_s⊥` matter.

## 2. Band broadening: the isolated-cluster error does not translate to the solid

The isolated-dimer gap contains **no inter-dimer band broadening**, so it is an *upper bound* on the
solid gap; Hubbard-I embedded in the full dispersion (cluster-DMFT) includes broadening and lies
lower. They bracket the true gap from opposite sides. Restoring coordination for Nb₃I₈:

- **Coordination scan** (central dimer + `z` out-of-plane weak-link neighbours, exact FCI):
  `z=0` 842 → `z=1` 797 → `z=2` 747 → `z=3` (≈ the paper's √3-of-3 weak neighbours) **650 meV**.
- **1-D SSH chain → thermodynamic limit** (DMRG/block2, chain length 8→20 sites, monotone
  730 → 709): **~708 meV** for the out-of-plane stacking alone.
- With realistic 3-D coordination (out-of-plane `z≈3` plus in-plane `t_∥`), the exact solid gap is
  **~600–650 meV — close to the Hubbard-I / cluster-DMFT value (~599)**.

So the ~29 % isolated-cluster discrepancy is largely an artifact of the isolated-cluster
approximation; once broadening is restored, exact diagonalization and Hubbard-I converge. (A
one-nearest-neighbour "bath bound" gives only a ~5 % shift and is *misleading* — it under-samples the
coordination.)

## 3. Conclusion

- **Reference numbers:** exact charge gaps for the Nb₃X₈ bilayer clusters (Section 1) — possibly a
  useful cross-check for the downfolded models.
- **Method assessment:** for the Nb₃X₈ solid gaps, Hubbard-I (in cluster-DMFT, with the full
  dispersion) appears **robust** — the exact solid-gap estimate lands near it. This is consistent with
  the paper's own use of Hubbard-I at integer filling and its switch to CTHYB for the harder/doped
  cases.
- Net: a null / confirming result. The interesting turn was methodological — the isolated-cluster gap
  looked like a 29 % Hubbard-I failure until band broadening was put back, at which point it
  dissolved.

## Caveats

- Two-orbital density-density model (the paper reports non-density-density terms of a few meV); no
  phonons or long-range Coulomb tail beyond the strong inter-layer term.
- The coordination / 1-D-chain estimates are proxies for the full 3-D lattice, not a self-consistent
  cluster-DMFT; they bound the solid gap and its trend, they do not reproduce the paper's DMFT.
- DMRG cross-checked against FCI where FCI is tractable (and it corrected an FCI convergence failure at
  the L=12 half-filled chain).

## Method / reproducibility

Two-orbital cluster: `h1 = [[0, t_s⊥],[t_s⊥, 0]]`; density-density ERIs `U₀` on-site, `U_s⊥`
inter-site; exact charge gap by full CI in the fixed particle-number sectors; Hubbard-I from the
atomic self-energy embedded in the dimer dispersion. Coordination scan and 1-D SSH chain: the same
construction on larger clusters (FCI up to ~L=12; DMRG/block2 to L=20). Parameters verbatim from
Table I (LT) and Table IV (HT) of arXiv:2501.10320. Code and test-gated derivation (including the
`t → 0 → U₀` validation and the coordination-collapse check) are in `nb3x8_gaps.py` /
`tests/test_nb3x8_gaps_spec.py`; the numbers regenerate with `python nb3x8_gaps.py`.
