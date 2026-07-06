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

## 3. Beyond the charge gap: optical, excitonic, and magnetic constants of the same clusters

The same exactly-solvable cluster fixes three more observables the paper did not report, two of
them in closed form. On the inversion-symmetric dimer the polarization `P = n₁ − n₂` is odd and
the singlet sector contains exactly **one** odd state (the ionic-odd combination, at energy
exactly `U₀`), so the entire optical spectrum is a single line and

```
Δ_opt = U₀ − E₀ ,          E₀ = (U₀+U_s⊥)/2 − √( ((U₀−U_s⊥)/2)² + 4 t_s⊥² )
J     = √( ((U₀−U_s⊥)/2)² + 4 t_s⊥² ) − (U₀−U_s⊥)/2        (singlet–triplet splitting)
```

with the charge gaps of Section 1 this gives the **exciton binding** `E_b = Δ_c − Δ_opt` and the
**inter-layer exchange** `J` for every parameter set (meV; `⟨(S₁ᶻ−S₂ᶻ)²⟩` is the local-moment
fraction, `‖P|ψ₀⟩‖²` the total oscillator weight of the one bright line):

| set             | U₀/\|t_s⊥\| | Δ_opt | E_b = Δ_c − Δ_opt | E_b/U_s⊥ | J | 4t²/(U₀−U_s⊥) err | moment | osc. wt |
|-----------------|-----------:|------:|------:|------:|--------:|-------:|------:|--------:|
| Nb₃I₈  LT-bulk  |   3.6 |  774 |   68 | 0.26 | 245.9 | +46.5 % | 0.759 | 0.96 |
| Nb₃Br₈ LT-bulk  |   7.0 |  964 |  122 | 0.36 | 119.1 | +14.1 % | 0.890 | 0.44 |
| Nb₃I₈  LT-bil   |   8.8 |  814 | 1147 | 0.85 | 234.4 | +40.4 % | 0.776 | 0.89 |
| Nb₃Cl₈ LT-bulk  |  10.7 | 1118 |  194 | 0.49 |  66.2 |  +6.3 % | 0.944 | 0.22 |
| Nb₃Br₈ LT-bil   |  14.2 | 1025 | 1256 | 0.85 | 111.7 | +12.2 % | 0.902 | 0.39 |
| Nb₃Cl₈ LT-bil   |  19.8 | 1189 | 1361 | 0.87 |  62.4 |  +5.5 % | 0.950 | 0.20 |
| Nb₃Br₈ HT-bulk  |  54.9 |  855 |  237 | 0.86 |   2.0 |  +0.2 % | 0.998 | 9.2e-3 |
| Nb₃Cl₈ HT-bulk  |  81.9 | 1065 |  304 | 0.90 |   1.1 |  +0.1 % | 0.999 | 4.1e-3 |
| Nb₃F₈  LT-bulk  |   529 | 1876 |  705 | 0.99 |  0.05 |  +0.0 % | 1.000 | 1.1e-4 |
| Nb₃F₈  LT-bil   |   798 | 2001 | 1978 | 0.99 |  0.05 |  +0.0 % | 1.000 | 1.0e-4 |

Observations (cluster-level, same standing as Section 1):

- **The exciton unbinds with hopping in the LT-bulk series** — `E_b/U_s⊥` falls 0.99 → 0.26 from
  F to I — but is *not* a single-ratio law either: the **bilayers** keep `E_b/U_s⊥ ≈ 0.85` even
  for the iodide, because their inter-site `U_s⊥` is much larger. The bilayer excitons are
  eV-scale (1.1–2.0 eV).
- **The inter-layer exchange is large for the iodides** (`J ≈ 235–246 meV`) and the perturbative
  superexchange formula `4t²/(U₀−U_s⊥)` overestimates it by **40–47 %** there — the iodide dimer
  is well beyond the Heisenberg regime — while it is essentially exact for the HT phases and the
  fluorides. Correspondingly the local moment is reduced by ~24 % in the iodides.
- **A selection-rule statement that is exact for the cluster:** all optical weight sits in one
  line (the odd singlet); the triplet is dipole-dark and appears only in the spin channel. The
  oscillator weight itself spans four orders of magnitude across the family — the iodide dimer
  is ~10⁴× more polarizable than the fluoride.

## 4. Finite-temperature magnetic susceptibility χ(T) of the same clusters

The same N=2 spectrum gives the magnetic susceptibility by the Van Vleck trace
`χ(T) = ⟨S_z,tot²⟩_thermal / T` (reduced units, k_B=1, meV; g=2 for the emu conversion
`χ[emu/mol] = 0.12931·χ_reduced[meV⁻¹]`). This is textbook magnetochemistry — the isotropic
S=½ dimer is the **Bleaney–Bowers** system `χ_BB(T) = (2/T)/(3+e^{J/T})` (Bleaney & Bowers, Proc.
R. Soc. A 214, 451, 1952), with the finite-T deviations being the exact two-site **Hubbard-dimer**
thermodynamics (Carrascal et al., arXiv:1502.05038; Anderson, Phys. Rev. 115, 2, 1959). It is
**not** a new phenomenon, and χ(T) of Nb₃Cl₈ is already measured (Sheckelton et al., Inorg. Chem.
Front. 4, 481, 2017, incl. the ~90 K singlet transition; Haraguchi et al., Inorg. Chem. 56, 3483,
2017) and modeled (Grytsiuk/Katsnelson/van Loon/Rösner, arXiv:2305.04854). What the exact cluster
adds is the **family-wide table from the ab-initio-downfolded (t, U₀, U_s⊥)**, and a clean
statement of *where the pure-spin picture breaks*.

| set             | J (meV) | E_s (meV) | E_s/J | θ_CW (meV) | χ(300K) emu/mol | μ_eff(300K) | T₅%/E_s | T₅%/J |
|-----------------|--------:|----------:|------:|-----------:|----------------:|------------:|--------:|------:|
| Nb₃I₈  LT-bulk  |   245.9 |   774.4 |   3.1 | −61.5 | 7.4e-7 | 0.021 | 0.407 |   1.28 |
| Nb₃Br₈ LT-bulk  |   119.1 |   963.7 |   8.1 | −29.8 | 9.7e-5 | 0.241 | 0.430 |   3.48 |
| Nb₃I₈  LT-bil   |   234.4 |   814.1 |   3.5 | −58.6 | 1.2e-6 | 0.026 | 0.411 |   1.43 |
| Nb₃Cl₈ LT-bulk  |    66.2 |  1117.5 |  16.9 | −16.6 | 6.3e-4 | 0.613 | 0.438 |   7.38 |
| Nb₃Br₈ LT-bil   |   111.7 |  1025.4 |   9.2 | −27.9 | 1.3e-4 | 0.277 | 0.432 |   3.96 |
| Nb₃Cl₈ LT-bil   |    62.4 |  1189.1 |  19.1 | −15.6 | 7.1e-4 | 0.651 | 0.438 |   8.35 |
| Nb₃Br₈ HT-bulk  |     2.0 |   854.6 |   432 |  −0.50 | 2.5e-3 | 1.213 | 0.444 |    192 |
| Nb₃Cl₈ HT-bulk  |     1.1 |  1065.3 |   969 |  −0.28 | 2.5e-3 | 1.218 | 0.444 |    430 |
| Nb₃F₈  LT-bulk  |    0.05 |  1876.0 | 36643 | −0.013 | 2.5e-3 | 1.224 | 0.444 |  16276 |
| Nb₃F₈  LT-bil   |    0.05 |  2001.1 | 40046 | −0.012 | 2.5e-3 | 1.224 | 0.444 |  17788 |

(θ_CW = −J/4 exactly, from the high-T expansion; E_s = the first ionic singlet above the triplet;
T₅% = the temperature where the exact χ departs from Bleaney–Bowers by 5 %. χ(300 K) here in
emu/mol.)

Observations (cluster-level):

- **The pure-spin (Bleaney–Bowers) description breaks at the charge scale, not the exchange
  scale.** The 5 %-deviation temperature is `T₅% ≈ 0.40·E_s` across the *entire* family
  (0.407–0.444), i.e. set by the first ionic singlet E_s, essentially independent of J. So the
  Heisenberg dimer is a good model for χ(T) whenever `T ≲ 0.4·E_s` — and since every E_s ≥ 774 meV
  (≫ 300 K = 26 meV), **all members are in the good-Bleaney–Bowers regime at and well above room
  temperature.** This supports using a spin-only dimer model to interpret Nb₃X₈ susceptibility.
- **In *reduced* (T/J) terms the iodides are the marginal case:** E_s/J ≈ 3 means charge and spin
  scales are least separated, so Nb₃I₈ departs from Bleaney–Bowers at only T₅/J ≈ 1.3, versus
  ~3–8 for the bromides/chlorides and ~10⁴ for the fluorides — the same "iodides are the most
  correlated / least Heisenberg-like" ordering seen in the exchange-J and exciton-binding trends
  (Section 3).
- **Room-temperature moments span the correlation crossover:** the weakly-coupled fluorides and
  HT phases sit near the free-coupled-pair Curie value (μ_eff ≈ 1.22, χ ≈ 2.5×10⁻³ emu/mol),
  while the strongly-dimerized iodides are essentially diamagnetic-looking singlets at 300 K
  (μ_eff ≈ 0.02) — directly comparable to SQUID data, and consistent with the Nb₃Cl₈ singlet
  behaviour reported by Sheckelton and Haraguchi.

## 5. Conclusion

- **Reference numbers:** exact charge gaps for the Nb₃X₈ bilayer clusters (Section 1) — possibly a
  useful cross-check for the downfolded models.
- **Method assessment:** for the Nb₃X₈ solid gaps, Hubbard-I (in cluster-DMFT, with the full
  dispersion) appears **robust** — the exact solid-gap estimate lands near it. This is consistent with
  the paper's own use of Hubbard-I at integer filling and its switch to CTHYB for the harder/doped
  cases.
- Net: a null / confirming result. The interesting turn was methodological — the isolated-cluster gap
  looked like a 29 % Hubbard-I failure until band broadening was put back, at which point it
  dissolved.
- **New reference constants (Section 3):** exact optical gaps, exciton bindings, and inter-layer
  exchange couplings for all ten parameter sets, two of them in closed form — cheap to cross-check
  against any future downfolding revision, and possibly useful anchors for optics/magnetism work
  on these materials.

## Caveats

- Two-orbital density-density model (the paper reports non-density-density terms of a few meV); no
  phonons or long-range Coulomb tail beyond the strong inter-layer term. In particular the
  fluorides' `J ≈ 0.05 meV` sits *below* those neglected terms and should be read as "≈ 0".
- Section 3's optical/exciton/exchange constants are **isolated-cluster** numbers like Section 1's
  gaps: in-plane (kagome) exchange and inter-dimer screening/broadening are absent, and `P = n₁−n₂`
  stands in for the dipole of a geometry-free model (trends physical, absolute intensities
  model-defined). They are reference values for the downfolded Hamiltonian, not solid-state
  predictions.
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
`tests/test_nb3x8_gaps_spec.py`; the numbers regenerate with `python nb3x8_gaps.py`. Section 3:
closed forms + test-gated derivations in `odmd_optical.py` / `odmd_spin.py`
(`tests/test_odmd_optical_spec.py`, `tests/test_odmd_spin_spec.py`, including the selection-rule
and atomic-limit gates); the tables regenerate with `python odmd_optical.py` and
`python odmd_spin.py`.
