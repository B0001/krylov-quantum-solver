# SPEC: H₁₀ equation of state — reproduce the published 10-point R grid (FCI check)

**Status:** DONE (2026-09-27, bead chem-4y9). G1 passes at all 10 grid points; max |Δ| = 5.25e-6 Ha,
20x inside the 1e-4 Ha gate.

---

## 1. Goal

`SPEC_hchain_tdl.md` and `SPEC_hchain_largen2.md` have only ever validated the H_n geometry
convention (uniform spacing, Bohr→Angstrom conversion, STO-6G placement) at one bond length,
R = 1.8 Bohr — the cohesive minimum. This spec closes `specs/BACKLOG.md`'s "one fixed R was the
easy point" part (a): run FCI at n=10 across the full **Simons Collaboration 10-point R grid**
(Motta et al., *PRX* 7, 031059, 2017 — R = 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6 Bohr,
spanning metallic to Mott-insulating) and check agreement with the paper's own published FCI
column to < 1e-4 Ha at every point. The claim is false if any point misses.

Part (b) of the same backlog item ("is D set by R, not n") was answered separately by bead
chem-1uu ([`SPEC_hchain_largen2.md`](SPEC_hchain_largen2.md) §12); this spec is part (a) only.

## 2. Background and honest framing

- **Reference.** Motta et al., *PRX* 7, 031059 (2017), arXiv:1705.01608. Table II ("Potential
  energy curve of H₁₀ with the minimal (STO-6G) basis") lists FCI as its own column at all 10 R
  values, energy **per atom**, to 6 decimals. The paper states this table is exact: "DMET[5], MRCI
  and MRCI+Q energies coincide with FCI to within 10⁻⁶" — i.e. four independent method families
  agree with FCI at this system size, so the table is a solid ground truth, not a single method's
  output. §IV.A additionally states the exact equilibrium values R_e = 1.786 Bohr,
  E₀ = −0.542457 Ha (per atom) — consistent with Table II's R=1.8 row (−0.542439, the same
  curve one step off the minimum), which is the cross-check that first confirmed the extracted
  table was read correctly.
- **What we can claim if G1 passes:** our `integrals()` geometry/basis convention
  (`benchmark_hchain_tdl.py`: uniform spacing along z, Bohr→Angstrom via CODATA 2018,
  STO-6G, canonical RHF→FCI) reproduces the published curve's *entire* R range, not just the one
  point every other H-chain spec in this repo assumes — including deep into the Mott crossover
  (R=3.6) where entanglement and multireference character are largest.
- **What we cannot claim:** anything about DMRG accuracy or the n→∞ limit — this spec is FCI-only,
  n=10, and exists purely to pin the geometry/basis convention. It is also not novel: this is a
  reproduction of a 2017 published benchmark table, to the precision that table was printed at
  (6 decimals ⇒ ~5e-7 Ha/atom rounding noise, ~5e-6 Ha total for n=10).

## 3. Approach

`benchmark_hchain_tdl.integrals(n=10, localize=False, R_bohr=R)` builds canonical-RHF full-space
STO-6G integrals for H₁₀ at each R in the grid; `hybrid_quantum_solver.dmrg_reference.fci_energy`
solves for the exact total energy (block2/DMRG plays no role — this is a pure-PySCF FCI check, so
the gate test needs no process isolation). The reference is
`specs/hchain_R_grid_motta_fci.csv`, hand-transcribed from Table II of the paper (per-atom FCI
column), with the total (× n=10) computed alongside it for direct comparison to `fci_energy`'s
return value.

## 4. Public interface

No new production code — this spec is a validation gate over the existing public interface:

```
benchmark_hchain_tdl.integrals(n, localize=False, R_bohr=R) -> (h1, eri, ne, e_core, e_hf)
hybrid_quantum_solver.dmrg_reference.fci_energy(h1, eri, ne, e_core) -> float
```

plus the vendored table `specs/hchain_R_grid_motta_fci.csv` (columns: `r_bohr`,
`fci_per_atom_ha`, `fci_total_ha_n10`).

## 5. Acceptance criteria (validation gates)

- **G1 — grid agreement (definition of done).** For every R in
  {1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6} Bohr, at n=10, canonical orbitals:
  `|fci_energy(*integrals(10, R_bohr=R)[:4]) - fci_total_ha_n10(R)| < 1e-4 Ha`.
- **G2 — canonical/localized invariance holds across the grid, not just at R=1.8.** For the same
  10 R values, `|fci_energy(canonical) - fci_energy(localized)| < 1e-8 Ha` (a unitary rotation of
  the full space; this is `test_hchain_tdl_spec.py::test_G6` generalized from 3 R values to all 10
  — cheap, and it would have caught a geometry bug that only one orbital basis happened to hide).

## 6. Implementation plan (test-first)

1. Vendor `specs/hchain_R_grid_motta_fci.csv` from the paper's Table II (already done — see §2).
2. Write `tests/test_hchain_R_grid_spec.py` encoding G1-G2 (pure PySCF, no block2 — no process
   isolation needed, unlike the DMRG gates in `test_hchain_tdl_spec.py`).
3. No new implementation code: `integrals()` and `fci_energy()` already exist and already take
   `R_bohr` (added by chem-1uu). Run to green.

## 7. Out of scope

- **DMRG at these R values, or the n→∞ limit at R != 1.8** — `SPEC_hchain_largen2.md` §12 already
  covers D-vs-R mechanism at n=20 for R ∈ {1.8, 2.4, 3.6}; extending the *headline* TDL curve to
  every R in this grid is a much larger, separate compute effort (not this spec).
  Cf. `specs/BACKLOG.md`'s note on chem-1uu: "the R=3.6 localized point may be closer to a trivial
  near-atomic limit than a hard Mott benchmark point."
- **The published table's non-FCI columns** (AFQMC, CCSD variants, embedding methods, etc.) —
  irrelevant to this repo's own solver, not transcribed.
- **N != 10** — the paper's Table II is only for the N=10 chain; larger N in the paper appear only
  as figures/CBS-extrapolated fits, not machine-checkable tables.

## 8. Caveats and risks

- **R1 (transcription risk):** the reference table was extracted from the arXiv PDF (1705.01608)
  by hand/script, not copy-pasted from a machine-readable source (the paper's own supplementary
  data repository, ref. [11], was listed in the preprint as "link to be finalized" and was not
  independently located). Mitigated by: (a) the internal cross-check against the paper's stated
  exact equilibrium value in prose (R_e=1.786, E0=-0.542457 vs the table's R=1.8 row -0.542439 —
  consistent, since 1.8 is one step off the true minimum at 1.786); (b) our own R=1.8 total energy
  (-5.424385...) already independently matches `data/hchain_tdl.csv`'s pre-existing FCI value at
  n=10, computed for `SPEC_hchain_tdl.md` long before this spec existed.
- **Minimal basis** (STO-6G): this validates the geometry/basis *convention*, not the physical
  H-chain — same caveat as every other H-chain spec in this repo.

## 9. Deliverables

- `specs/hchain_R_grid_motta_fci.csv` — vendored reference table (Motta et al. Table II, FCI
  column).
- `tests/test_hchain_R_grid_spec.py` — gates G1-G2.
- This spec, recording the result (§ Status) and the transcription caveat (R1).
