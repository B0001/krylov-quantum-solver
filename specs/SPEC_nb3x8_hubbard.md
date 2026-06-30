# SPEC: Nb3X8 / Hubbard model loader into the universal integral interface

**Status:** OPEN — gates G1–G5 in `tests/test_nb3x8_hubbard_spec.py`, implemented in
`hybrid_quantum_solver/model_hamiltonians.py`. Derived from `CLAUDE_CODE_HANDOFF.md` (the Nb3X8
Model-Database validation shortcut).

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

Add a loader that maps a tight-binding **hopping matrix** + a **Hubbard/cRPA interaction** into the
project's universal active-space tuple `(h1, eri, e_core, nelec, norb)`, so the *same* validated
solver stack runs on model Hamiltonians with no heavy-atom PySCF/CASCI. Claim: for the half-filled
two-site Hubbard dimer the number-conserving solver reproduces the **exact analytic** ground-state
energy `E0 = (U − √(U² + 16t²))/2` across the full weak-to-strong-correlation range, and the
general (rank-4 cRPA) mapping reproduces PySCF FCI in the correct particle-number sector.

## 2. Background and honest framing

`CLAUDE_CODE_HANDOFF.md` (literature + code-review session) identifies the **Nb3X8 Model Database**
(Aretz, Grytsiuk, Strand, van Loon, Rösner 2025; ref [86] of `arXiv:2501.10320`, PRX DOI
10.1103/wr7w-nfhg) as the cleanest validation target: it publishes hopping matrices and the full
rank-4 cRPA Coulomb tensors `U_ijkl` per compound × structure, which drop straight onto our
`(h1, eri)` interface — sidestepping open-shell Nb SCF entirely. The bulk low-T electronic structure
reduces to a generalized **Hubbard dimer** (two dimerized trimer MOs, singlet ground state), whose
ground-state energy is analytic — a reference that can embarrass the loader.

- **What we can claim if the gates pass:** the loader correctly maps lattice/cRPA parameters into the
  vetted Jordan–Wigner interface — checked against an analytic Hubbard reference *and* PySCF FCI — and
  the existing number-conserving Krylov/SKQD solvers run on it unchanged. It also pins a real pitfall
  (full-Fock-space diagonalization gives the wrong filling).
- **What we cannot claim:** this is **not** a validated-against-DMFT materials result. The database's
  own DMFT/Hubbard-I gaps are not bundled (the DB file format is not in the repo — contact
  m.roesner@science.ru.nl). A named compound's reported singlet-triplet gap is a *Hubbard-dimer model
  prediction*, not a benchmarked number. No novelty, no quantum advantage at this scale (2–N sites,
  FCI-tractable). Per the handoff, Nb3I8 is the **weakly** correlated family member (U/t ≈ 4), kept
  here as the easy analytic anchor, not a "strong correlation" headline.

## 3. Approach

`hubbard_integrals(hopping, interaction, nelec, units=…)` → `ModelIntegrals(h1, eri, e_core, nelec,
norb)`. The on-site Hubbard `U` maps to the chemist-notation diagonal `eri[i,i,i,i] = U` (verified:
this reproduces `U n_{i↑}n_{i↓}` exactly); a full rank-4 tensor is used as `eri` verbatim. Units
(`Ha`/`eV`/`meV`) are scaled to Hartree. `load_from_nb3x8_database` wraps it: database mode (parsed
`hopping`+`coulomb`) or named-compound mode (published bulk-LT dimer params). References: the analytic
dimer formula and `dmrg_reference.fci_energy` (in-sector FCI). Reuse `build_hamiltonian_from_integrals`
and `QuantumKrylovSolver` verbatim — no new physics core.

## 4. Public interface

```
hybrid_quantum_solver.model_hamiltonians
  hubbard_integrals(hopping, interaction, nelec, *, e_core=0.0, units="Ha") -> ModelIntegrals
  ModelIntegrals(h1, eri, e_core, nelec, norb).as_tuple() / .to_hamiltonian()
  fixed_filling_energy(model) -> float                       # in-sector PySCF FCI reference
  hubbard_dimer_energy(t, U) -> float                        # analytic (U − √(U²+16t²))/2
  hubbard_dimer_gap(t, U) -> float                           # singlet-triplet gap (triplet at 0)
  load_from_nb3x8_database(compound="Nb3I8", *, hopping=None, coulomb=None, nelec=None,
      e_core=0.0, units="meV") -> ModelIntegrals
  NB3X8_BULK_DIMER_PARAMS                                     # published bulk-LT dimer params (meV)
```

## 5. Acceptance criteria (validation gates)

Gates live in `tests/test_nb3x8_hubbard_spec.py` (small/fast; pyscf+qiskit, no block2).

- **G1 — analytic dimer, all U/t.** For the half-filled 2-site Hubbard dimer at `t=1` and
  `U ∈ {0, 2, 4, 8, 20}`, the number-conserving `QuantumKrylovSolver` energy matches
  `(U−√(U²+16t²))/2` to `< 1e-6 Ha`. Covers the U=0 (`−2t`) and deep-Mott (`−4t²/U`) limits.
- **G2 — variational floor.** At every U, `E_solver ≥ E_analytic − 1e-9` (a number-conserving
  Rayleigh–Ritz estimate cannot beat the in-sector ground state).
- **G3 — rank-4 / cRPA mapping vs FCI.** For a 3-site model with a *full* rank-4 interaction tensor
  (off-diagonal density-density `U_ij` included) and a generic hopping matrix, the loader's
  `(h1, eri, nelec)` reproduces PySCF FCI to `< 1e-9 Ha` — the mapping is faithful, not just the
  diagonal special case.
- **G4 — fixed-filling pitfall (a finding, asserted).** The full-Fock-space
  `MolecularHamiltonian.ground_state_energy()` does **not** equal the half-filled energy once `U/t`
  is large (it returns the 1-electron bonding state `−t`), whereas `fixed_filling_energy` and the
  number-conserving solver **do**. Encodes SKQD checklist item 5 (U(1) conservation) as a guardrail.
- **G5 — units + Nb3I8 anchor.** `load_from_nb3x8_database("Nb3I8")` (meV inputs) round-trips through
  the eV→Ha scaling: its solver energy equals `hubbard_dimer_energy(t, U)` (computed in meV, scaled
  to Ha) to `< 1e-9 Ha`, and the singlet-triplet gap is `(√(U²+16t²)−U)/2 ≈ 194 meV` (> 0, singlet
  below triplet — qualitatively matches the published "singlet ground, triplet first excited").

## 6. Implementation plan (test-first)

1. `tests/test_nb3x8_hubbard_spec.py` encoding G1–G5 (initially RED).
2. `hybrid_quantum_solver/model_hamiltonians.py` composing `build_hamiltonian_from_integrals`,
   `QuantumKrylovSolver`, and `fci_energy`.
3. Green via `make gates` (own process; shares pyscf, stays out of the block2 group).

## 7. Out of scope

- Parsing the actual Nb3X8 Model-Database file format (not bundled; loader takes already-parsed
  arrays). When the files are obtained, add a thin format parser + a per-compound DMFT-gap
  regression — its own spec.
- Periodic/lattice (k-space) Hubbard at thermodynamic-limit sizes; this is finite-cluster only.
- Any DMFT/Hubbard-I comparison (needs the database; see §2).
- Retargeting to Nb3Cl8 (the genuine single-orbital Mott target) — a follow-up once the DB anchor is
  in (handoff §Physics).

## 8. Caveats and risks

- **R1 — full-Fock-space trap.** Diagonalizing the qubit Hamiltonian over all particle numbers gives
  the wrong state for a Hubbard model; G4 pins this. Always use the number-conserving path.
- **R2 — eri notation.** `eri[i,i,i,i]=U` is chemist notation `(ii|ii)`; verified to reproduce
  `U n↑n↓` (G1/G3 would fail otherwise). A physicist-notation tensor must be reordered before load.
- **R3 — sign of t.** Only `t²` enters the analytic energy, so the hopping sign is irrelevant to the
  gap; it matters for which orbital is bonding (HF reference). The gates use `−t` off-diagonals.

## 9. Deliverables

- `hybrid_quantum_solver/model_hamiltonians.py` — the loader + analytic references.
- `tests/test_nb3x8_hubbard_spec.py` — gates G1–G5.
- Backlog entry + `__init__` exports; results summary (with §2 caveats) in the PR description.
