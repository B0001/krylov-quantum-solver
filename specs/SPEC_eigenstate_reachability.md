# SPEC: a per-eigenstate symmetry decision for HF reachability — killed as a reference veto, kept as a diagnostic

**Status:** REVISED after merge (chem-obf fix-forward, 2026-10-03). PR #62 used this decision to
*veto* the levels it judges symmetry-forbidden from every reference. A post-merge review showed that
QKSD from |HF⟩ converges to exactly the level it vetoed, so **the veto premise is KILLED for
references** (§0, gate G9): a reference must describe what the data targets. `reachable_mask` is
the population cut again; the decision survives only as a best-effort **diagnostic** that warns,
never changes a mask and never raises. G1–G8 are re-pinned to those semantics and G9–G10 added —
pre-registered in this revision, committed before the revised gate file ran. Results: §9.

---

## 0. The kill (post-merge review, 2026-10-03)

What PR #62 claimed holds as *physics*: the residue level is symmetry-forbidden to the exact Ag HF
determinant, and the decision identifies it against an independent symmetry-adapted FCI (G1, G2).
What it got wrong is what a reference is for. The qubit Hamiltonian is built from the *loose* SCF's
MOs. In that basis |HF⟩ really carries the residue amplitude, exact evolution conserves it, and the
Krylov space amplifies it like any other component until the solver lands on the lowest level
|HF⟩ touches. Scout (macOS, `QuantumKrylovSolver(mh)` defaults, electronic energies; script not
committed — G9 re-derives what it asserts):

| witness | residue pop | first M with E(M) − E_B1g < 1e-2 Ha | E − E_B1g at M = 24 | max over M ∈ [24, 40] |
|---|---|---|---|---|
| sq-H₄ 1.10, `conv_tol` 1e-6 | 2.66e-6 | 13 | 1.3e-4 | 3.8e-4 |
| sq-H₄ 1.19, `conv_tol` 1e-6 | 5.48e-8 | 21 | 1.5e-4 | 1.1e-3 |
| sq-H₄ 1.19, builder default 1e-9 | 1.40e-8 | 22 | 1.5e-4 | 2.0e-3 |
| sq-H₄ 1.10 or 1.19, `conv_tol` 1e-13 | ~4e-30 | never (M ≤ 40) | 0.15 | — (stays on Ag) |

E_B1g = −4.556211 / −4.374274 Ha and E_Ag = −4.404505 / −4.224669 Ha at a = 1.10 / 1.19 (PySCF
symmetry-adapted FCI). PR #62's reference sat ~150 mHa *above* the energy the solver returns: any
"variational floor" built on it is violated by the data it is meant to check.

The review also confirmed two crashes, fixed here: `orbital_parity_bounds` raised on geometries
PySCF half-detects (sq-H₄ a = 1.1 with one atom moved +1.5e-5 Å along x → `IndexError` inside
`symm_adapted_basis`, +2e-5 Å → `PointGroupSymmetryError`), which took down every rewired site; and
a tuple `active_electrons` (odd-electron active spaces, e.g. OH doublet `(3, 2)`) raised
`TypeError` in `_driver_orbitals`. Both reproduced by probe on macOS. The other review findings only
matter for a veto; they are the diagnostic's known limitations (§8).

**The fix (forward, not a revert).**
`reachable_mask(mh, w, vecs, pops, tol) = pops > tol` — the cut the sites made before PR #62. The
symmetry decision runs inside it as a diagnostic: a population-admitted level it judges forbidden
triggers a `RuntimeWarning` naming the level's energy and population, calling the population SCF
residue and the remedy a rebuild at `conv_tol=TIGHT_SCF_CONV_TOL`. Any exception inside the
diagnostic means "no diagnostic". A tuple `active_electrons` is summed, as `ActiveSpaceTransformer`
does. Kept from PR #62: every site's routing through `_dense_hf_projection` (the macOS ZHEEVD
fallback that fixed five gate files) and every signature. Where the residue really belongs — a
tight or symmetric SCF in the builder — is a follow-up.

## 1. Goal (original, PR #62)

`SPEC_reachability_tolerance` §2b: no fixed population threshold separates a physical HF overlap
from SCF-convergence residue on a symmetry-forbidden level (square H₄: the residue sits at 5e-10 at
a = 1.1 Å and at 1.4e-8 at a = 1.19 Å, while LiH carries *physical* levels between 1e-10 and 1e-8).
**Claim:** an exact symmetry decision per eigenvector — *is the majority of ψ_k's weight in the HF
determinant's symmetry sector?* — removes every residue level the population cut admits on those
witnesses, agrees with an independent symmetry-adapted FCI, and changes **nothing** on the ordinary
gated systems (the bead's prediction). *Since the revision:* the identification and "nothing
changes on ordinary systems" halves survive (G1–G3); "removes" is what §0 kills.

## 2. Background and honest framing

The premise "SPEC_symmetry_reachability built a symmetry-aware decision" is only half true: it
built `symmetry_filter_available(atom)` (one bool per *system*). Nothing labelled eigenstates, and
`MolecularHamiltonian` carried no geometry, so no site could have consumed it. This spec builds the
missing piece first.

**Why the residue exists at all.** With a loose SCF the driver's MOs are not exactly symmetry
adapted, so the qubit Hamiltonian's eigenvectors (exact symmetry eigenstates of the *physical* H)
overlap the slightly asymmetric HF bitstring by a residue amplitude. The residue is
platform-dependent — measured here on macOS: square H₄ a = 1.1 Å carries **2.7e-6** at
`conv_tol=1e-6` and **none** at the default 1e-9, where Linux records 5.07e-10 (this is why
`test_reachability_tolerance_spec` G1–G3 fail on macOS).

**Design-time probe numbers** (macOS, before the gate file was written; scripts not committed — the
gate file re-derives every number it asserts):

| system | MO impurity of the residue bit | β (bound, §3) | population cut admits | decision |
|---|---|---|---|---|
| sq-H₄ 1.10, ct 1e-6 | 1.4e-6 | 9.5e-3 | B1g @ −4.556211, p = 2.7e-6 | flagged; allowed set = Ag FCI (10 states, 6e-15) |
| sq-H₄ 1.19, ct 1e-9 | 7.6e-9 | 7.0e-4 | B1g @ −4.374274, p = 1.4e-8 | flagged; allowed set = Ag FCI |
| sq-H₄ 1.20 (broken RHF) | 0.5 | 6.3 | — | bit excluded; naive labels would flag p = 3.7e-2 |
| LiH (π pair rotated by the free SCF) | 0.05–0.16 | 1.7–3.3 | — | σᵥ bit excluded, C₂ bit used; nothing flagged |

## 3. Approach

1. **Geometry.** `MolecularHamiltonian.build_args` (new, default `None`) records the builder's
   arguments. The decision re-runs that exact PySCF driver (+ active-space transform + JW map) and
   uses its MO coefficients **only if the rebuilt operator reproduces `mh.qubit_hamiltonian`**
   (non-identity Pauli coefficients equal to 1e-12 — a reproducibility check of a deterministic
   computation, not a physics threshold; the identity term is ignored so a centered frame still
   verifies). Otherwise no spatial symmetry is used.
2. **Parity bits.** For each bit b of PySCF's abelian irrep ids (D2h subgroup; linear groups via
   `id % 10`), the symmetry operation g_b in the MO basis is `R_b = Σ_Γ χ_b(Γ) Π_Γ`
   (S-metric irrep projectors from `symm.symm_adapted_basis` in the driver's own frame). Its
   Z-string surrogate uses the labels `d = sign(diag R_b)`. Bound: `||Γ(D_b) − Γ(R_b)||` on
   N-electron states ≤ β_b = sum of the N largest spin-orbital `|1 − λ|` over eigenvalues λ of
   `D_b R_b` (rigorous for the full space, where `D_b R_b` is orthogonal; a heuristic for active
   spaces — and a wrong one when the window splits a degenerate shell, §8 L1).
3. **Sector.** HF's exact sector σ = same (N_α, N_β) as the HF bitstring, plus the parity of every
   bit in the greedy set with **Σβ < 1** (smallest first). Σβ < 1 is not a tuned constant: the
   product projector then differs from the exact one by < 1/2 in norm, so the majority decision
   below cannot flip for any non-degenerate symmetry eigenstate.
4. **Decision — now a diagnostic.** Per cluster C of numerically degenerate eigenvalues (gap ≤ 1e-9
   Ha; `eigh` splits exact degeneracies by ~1e-14), on HF's basis-independent projection
   φ_C = V_C V_C†|HF⟩: `allowed_C ⟺ ||P_σ φ_C||² > ||φ_C||²/2` (the rounding of a 0/1 quantity; for
   a non-degenerate level it is the majority test on ψ_k's sector weight). It is defined only where
   HF's projection is above roundoff. *Revised after the first code review:* the first version
   tested each eigenvector's own sector weight, and `eigh` mixes exactly degenerate open-shell
   M_s = ±S pairs (different (N_α, N_β) sectors) — it dropped populated levels (H₄ triplet: 2 with
   p > 1e-3; OH doublet 38 → 22 kept at 1e-8). Gated as G8. *Revised after merge (§0):* `reachable_k ⟺ pop_k > tol`, and the
   diagnostic warns on `pop_k > tol ∧ ¬allowed_C(k)`.

**Reference:** PySCF symmetry-adapted FCI (`fci` with `wfnsym`, tight symmetric SCF) — independent
of the qubit path — for the diagnostic's verdicts, and `QuantumKrylovSolver` from |HF⟩ for what the
data targets (G9).

## 4. Public interface

```
MolecularHamiltonian.build_args : dict | None          # set by build_molecular_hamiltonian
reachability.orbital_parity_bounds(mh) -> list[(odd: bool[n_orb], beta: float)] | None   # diagnostic
reachability.hf_symmetry_sector(mh)    -> bool[2**n] | None   # diagnostic; None: not a JW HF layout
reachability.symmetry_allowed(mh, w, vecs) -> bool[n_eig]      # diagnostic, per degenerate cluster
reachability.reachable_mask(mh, w, vecs, pops, tol) -> pops > tol   # + RuntimeWarning on a flagged level
reachability.reachable_eigenpairs(mh, tol)   # the population cut again (the pre-PR-#62 numbers)
```
`_dense_hf_projection`, `hf_population_spectrum` and the signatures above are unchanged.

## 5. Acceptance criteria (validation gates)

`tests/test_eigenstate_reachability_spec.py`. **Revised and pre-registered 2026-10-03**, after the
§0 scout and before the revised gate file ran; thresholds kept from PR #62 wherever its claim
survives. Witnesses: sq-H₄ a = 1.10 and 1.19, **both at `conv_tol=1e-6`**. (a = 1.19 moves off the
builder default, where its 1.40e-8 residue cleared the 1e-8 site cut by only 1.4×; at 1e-6 the scout
measured 5.48e-8, unchanged for `conv_tol` 1e-5…1e-8. *Revised after code review, before merge:*
5.48e-8 is still only 5.5× above 1e-8 and macOS-measured, so the checks at the 1e-8 cut use
a = 1.10 only — G1 gates a = 1.19 at 1e-10 (548×), and G7 compares each site with the lowest level
above its own cut instead of naming B1g at 1e-8.)

- **G1 — the witnesses: identified, warned, kept.** At tol ∈ {1e-8, 1e-10} (a = 1.19: 1e-10 only,
  see above): `reachable_mask == (pops > tol)`; the call emits a `RuntimeWarning` naming SCF
  residue; the lowest populated level `symmetry_allowed` rejects is the B1g FCI ground state, and
  the lowest populated level it allows is the Ag one (both |ΔE| < 1e-8 Ha).
- **G2 — brute force vs an independent reference, every populated level.** sq-H₄ a ∈ {1.05, 1.10,
  1.19, 1.35} × `conv_tol` ∈ {1e-6, 1e-9}, linear H₄ and H₂: the eigenvectors in HF's sector
  (per-vector weight > 1/2) reproduce PySCF's D2h Ag FCI spectrum (same count, max |ΔE| < 1e-8 Ha),
  **and** `symmetry_allowed` is True exactly on the eigenvectors with population > 1e-20 whose
  energy is in that spectrum (|ΔE| < 1e-8). The floor 1e-20 (amplitude 1e-10, six orders above
  double-precision roundoff) is where HF's projection, hence the decision, is defined; the scout
  gave the same verdicts at floors 1e-30…1e-16 and up to 6 mismatches per case at floor 0 (levels whose
  HF overlap is roundoff, where the decision is undefined).
- **G3 — no false alarm on ordinary systems.** H₂ (0.74, 2.0), linear H₄, HeH⁺, LiH (12 qubits),
  LiH CAS(2,5), N₂ CAS(6,6), at tol ∈ {1e-8, 1e-10}: `symmetry_allowed` is True on every populated
  level, and `reachable_mask` returns the cut with no warning.
- **G4 — the bound is load-bearing.** sq-H₄ a ∈ {1.20, 1.40} (broken RHF): some bit has β ≥ 1 and is
  excluded; nothing populated is flagged; at a = 1.20 the naive labels (every bit, no bound) would
  flag a level with population > 1e-3. LiH: the C₂ bit is exact (β < 1e-10). *Dropped:* PR #62's
  "the σᵥ bit has β ≥ 1" — it measures the arbitrary rotation SCF/LAPACK leave inside LiH's
  degenerate π pair, not LiH: β = 3.2 at the default `conv_tol`, 0.646 at 1e-7, 0.606 with the atom
  order swapped (probe; matches the review).
- **G5 — fallback.** No `build_args` (integrals path): only (N_α, N_β) is used, no warning, the cut
  is returned. A `build_args` that does not reproduce the operator (a different geometry swapped
  in) disables the spatial bits; a constant shift does not.
- **G6 — never near the boundary, and the decision is its rounding.** On sq-H₄ 1.10 and 1.19 (at
  1e-6), 1.20, linear H₄, LiH and N₂ CAS(6,6): eigenvectors with population > 1e-10 have sector
  weight within 1e-3 of 0 or 1, **and** `symmetry_allowed` equals (weight > 1/2) on them.
- **G7 — the sites use the population cut.** At the witnesses: the centered frames of `odmd`, `msd`,
  `trotter_odmd`, `device_odmd`, `trotter_resolution_floor` equal the 1e-8 population-cut frame
  (μ abs 1e-10, τ rel 1e-10); `reachable_eigenpairs(mh)[0][0]` is the B1g FCI energy (1e-8); and
  `hf_overlap_certificate.exact_reachable_overlap` and `hf_overlap_subspace.exact_hf_subspace_overlap
  (mh, 1)` equal √pop of the lowest level above their own cuts, 1e-10 and 1e-8 (rel 1e-8).
- **G8 — open shells, basis invariance.** H₄ triplet and OH doublet: `symmetry_allowed` is True on
  every populated level at both tols; on the triplet it stays so in a deliberately M_s-mixed
  degenerate eigenbasis, in which the per-vector rule flags a level with p > 1e-3 (non-vacuous).
- **G9 — the falsification (DEFINITION OF DONE of this revision).** From |HF⟩, with
  `QuantumKrylovSolver(mh)` defaults, at sq-H₄ 1.10 / 1e-6, 1.19 / 1e-6 and 1.19 at the builder's
  default `conv_tol`: for **every M ∈ [28, 32]**, −1e-9 ≤ E(M) − E_B1g < **1e-2 Ha** (electronic;
  the lower bound is the variational floor), i.e. > 0.1 Ha below the Ag level the veto kept; and
  `reachable_eigenpairs(mh)[0][0]` is E_B1g (1e-8) — the reference keeps the level the solver
  reaches. Control, a ∈ {1.10, 1.19} at `TIGHT_SCF_CONV_TOL`: E(M) ≥ E_Ag − 1e-8 for every M ≤ 32,
  and `reachable_eigenpairs` starts at E_Ag (1e-8). The [28, 32] window sits ≥ 6 steps past the
  scout's first convergence (≤ 22) and 1e-2 Ha is ≥ 5× the scout's worst value on it.
- **G10 — the diagnostic never raises.** `reachable_mask` returns the cut, without raising, on the
  §0 crash geometries (sq-H₄ 1.1 with the second atom at x = 1.100015 and 1.10002 Å) and when
  `symmetry_allowed` itself raises (monkeypatched — the non-vacuous check). OH doublet,
  `active_electrons=(3, 2), active_orbitals=4`: `orbital_parity_bounds` verifies the rebuilt
  operator (not None) and nothing populated is flagged (the scout: 18 levels with p > 1e-10, none
  flagged).

## 6. Implementation plan (test-first)

1. This spec, then the gate file (fails: functions missing).
2. `build_args` on `MolecularHamiltonian`; the decision in `reachability.py`.
3. Rewire: `reachable_eigenpairs` (→ certified_gaps, hf_overlap_certificate, certified_dipole,
   certified_noise, excited_bounds), `hf_overlap_subspace`, `odmd`, `msd`, `trotter_odmd`,
   `device_odmd`, `trotter_resolution_floor` — each through `_dense_hf_projection` (which also
   removes their macOS ZHEEVD crash at 12 qubits).
4. Re-run every affected gate; record any moved number below.
5. *Revision:* this spec (§0, §5) committed first; then `reachable_mask` → population cut +
   diagnostic, the tuple fix, docstrings, the re-pinned gate file; then the before/after gate run.

## 7. Out of scope — sites deliberately NOT rewired, and why (part of the finding)

*Since the revision every site is on the population cut, so this section now only says where the
diagnostic does not run.* The symmetry argument holds where a site thresholds HF populations **in the
exact-H eigenbasis** (exact evolution conserves them). It does not transfer to:

- `trotter_odmd.select_ground_eigenphase` (`pop_cut`) and `trotter_resolution_floor._eigenphase_energy`
  (`pops_u > 1e-8`): populations in the **Trotter circuit's** eigenbasis. A product of single-Pauli
  exponentials conserves the Z-string parities but **not N or S_z** (individual JW terms do not), so
  the HF sector is not invariant there — that leakage is real signal the probes measure.
- `odmd_spectral.reference_signal` (:79): a **kicked** reference in a different particle-number
  sector; HF's sector is the wrong sector by construction.
- Arbitrary-reference / model-Hamiltonian sites (`odmd_spectral` :120/:153, `odmd_optical` :58/:92,
  `odmd_spin` :79, `visibility_law` :110, `nb3x8_alloy` :105): built from integrals, no geometry —
  the decision reduces to N/S_z, which their populations already respect. Unchanged.
- `rodeo.py` :71/:87 (literal 1e-8; `overlap_tol` at :49 unused): molecular but only H₂/H₄ (no
  residue); left for a follow-up rather than widening this diff.

## 8. Known limitations of the diagnostic, and other caveats

The diagnostic can warn wrongly or stay silent; since the revision it never changes a number.

- **L1 — compressed R_b (review).** β is rigorous only where `D_b R_b` is orthogonal (the full
  orbital space). An active window that splits a degenerate shell compresses R_b far from orthogonal
  and β under-reports: symmetric CH₄ (Td, C–H 1.087 Å) CAS(2,2) flags a populated level (probe:
  −0.669500 Ha electronic, p = 3.3e-6) that `QuantumKrylovSolver.solve_excited(4)` returns as θ₁ —
  under the PR #62 veto `reachable_gap` was off by 0.81 Ha (review). The review also reports NH₃
  CAS(4,4) flagging populations up to 8.1e-4; not reproduced at our probe geometry (PySCF found only
  Cs there). Upgrade path: an orthogonality check on R_b.
- **L2 — symmetry PySCF detects only within its geometric tolerance is treated as exact (review).**
  sq-H₄ 1.1 with the fourth atom moved 5e-6 / 1e-5 Å along x, `conv_tol=1e-13`: a level whose
  population is real distortion, 2.1e-9 / 5.8e-9 (probe), is flagged as residue. Upgrade path:
  tolerance-aware symmetry detection.
- **L3 — mapper not recorded (review).** `hf_symmetry_sector`'s Jordan–Wigner layout guard checks
  only the HF index, which an `InterleavedQubitMapper` operator can share (linear H₄, spin = 2), and
  `build_args` does not record the mapper. No in-repo caller. Upgrade path: record the mapper.
- **L4 — accidental degeneracy.** Levels of different symmetry closer than 1e-9 Ha form one cluster
  and share its majority verdict: a missed or a spurious warning for one of them. Never seen on a
  gated system.
- **R2 — residual symmetry outside PySCF's abelian frame.** At broken-RHF geometries a level can
  still carry convergence-dependent population (sq-H₄ 1.20: 1.5e-9 → 6.2e-9 between `conv_tol` 1e-9
  and 1e-6) that the D2h bits cannot see; the diagnostic is silent there.
- **R3 — cost.** One extra SCF + JW map per Hamiltonian (cached on the instance), as in PR #62.
- **R4 — Linux.** With the population cut restored, the Linux-only failures PR #62 predicted for
  `test_reachability_tolerance_spec` G2 and G5 (`constant_governs_the_certified_reference`) should
  not happen: those gates pin the residue the cut keeps. Not verifiable on this macOS host.

## 9. Results (macOS 27, Apple M3, 2 BLAS/OMP threads)

**Revision (2026-10-03):**

- **This gate file:** 45 passed in 35 s (G1–G10; 46 before the code-review change to G1/G7 in §5).
- **G9, the kill:** max over M = 28…32 of E(M) − E_B1g = 3.0e-6 / 1.1e-3 / 1.1e-3 Ha (1.10 at
  1e-6, 1.19 at 1e-6, 1.19 at the default) — the solver sits on B1g, ~0.15 Ha below Ag. Control at
  `TIGHT_SCF_CONV_TOL`: min over M ≤ 32 of E(M) − E_Ag = −4.6e-14 / −8.0e-15 (stays on Ag).
- **Mutation check** (temporary pytest plugin, not committed): `symmetry_allowed` → all False fails
  32/45, including every G2 (10/10) and G6 (6/6) case — PR #62's G2 and G6 passed it (16/16, as
  the review said: they never called the decision). → all True fails 11/45 (G1 3/3, G2 on the 6
  residue cases, G6 on both witnesses).
- **24 affected gate files** (one process each, `GATE_NO_CACHE=1`), origin/main fabaa7d → this
  revision: `test_scf_conv_tol_spec` 12/14 → **14/14** (its G2 [H4 square 1.35] and [1.19] pin
  the default-`conv_tol` residue as the reachable reference; PR #62 broke them). Apart from this
  gate file (35 → 46 tests at that run), every other file unchanged, including the five macOS ZHEEVD fixes (`adaptive_shots` 5/5, `centered_pds` 6/6,
  `msd_sampling` 4/4, `odmd` 4/4, `odmd_uq` 4/4). Still failing with the same test IDs before and
  after, all raw-`eigh` ZHEEVD/`LinAlgError` on macOS: `test_chained_overlap_spec` (9),
  `test_odmd_excited_spec` (1), `test_subspace_floor_resolvability_spec` (13 errors).
  `test_reachability_tolerance_spec`: 9 passed, 3 skipped (macOS) both times.
- **Recorded numbers that moved:** at the residue witnesses, the references are the population
  cut's again (the pre-PR-#62 values); nothing moved on the ordinary systems (G3).

**PR #62 (2026-10-02; superseded):**

- **This gate file:** 35 passed in 46 s (G1–G8).
- **Recorded numbers that moved:** none on the ordinary systems — G3 shows identical masks, and the
  rewired sites call the same `eigh` on the same matrix. The sq-H₄ witness numbers move by design
  (G1, G7). 43 affected gate files re-run before/after: **no new failure**; five files that crashed
  on macOS's ZHEEVD now pass because their dense `eigh` goes through `_dense_hf_projection`:
  `adaptive_shots` 1/5 → 5/5, `centered_pds` 1/6 → 6/6 (chem-aj1), `msd_sampling` 1/4 → 4/4,
  `odmd` 0/4 → 4/4, `odmd_uq` 0/4 → 4/4. Still failing, unchanged and outside this change: raw
  `np.linalg.eigh` inside `test_chained_overlap_spec` (10), `test_odmd_excited_spec` (1),
  `test_subspace_floor_resolvability_spec` (13 errors), and `test_reachability_tolerance_spec`
  G1–G3 (no residue on macOS at the default `conv_tol`).
- **Unit-1 interplay:** `map_canonical` (chem-8tl) and the plain JW map differ by ≤ 1.8e-15 per
  coefficient (LiH), far inside the 1e-12 rebuild check.
- **Cross-PR effect missed (review):** `test_scf_conv_tol_spec` G2 [H4 square 1.35] and [1.19]
  pin the default-`conv_tol` residue as the reachable reference (overlap < 1e-3); with the veto it
  became 0.613 / 0.652 and both failed on main (reproduced on this host's origin/main run).

## 10. Deliverables

- `reachability.py` — the population cut, the diagnostic, the tuple fix; `molecular_hamiltonian.py`
  — `build_args` only.
- Rewired sites listed in §6.3 (on `_dense_hf_projection` + `reachable_mask`);
  `tests/test_eigenstate_reachability_spec.py` — G1–G10.
