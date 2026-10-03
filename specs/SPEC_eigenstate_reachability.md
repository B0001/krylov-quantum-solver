# SPEC: a per-eigenstate symmetry decision for HF reachability, wired into the threshold sites

**Status:** IMPLEMENTED once gates green (chem-obf). Builds the per-eigenstate decision that
[`SPEC_symmetry_reachability`](SPEC_symmetry_reachability.md) did not deliver (it decided only
*whether* a filter is possible), then rewires the `|⟨HF|ψ_k⟩|² > tol` sites onto it.

---

## 1. Goal

`SPEC_reachability_tolerance` §2b: no fixed population threshold separates a physical HF overlap
from SCF-convergence residue on a symmetry-forbidden level (square H₄: the residue sits at 5e-10 at
a = 1.1 Å and at 1.4e-8 at a = 1.19 Å, while LiH carries *physical* levels between 1e-10 and 1e-8).
**Claim:** an exact symmetry decision per eigenvector — *is the majority of ψ_k's weight in the HF
determinant's symmetry sector?* — removes every residue level the population cut admits on those
witnesses, agrees with an independent symmetry-adapted FCI, and changes **nothing** on the ordinary
gated systems (the bead's prediction).

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
| sq-H₄ 1.10, ct 1e-6 | 1.4e-6 | 9.5e-3 | B1g @ −4.556211, p = 2.7e-6 | rejected; allowed set = Ag FCI (10 states, 6e-15) |
| sq-H₄ 1.19, ct 1e-9 | 7.6e-9 | 7.0e-4 | B1g @ −4.374274, p = 1.4e-8 | rejected; allowed set = Ag FCI |
| sq-H₄ 1.20 (broken RHF) | 0.5 | 6.3 | — | bit excluded; naive labels would veto p = 3.7e-2 |
| LiH (π pair rotated by the free SCF) | 0.05–0.16 | 1.7–3.3 | — | σᵥ bit excluded, C₂ bit used; no veto |

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
   spaces, whose compression is only measured orthogonal to ~1e-14).
3. **Sector.** HF's exact sector σ = same (N_α, N_β) as the HF bitstring, plus the parity of every
   bit in the greedy set with **Σβ < 1** (smallest first). Σβ < 1 is not a tuned constant: the
   product projector then differs from the exact one by < 1/2 in norm, so the majority decision
   below cannot flip for any non-degenerate symmetry eigenstate.
4. **Decision.** `allowed_k ⟺ Σ_{x∈σ} |V[x,k]|² > 1/2` (the rounding of a 0/1 quantity), and
   `reachable_k ⟺ allowed_k ∧ pop_k > tol`. The site's own `tol` is kept for allowed levels: below
   it they are physically populated but faint (LiH has 4 such levels in (1e-10, 1e-8]); that
   remaining 1e-8/1e-10 difference is the sites' visibility choice, not the artifact.

**Reference:** PySCF symmetry-adapted FCI (`fci` with `wfnsym`, tight symmetric SCF) — independent
of the qubit path — and the dense population spectrum for the "no change" controls.

## 4. Public interface

```
MolecularHamiltonian.build_args : dict | None          # set by build_molecular_hamiltonian
reachability.orbital_parity_bounds(mh) -> list[(odd: bool[n_orb], beta: float)] | None
reachability.hf_symmetry_sector(mh)    -> bool[2**n] | None   # None: not a JW HF layout
reachability.symmetry_allowed(mh, vecs) -> bool[n_eig]         # majority weight in the sector
reachability.reachable_mask(mh, vecs, pops, tol) -> bool[n_eig]
reachability.reachable_eigenpairs(mh, tol)   # DELIBERATE CHANGE: now applies reachable_mask
```
`_dense_hf_projection`, `hf_population_spectrum` and the signatures above are unchanged.

## 5. Acceptance criteria (validation gates)

`tests/test_eigenstate_reachability_spec.py`. Thresholds pre-registered before the gate ran.

- **G1 — the witnesses (DEFINITION OF DONE).** sq-H₄ a = 1.10 / `conv_tol=1e-6` at tol 1e-10, and
  a = 1.19 / default `conv_tol` at tol 1e-8: the population cut admits ≥ 1 level the decision
  rejects; the lowest population-reachable energy is the B1g FCI ground state and the lowest
  decision-reachable one is the Ag FCI ground state (both |ΔE| < 1e-8 Ha).
- **G2 — brute force vs an independent reference.** sq-H₄ a ∈ {1.05, 1.10, 1.19, 1.35} ×
  `conv_tol` ∈ {1e-6, 1e-9}, linear H₄ and H₂: the multiset of allowed eigenvalues equals PySCF's
  D2h symmetry-adapted FCI spectrum in the HF irrep (Ag, S_z = 0) — same count, max |ΔE| < 1e-8 Ha.
- **G3 — the bead's prediction: no recorded number moves.** H₂ (0.74, 2.0), linear H₄, HeH⁺, LiH
  (12 qubits), LiH CAS(2,5), N₂ CAS(6,6): `reachable_mask == (pops > tol)` for tol ∈ {1e-8, 1e-10}.
  One gained or lost level kills it (and becomes the finding).
- **G4 — the bound is load-bearing.** sq-H₄ a ∈ {1.20, 1.40} (broken RHF) and LiH: some bit has
  β ≥ 1 and is excluded; the decision then equals the population cut; and at a = 1.20 the naive
  labels (every bit, no bound) would veto a level with population > 1e-3.
- **G5 — fallback is the old behaviour.** No `build_args` (integrals path): only N/S_z is used
  and the decision equals the population cut. A `build_args` that does not reproduce the operator
  (a different geometry swapped in) disables the spatial bits; a constant shift does not.
- **G6 — never near the boundary.** On every G1–G4 system, eigenvectors with population > 1e-10
  have sector weight within 1e-3 of 0 or 1.
- **G7 — the sites consume it.** At the two witnesses, `reachable_eigenpairs`,
  `hf_overlap_certificate.exact_reachable_overlap`, `hf_overlap_subspace.exact_hf_subspace_overlap`
  and the centered frames of `odmd`, `msd`, `trotter_odmd`, `device_odmd`,
  `trotter_resolution_floor` are built from the decision's reachable set (μ, τ, overlaps), not the
  population cut's.

## 6. Implementation plan (test-first)

1. This spec, then the gate file (fails: functions missing).
2. `build_args` on `MolecularHamiltonian`; the decision in `reachability.py`.
3. Rewire: `reachable_eigenpairs` (→ certified_gaps, hf_overlap_certificate, certified_dipole,
   certified_noise, excited_bounds), `hf_overlap_subspace`, `odmd`, `msd`, `trotter_odmd`,
   `device_odmd`, `trotter_resolution_floor` — each through `_dense_hf_projection` (which also
   removes their macOS ZHEEVD crash at 12 qubits).
4. Re-run every affected gate; record any moved number below.

## 7. Out of scope — sites deliberately NOT rewired, and why (part of the finding)

The symmetry argument holds where a site thresholds HF populations **in the exact-H eigenbasis**
(exact evolution conserves them). It does not transfer to:

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

## 8. Caveats and risks

- **R1 — degenerate clusters.** The majority rule is per eigenvector. An HF-sector level exactly
  degenerate with other-sector levels can be split by `eigh` so that no vector has weight > 1/2.
  Not seen on any gated system (G6); spin multiplets do mix this way but carry zero HF population.
- **R2 — residual symmetry outside PySCF's abelian frame.** At broken-RHF geometries a level can
  still carry convergence-dependent population (sq-H₄ 1.20: 1.5e-9 → 6.2e-9 between `conv_tol` 1e-9
  and 1e-6) that the D2h bits cannot see; there the 1e-8/1e-10 split still bites.
- **R3 — cost.** One extra SCF + JW map per decision (validation scale, like the dense `eigh`).
- **R4 — Linux numbers.** Where Linux's default-`conv_tol` SCF leaves residue (sq-H₄ 1.1 Å,
  5.07e-10), `exact_reachable_overlap` now returns the Ag overlap, so
  `test_reachability_tolerance_spec` G2 and G5 (`constant_governs_the_certified_reference`), which
  pin the *bug*, will fail on Linux. That file is owned elsewhere; recorded as a follow-up.

## 9. Results (macOS 27, Apple M3, 2 BLAS/OMP threads, 2026-10-02)

- **This gate file:** 33 passed in 27 s (G1–G7).
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

## 10. Deliverables

- `reachability.py` — the decision; `molecular_hamiltonian.py` — `build_args` only.
- Rewired sites listed in §6.3; `tests/test_eigenstate_reachability_spec.py` — G1–G7.
