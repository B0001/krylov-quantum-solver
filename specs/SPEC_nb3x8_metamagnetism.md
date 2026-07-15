# SPEC: Nb3X8 metamagnetism — the field-driven singlet-to-triplet crossing

**Status:** IMPLEMENTED — gates G1–G4 green. Definition-of-done gate is G1.

---

## 1. Goal

Add a Zeeman term to the same exactly-diagonalizable Nb₃X₈ dimer used by
`nb3x8_susceptibility.py`/`odmd_spin.py`. Claim: the isolated dimer's ground state undergoes a
field-driven level crossing from the non-magnetic singlet to the fully polarized (Sz=+1) triplet
member at an exactly closed-form critical field h_c = J (the same interlayer exchange already
measured in `SPEC_odmd_spin`), reproduced by **direct diagonalization of the full field-augmented
Hamiltonian H0 − h·Ŝz** (not assumed from block symmetry). Falsifiable: any numeric crossing field
that disagrees with J, or a fractional/smeared magnetization step, kills the claim.

## 2. Background and honest framing

The singlet ground state and 3-fold-degenerate triplet at energy J (`n2_spectrum`,
`dimer_exchange_analytic`) are already gated in `SPEC_odmd_spin`/`SPEC_nb3x8_susceptibility`. A
uniform Zeeman field −h·Ŝz,tot couples only to the total spin projection; because [H0, Ŝz] = 0,
the fully-polarized Sz=±1 triplet members are single Slater determinants (the only two-same-spin
determinant in a 2-orbital space) whose energy is field-independent apart from the rigid −h·Sz
shift, while the Sz=0 singlet ground state is untouched by the field. This predicts an exact,
temperature-free metamagnetic step: M(h) = 0 for h < J, M(h) = 1 for h > J (units of ħ per dimer).

This is a genuinely new observable for this family (not covered by the existing χ(T), thermo, or
strain specs, all of which are zero-field) and turns the already-measured J into a **quantitative,
falsifiable magnet-field prediction** — expressed both in energy units and, via B = h/(gμ_B),
Tesla — with a feasibility read against real pulsed-field magnet records.

- **What we can claim:** the exact closed-form critical field h_c = J, verified against direct
  numerical diagonalization of the full field-augmented Hamiltonian (not merely asserted from
  symmetry); the resulting Tesla-scale prediction for the whole halide family; an honest
  feasibility comparison against the best documented pulsed-field records (non-destructive:
  100.75 T, Los Alamos, March 2012; destructive indoor flux-compression: 1200 T, U. Tokyo, 2018 —
  see Sources).
- **What we cannot claim:** an experimentally *observed* transition — this is a prediction from
  the isolated-dimer model, not a match to measured magnetization data (no such field-dependent
  data is cited by `arXiv:2501.10320` or the χ(T) experimental references already in this repo).
  `SPEC_nb3x8_magnetometry`'s own finding — that the isolated dimer *overpredicts* the ordering
  temperature by 2.3–5.3× — is a live reason to expect the real lattice's critical field to sit
  below this isolated-cluster number too; this spec does not attempt that renormalization.
  Nb₃F₈'s h_c (≈0.44 T) is not a real prediction: its J (0.051 meV) is already flagged in
  `SPEC_odmd_spin` as below the model's own neglected non-density-density terms — quote as "≈0",
  not a physical claim, same as the rest of the suite.

## 3. Approach

**References:**
1. Closed form: `dimer_exchange_analytic(U0, t, Us)` (already validated in `SPEC_odmd_spin`) is
   the exact singlet-triplet gap J = h_c, since the fully-polarized branch sits at the same
   zero-field energy as the other triplet members.
2. Direct numerical check: build the (U0,t,Us) dimer Hamiltonian via `dimer_cluster_integrals`
   (validated in `SPEC_nb3x8_gaps`/`SPEC_nb3x8_susceptibility`), add −h·Ŝz using the same Ŝz
   operator construction as `nb3x8_susceptibility.py`, and diagonalize the full N=2-sector matrix
   directly at each h (no assumed block structure) — bisect on the ground state's ⟨Ŝz⟩ to find the
   numeric crossing field, and compare the full ground-state energy curve against the closed-form
   prediction E(h) = min(E_singlet, E_triplet − h) across a grid.
3. Feasibility reference: two documented pulsed-field records (cited above), used only as a
   Tesla-scale sanity comparison, not a source of a physics gate.

## 4. Public interface

Reuses `dimer_cluster_integrals` (`nb3x8_gaps.py`), `dimer_exchange_analytic` (`odmd_spin.py`),
`n2_spectrum` (`nb3x8_susceptibility.py`). New module `nb3x8_metamagnetism.py`:

```
nb3x8_metamagnetism.G_MU_B                                      # g=2 Bohr magneton, meV/T
nb3x8_metamagnetism.zeeman_ground_state(U0, t, Us, h) -> (energy: float, sz: float)
nb3x8_metamagnetism.magnetization(U0, t, Us, h) -> float         # <Sz> of the ground state at field h (meV)
nb3x8_metamagnetism.critical_field_numeric(U0, t, Us, tol=1e-9) -> float   # bisection, meV
nb3x8_metamagnetism.critical_field_tesla(U0, t, Us) -> float     # B_c = J / (g*mu_B), Tesla
nb3x8_metamagnetism (CLI __main__)                               # family table + feasibility finding
```

## 5. Acceptance criteria (validation gates)

`tests/test_nb3x8_metamagnetism_spec.py` (test-first).

- **G1 — closed form == direct field-augmented ED (DEFINITION OF DONE).** For Nb₃Cl₈/Br₈/I₈:
  `|critical_field_numeric(U0,t,Us) - dimer_exchange_analytic(U0,t,Us)| < 1e-6` meV (bisection on
  the full Hamiltonian, no block-diagonal shortcut assumed).
- **G2 — rigid-shift identity.** Across a grid of h in [0, 1.5·J], the full field-augmented
  ground-state energy from direct diagonalization matches `min(E_singlet, E_triplet - h)`
  (referenced to the same absolute zero-field ground energy) to < 1e-8 meV, for all 4 materials —
  this is what makes the closed form more than a coincidence at one point.
- **G3 — clean magnetization step, no fractional plateau.** `magnetization(h_c - eps) < 0.01` and
  `magnetization(h_c + eps) > 0.99` (eps = 1e-3·J) for Cl/Br/I; and the first ionic singlet
  (`ionic_singlet_energy`, `SPEC_nb3x8_susceptibility`) sits above J for every material (so no
  other level intervenes before the crossing — the grid in G2 is a fair test of the two-level
  picture).
- **G4 — Tesla-scale honest feasibility boundary.** `critical_field_tesla` for Cl and Br exceeds
  the non-destructive pulsed-field record (100.75 T) but is below the destructive flux-compression
  record (1200 T); Nb₃I₈'s exceeds even the destructive record. Nb₃F₈ excluded from this gate (its
  J is below the model's own noise floor, per §2).

> G1 is the definition of done: the exact closed form derived from symmetry must survive an
> unassuming direct diagonalization of the full field-augmented matrix.

## 6. Implementation plan (test-first)

1. `tests/test_nb3x8_metamagnetism_spec.py` encoding G1–G4 (initially failing — no module).
2. `nb3x8_metamagnetism.py` composing the existing dimer/Ŝz primitives plus the Zeeman term.
3. `make gates` (own process; no block2/qiskit-aer conflict — this module only uses qiskit-nature
   mapping + numpy, same footprint as `nb3x8_susceptibility.py`).

## 7. Out of scope

- Any renormalization from the isolated cluster to the real lattice (the `SPEC_nb3x8_magnetometry`
  2.3–5.3× overcoupling finding already flags this as the likely direction of a correction).
- In-plane kagome exchange, anisotropy, or a full field-angle-dependent Hamiltonian — uniform
  field along the same Ŝz axis as `SPEC_odmd_spin` only.
- Finite-temperature magnetization (this is the T=0 ground-state crossing only; combining with the
  existing χ(T)/thermo machinery for M(h,T) is a natural follow-up, not attempted here).
- Verifying the disputed ~3000 T claims — only the two records with independent confirmation
  (Los Alamos 2012, U. Tokyo 2018) are used as the feasibility yardstick.

## 8. Caveats and risks

- **R1:** the g=2 (spin-only) assumption matches the existing susceptibility module's own
  documented convention (`nb3x8_susceptibility.py`'s EMU_PER_REDUCED derivation) — a real g-factor
  for this correlated Nb 4d system may differ, which would rescale B_c linearly. Not re-derived
  here; flagged, not corrected.
- **R2:** density-density-only Hamiltonian (same limitation as the rest of the Nb₃X₈ thread) — the
  same caveat that made Nb₃F₈'s J itself unreliable applies identically to its h_c.
- **R3:** isolated single dimer — no inter-dimer coupling, so this is a prediction for a
  hypothetical decoupled cluster, not the bulk (see §2's link to the magnetometry overcoupling
  finding).

## 9. Deliverables

- `nb3x8_metamagnetism.py` — new module (Zeeman-augmented dimer diagonalization).
- `tests/test_nb3x8_metamagnetism_spec.py` — gates G1–G4.
- `specs/BACKLOG.md` — recorded finding.
- Results summary (with the §2/§7 caveats) in the PR description.

## Sources (feasibility yardstick, §2/G4)

- Los Alamos National Laboratory, 100.75 T non-destructive pulsed magnet, March 2012 —
  <https://phys.org/news/2011-08-los-alamos-world-record-pulsed-magnetic.html>,
  <https://nationalmaglab.org/about-the-maglab/around-the-lab/meet-the-magnets/meet-the-100-tesla-pulsed-magnet/>.
- University of Tokyo, 1200 T indoor electromagnetic flux-compression record, 2018 —
  <https://www.sciencedaily.com/releases/2018/09/180917135933.htm>,
  <https://phys.org/news/2018-09-world-magnetic-field.html>.
