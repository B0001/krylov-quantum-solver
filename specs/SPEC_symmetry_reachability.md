# SPEC: the symmetry filter is available exactly when the artifact is present

**Status:** IMPLEMENTED once gates green. Closes the backlog hypothesis *"A symmetry/sector-aware
reachability test — the replacement for the constant"*, filed by
[`SPEC_reachability_tolerance.md`](SPEC_reachability_tolerance.md) §2b.

---

## 1. Goal

`SPEC_reachability_tolerance` §2b proved that **no fixed amplitude threshold** separates a physical
HF overlap from an SCF convergence residue on a symmetry-forbidden state: at square H₄ a = 1.1 Å the
level admitted at `tol=1e-10` is B1g with FCI coefficient exactly zero, and at a = 1.190 Å the same
residue reaches 1.4e-8, clearing the looser threshold too. The proposed replacement was a spatial
**irrep filter** — keep the eigenstates whose symmetry matches the HF determinant, no threshold.

When filing it I named its likely killer: *"broken-symmetry SCF solutions are exactly the case PySCF
cannot symmetry-adapt, and the 0.47-amplitude geometries show the driver lands there routinely."*

**The claim this spec tests: that killer is self-consistent rather than fatal — the filter is
available exactly when the artifact is present.**

## 2. Background and honest framing

Measured across the square-H₄ family, 9 geometries, **9/9**:

| regime | p₀ | irrep filter | why |
|---|---|---|---|
| symmetric SCF (a = 1.05, 1.10, 1.35) | 1.9e-29 … 6.1e-9 | **available** | forbidden level carries only SCF residue — the artifact, and the filter can remove it |
| broken-symmetry SCF (six geometries) | 0.44 … 0.48 | **unavailable** | reference is genuinely lower (~0.08 Ha) and genuinely overlaps that level — nothing to remove |

The two regimes are separated by **seven orders of magnitude** with nothing in between (G2), so this
is a clean dichotomy, not a graded continuum where a filter would have to make a judgement call.

**What we can claim:** where the amplitude threshold fails, the symmetry filter works; where the
symmetry filter refuses, the amplitude is unambiguous anyway. So a reachability test built on
"symmetry when available, amplitude otherwise" has no gap — on this family.

**What we cannot claim.** This is **one family** (square H₄) plus four ordinary systems as controls.
It is an existence-and-mechanism result, not a proof that the dichotomy holds universally; a system
with a genuinely intermediate p₀ *and* a symmetry-broken reference would break it, and this spec
does not rule that out. It also does **not** rewire the ~13 threshold sites in the repo — that is a
separate change with its own blast radius, and the module is deliberately standalone until then.

**Not claimed:** that RHF symmetry breaking on square H₄ is a defect. It is textbook behaviour for a
strongly-correlated square, and the broken solution is variationally *better*. The finding is about
what a reachability test may conclude, not about the chemistry.

## 3. Approach

Two questions, decided separately and for different reasons:

1. **Is the reference symmetry-broken?** Compare unsymmetrized RHF against symmetry-enforced RHF by
   **energy**. Broken ⟺ the free solve found a strictly lower determinant.
2. **Can the orbitals be labelled?** Only asked when (1) says not broken.

**The distinction LiH forced.** An earlier version decided availability by *whether labelling threw*,
and LiH failed — not because anything was broken (`dE ≈ 0`) but because the free solve returns an
arbitrary rotation within the degenerate E1x/E1y pair, which cannot be labelled column-by-column.
Same determinant, different basis. Deciding by energy separates *genuine breaking* (refuse) from
*degenerate-block rotation* (relabel via the symmetry-adapted solve). Without that split the filter
refuses ordinary closed-shell molecules and is useless in practice — G4 exists to catch exactly that.

**Reference:** symmetry-enforced RHF energy (for the break test) and `pyscf.symm.label_orb_symm`
(for labelling); the artifact regime is identified by the dense HF population spectrum.

## 4. Public interface

```
reachability.scf_symmetry_status(atom, basis, conv_tol)   -> (is_broken: bool, dE: float)
reachability.hf_orbital_irreps(atom, basis, conv_tol)     -> list[str] | None
reachability.symmetry_filter_available(atom, ...)         -> bool
reachability.hf_population_spectrum(mh)                   -> np.ndarray     # dense, validation only
reachability.SYMMETRY_BREAK_TOL   = 1e-6
reachability.TIGHT_SCF_CONV_TOL   = 1e-13
```

## 5. Acceptance criteria (validation gates)

`tests/test_symmetry_reachability_spec.py`.

- **G1 — the correlation (DEFINITION OF DONE).** Across the 9-geometry square-H₄ sweep, filter
  availability equals artifact presence. **One mismatched geometry kills it.**
- **G2 — the regimes are genuinely distinct.** Artifact geometries have p₀ < 1e-7; broken-symmetry
  ones have p₀ > 0.4. Killed if anything lands between — the dichotomy would be a continuum.
- **G3 — the mechanism.** Unavailable ⟺ the free solve found a solution ≥ 0.05 Ha lower (a = 1.20,
  1.40); available ⟺ the two solves agree to within `SYMMETRY_BREAK_TOL` (a = 1.05, 1.10, 1.35).
- **G4 — not vacuous.** The filter must accept H₂ eq/stretched, linear H₄, and **LiH** — the systems
  the specs actually gate on. This is the gate that failed first and forced §3's distinction.
- **G5 — scope.** On H₂ and linear H₄ every population below the physical ground state is < 1e-20,
  so no threshold choice matters there. Bounds the whole reachability problem to systems like
  square H₄.
- **G6 — constants pinned.** `TIGHT_SCF_CONV_TOL == 1e-13` (the value at which the residue collapses)
  and `SYMMETRY_BREAK_TOL == 1e-6`.

## 6. Implementation plan (test-first)

1. `tests/test_symmetry_reachability_spec.py` encoding G1–G6 (initially failing).
2. `reachability.py` — the decision procedure, standalone.
3. `make gates`.

## 7. Out of scope

- **Rewiring the ~13 reachability-threshold sites.** Needs its own blast-radius analysis; filed
  separately. This module is the decision procedure those sites will consume, nothing more.
- Systems beyond the square-H₄ family and the four controls.
- Open-shell / UHF references, and any system where the HF determinant is not a single closed-shell
  RHF solution.

## 8. Caveats and risks

- **R1 — one family.** G1's 9/9 is existence and mechanism, not universality. A system with an
  intermediate p₀ *and* a broken reference would falsify the dichotomy; none was found, none is
  ruled out.
- **R2 — `symmetry=True` requires PySCF to detect the point group.** Systems it cannot symmetry-adapt
  at all fall through to `None`, which is a refusal, not an error. Callers must handle `None`.
- **R3 — the break test costs a second SCF.** Negligible at these sizes, not free at scale.
- The tight `conv_tol=1e-13` default matters: at PySCFDriver's default 1e-9 the residue is 5e-10 and
  the artifact is *present*, which is the whole reason this spec exists. Callers reusing a loosely
  converged reference get the old behaviour.

## 9. Deliverables

- `reachability.py` — the decision procedure.
- `tests/test_symmetry_reachability_spec.py` — G1–G6.
