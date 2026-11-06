# SPEC: the symmetry filter is available exactly when the artifact is present

**Status:** IMPLEMENTED once gates green. Closes the backlog hypothesis *"A symmetry/sector-aware
reachability test — the replacement for the constant"*, filed by
[`SPEC_reachability_tolerance.md`](SPEC_reachability_tolerance.md) §2b.

**Revised 2026-11-05 (chem-8ls, `SPEC_scf_stability_rule.md`):** §10: the artifact regime was the SCF stopping on a saddle point; G2 and G3 now assert the stable branch at every side.

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

## 10. chem-8ls revision (2026-11-05): the artifact regime was the SCF stopping on a saddle

**What changed.** Every builder now follows internal SCF instabilities to a stable minimum, and `scf_symmetry_status` reads the builder's own SCF (`SPEC_scf_stability_rule.md`).
- The D2h-symmetric RHF is a saddle point at every side of this sweep. Its lowest orbital-Hessian eigenvalue is −0.53 to −0.67 Ha.
- So the rule puts every side on the broken-symmetry minimum, on macOS arm64 and Linux x86_64 alike.
- The "symmetric SCF" rows of §2 were sides where the free SCF happened to stop on that saddle. Which sides those were depended on the platform: 1.05/1.10/1.35 in this sweep on macOS, 1.15/1.20/1.25/1.30 on Linux (`SPEC_scf_stability_rule` §9).

**Gates.**
- **G1** still holds 9/9, now as one regime: no side shows the artifact, and the filter is available at none.
- **G2 is re-pinned.** It was "both regimes occur, seven orders apart". It is now "every side has 0.4 < p₀ < 0.5, strictly decreasing in a (no branch jump)". The two-regime claim is retired: its small-p₀ regime only ever existed on the saddle.
- **G3, the "unavailable" half,** now runs at all nine sides; it was 1.20 and 1.40 only.
- **G3, the "available" half (1.05, 1.10, 1.35), is retired.** Those were the macOS saddle stops.
- **G4–G6 are unchanged.** The filter is still available on H₂, linear H₄ and LiH, whose RHF is stable, and the rule leaves those builds bit-identical (`SPEC_scf_stability_rule` G3b).

**Old values next to new.** The old values are from §2 (macOS, before). The new values were measured 2026-11-05 on macOS 27.2 arm64 and on Linux x86_64 (the `claude-chem` image), with pyscf 2.13.1 and qiskit-nature 0.8.0 on both, and are identical to 8 digits.

| a (Å) | 1.00 | 1.05 | 1.10 | 1.15 | 1.20 | 1.25 | 1.30 | 1.35 | 1.40 |
|---|---|---|---|---|---|---|---|---|---|
| p₀ before (§2) | 0.44–0.48 | ≤ 6.1e-9 | ≤ 6.1e-9 | 0.44–0.48 | 0.44–0.48 | 0.44–0.48 | 0.44–0.48 | ≤ 6.1e-9 | 0.44–0.48 |
| regime before | broken | symmetric | symmetric | broken | broken | broken | broken | symmetric | broken |
| p₀ now | 0.4759 | 0.4727 | 0.4692 | 0.4654 | 0.4612 | 0.4567 | 0.4517 | 0.4463 | 0.4405 |
| regime now | broken | broken | broken | broken | broken | broken | broken | broken | broken |

**What survives of the finding.** "Available exactly when the artifact is present" is still what G1 checks. But on square H₄ the builder no longer makes the artifact-present case. A reference that sits on the saddle (for example a symmetry-adapted SCF) would still carry a forbidden level, and could still be labelled and filtered. `SPEC_eigenstate_reachability`'s residue witness is such a reference, built on purpose.

**Gate run.** See `SPEC_scf_stability_rule` §9.
