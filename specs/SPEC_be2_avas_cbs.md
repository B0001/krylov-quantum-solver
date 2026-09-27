# SPEC: Be2 -- a basis-consistent (AVAS) active space, retested against the chem-mom CBS attribution

**Status:** CLOSED -- gates G1-G4 PASS (2026-09-26). `be2_avas_cbs.py` merged. **Finding: the fix
works exactly where it targets (the active space IS now basis-consistent, G1/G2), but that alone
does NOT rescue the "method, not basis" CBS attribution (G3 still fails) -- it shrinks the QZ->5Z
"CBS" swing from +1224 cm^-1 (chem-mom, canonical-orbital CAS) to +530 cm^-1 (AVAS CAS), still >>
the pre-registered 20 cm^-1 confirmation bar. The single-basis wells are still non-monotone in X,
but the decomposition (G4) shows why: cc-pVTZ is the outlier (its CASCI-only and NEVPT2-only pieces
both differ sharply from QZ/5Z, which agree much better with each other than either does with TZ) --
consistent with TZ being pre-asymptotic for this CAS(4,8), not with a basis-inconsistent space.**

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

`specs/SPEC_be2_cbs.md` G5 (bead chem-mom) found that adding cc-pV5Z to the Be2 CBS(TZ/QZ)+NEVPT2
composition moves the extrapolated well depth by +1224 cm^-1 against a pre-registered ~20 cm^-1
confirmation bar, and diagnosed the likely cause: CASCI(4,8) built from the 8 lowest-energy
*canonical HF virtuals* is a different physical active space at every basis (outermost active
orbital rms spread 7.1 -> 8.3 -> 9.4 bohr, orbital energy +0.150 -> +0.074 Ha, TZ->QZ->5Z). This
spec's claim: selecting the active space by AVAS-projecting onto a **fixed** Be 2s/2p reference
(independent of the target basis) removes that confound, so the TZ/QZ/5Z wells become basis-
comparable and the CBS(QZ/5Z) vs CBS(TZ/QZ) attribution test (chem-mom's G5) can be rerun on a
like-for-like active space. The claim is false if AVAS does not actually stabilize the active-space
character across bases, or if the CBS attribution still does not survive 5Z even once it does.

## 2. Background and honest framing

- **Prior art.** `specs/SPEC_be2_cbs.md` (G1-G4: CASCI(4,8)+NEVPT2, CBS(TZ/QZ), underbinds
  experiment by ~half); its G5/bead chem-mom (`be2_cbs_5z.py`, `results/be2_cbs_5z/`): adding 5Z
  falsifies the TZ/QZ "method, not basis" attribution, diagnosed to the active-space drift above.
  AVAS is Sayfutyarova, Sun, Chan & Knizia, *J. Chem. Theory Comput.* 13, 4063 (2017),
  arXiv:1701.07862 -- selects active orbitals by projecting the HF occupied+virtual space onto a
  set of AO labels in a reference minimal basis, keeping whatever has projector eigenvalue above a
  threshold.
- **What we can claim if gates pass.** The AVAS(Be 2s/2p, `minao='ano'`)-selected CAS(4,8) is the
  *same* active space (same size, same orbital character/spread) at cc-pVTZ, QZ and 5Z, unlike the
  canonical-orbital CAS(4,8) it replaces (G1, G2) -- this directly fixes the diagnosed confound.
- **What we cannot claim (stated up front).**
  1. Fixing the active-space confound does not, by itself, produce a basis-converged CBS well --
     the CASCI reference energy and the NEVPT2 dynamic correlation on top of it can still each have
     their own basis-set dependence even on an unchanging active space. **Materialized (G3):** it
     doesn't -- the QZ/5Z shift is still 27x the confirmation bar, just smaller than before.
  2. AVAS orbitals are not canonical HF eigenstates (they're an AO-projected rotation within
     occ+virt space), so "orbital energy" reported here is `<phi|F|phi>` (a Fock expectation value),
     not an eigenvalue -- a like-for-like comparison to the old canonical-orbital diagnostic, not a
     claim that these are stationary states of anything.
  3. `minao='minao'` (pyscf's AVAS default) is UNUSABLE for this system: Be's minimal ground-state
     basis is 1s^2 2s^2 with no 2p shell at all, so the default reference silently drops all 2p
     character from the active space (verified interactively: ncas=2, not 8). `minao='ano'` (the
     full ANO basis, which includes shells unoccupied in the free atom) is required and is asserted
     at runtime (G1), not assumed.
  4. This does not revisit CASSCF (still numerically unstable here, per `SPEC_be2_cbs.md` R1) or
     try a larger/different active space -- same CAS(4,8), same NEVPT2, only the orbital-selection
     rule changes.

## 3. Approach

`avas_casci_nevpt2_point(R, basis)`: RHF -> `pyscf.mcscf.avas.avas(mf, ('Be 2s','Be 2p'),
minao='ano', threshold=0.2, canonicalize=False)` to get `(ncas, nelecas, mo)` -> assert
`(ncas, nelecas) == (8, 4)` (hard failure otherwise, not a silent fallback) -> `CASCI(mf, 8,
4).kernel(mo)` -> `NEVPT(mc).kernel()` for the dynamic correlation (core unfrozen, same as
`be2_cbs.py`). Reuses `Be2Point`, `HA2CM`, `BASIS_CARDINAL`, `cbs_extrapolate_correlation` from
`be2_cbs.py` unchanged; only the orbital-selection step differs from
`be2_cbs.casci_nevpt2_point`. Reference: the same experimental D_e = 929.7(20) cm^-1, R_e =
2.4498(9) A (Merritt, Bondybey & Heaven, *Science* 2009) `SPEC_be2_cbs.md` checks against, plus
chem-mom's own pre-registered ~20 cm^-1 CBS-consistency bar as the reference for "did the confound
go away."

## 4. Public interface

```
be2_avas_cbs.avas_active_space(mf, aolabels=("Be 2s","Be 2p"), minao="ano", threshold=0.2,
                               cas_electrons=4, cas_orbitals=8) -> (ncas, nelecas, mo_coeff)
be2_avas_cbs.avas_casci_nevpt2_point(R, basis) -> Be2Point            # reuses be2_cbs.Be2Point
be2_avas_cbs.active_orbital_diagnostics(basis, R) -> (fock_energies, rms_spreads_bohr)
be2_avas_cbs.quartic_well(curve) / cbs_curve(pts, lo, hi) / summarize(pts)   # driver internals
```

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_be2_avas_cbs_spec.py`. Cheap grid `R in {2.4, 2.45, 2.6}` (well window) +
`{8.0}` (asymptote) at TZ/QZ/5Z, matching `test_be2_cbs_spec.py`'s convention -- ~1 min total
(5Z ~10 s/point dominates).

- **G1 -- the active space has the expected shape everywhere it's used.** `avas_active_space`
  returns `ncas=8, nelecas=4` at every (basis, R) checked (TZ/QZ/5Z x the gate grid); a mismatch is
  a hard `RuntimeError` from the production code itself, not just a test assertion -- this is the
  literal fix for the failure mode diagnosed in chem-mom (silently different-sized/shaped spaces).
- **G2 -- the active-space character is stable across basis (the actual fix, checked directly).**
  Per-orbital rms spread (sorted, since degenerate p-pairs can swap order) agrees between cc-pVTZ
  and cc-pV5Z to **within 1.0 bohr** at both R=2.5 A and R=8.0 A -- MEASURED max diff 0.004 bohr
  (R=2.5) / 0.0003 bohr (R=8.0), i.e. two orders of magnitude inside the bar, vs. the canonical-
  orbital baseline's 2.3 bohr (7.1->9.4) drift on the outermost orbital. This is the gate that
  would fail loudly if AVAS were not actually basis-consistent.
- **G3 -- CBS attribution, retested on the fixed active space (chem-mom's own pre-registered bar,
  fixed before this spec's numbers existed).** `|De_CBS(QZ/5Z) - De_CBS(TZ/QZ)| <= 20 cm^-1` ->
  CONFIRMED; else NOT CONFIRMED, record the residual. MEASURED (cheap 4-point grid):
  CBS(TZ/QZ) De=55.7 cm^-1 @ Re=2.60 A, CBS(QZ/5Z) De=592.3 cm^-1 @ Re=2.53 A, shift **+536.6
  cm^-1** -- NOT CONFIRMED. Fixing G1/G2 shrank the swing (chem-mom's canonical-orbital version:
  +1224 cm^-1) but did not close it. Gate pins the measured ranges as a regression.
- **G4 -- single-basis De monotone in X, or the failure mode identified (this spec's own
  criterion, since G3 already falsifies the "basis-consistent active space is sufficient"
  hypothesis).** MEASURED single-basis wells (cheap grid): TZ 730.6, QZ 304.7, 5Z 446.2 cm^-1 --
  still non-monotone. Decomposition (`De = De_casci_only + De_from_nevpt2`) shows QZ and 5Z agree
  with each other far better than either agrees with TZ (`|De_QZ - De_5Z| < |De_TZ - De_QZ|`,
  MEASURED 141.5 < 425.9 cm^-1) on BOTH pieces -- consistent with cc-pVTZ being pre-asymptotic for
  this CAS(4,8)+NEVPT2 recipe (a known small-basis artifact), not with the active space still being
  basis-inconsistent (which G1/G2 rule out directly). Gate checks this inequality plus pins the
  three single-basis wells as a regression.

> Definition of done: **G3**, which is the literal retest of chem-mom's pre-registered CBS
> criterion on the fixed active space -- it is the honest finding this spec closes on (NOT
> CONFIRMED, smaller residual, different diagnosed cause).

## 6. Implementation plan (test-first)

1. `tests/test_be2_avas_cbs_spec.py` encodes G1-G4 against `be2_avas_cbs.py` (initially failing).
2. `be2_avas_cbs.py`: `avas_active_space`, `avas_casci_nevpt2_point`,
   `active_orbital_diagnostics`, reusing `be2_cbs.py`'s `Be2Point`/`cbs_extrapolate_correlation`.
3. `make gates` on `test_be2_avas_cbs_spec.py` in isolation (PySCF only, no block2).

## 7. Out of scope

- Rerunning `study_be2.py`'s CAS(4,8)-on-canonical-orbitals FCI/Krylov/DMRG comparisons with AVAS
  orbitals -- that spec's claim is qualitative (HF unbound vs correlated bound) and untouched by
  this finding.
- A third active-space construction (valence natural orbitals from a correlated density) if AVAS
  itself turned out basis-inconsistent -- moot, since G1/G2 pass.
- Chasing the TZ outlier further (larger basis at low cardinal number, F12, etc.) -- a genuine
  follow-up if the CBS attribution is to be pursued further.
- Any change to `be2_cbs.py`'s own gates (G1-G5 stay pinned to the canonical-orbital numbers they
  measured; this is a new spec, not a revision of that one, since the canonical-orbital active
  space is a `study_be2.py`-validated convention this spec explicitly leaves alone elsewhere).

## 8. Caveats and risks

- **R1 (materialized):** `minao='minao'` (pyscf's AVAS default) silently drops Be 2p from the
  active space for this system (no 2p in Be's minimal ground-state basis) -- caught during
  development by checking `ncas` before writing any driver code, and enforced at runtime by G1's
  hard assertion, not just documented.
- **R2:** AVAS's `threshold=0.2` default was not scanned -- a different threshold could in
  principle select a different-sized space at one basis and not another; G1's runtime assertion
  would catch that (hard failure) rather than silently mis-comparing energies, which is why it's a
  gate and not just a diagnostic print.
- **R3:** the cheap 3-point quadratic well fit (R=2.4/2.45/2.6 vs 8.0) is used for the gates;
  cross-checked against the driver's 13-point wide-quartic fit (`results/be2_avas_cbs/wells.csv`)
  and agrees to within ~5 cm^-1 / 0.02 A -- informal, not gated, same convention as
  `SPEC_be2_cbs.md` R2.
- **R4:** this does not identify what actually causes the residual QZ/5Z CBS drift once the active
  space is fixed (basis-set superposition in NEVPT2 on an 8-orbital CAS, incompleteness of 2nd-order
  perturbation theory, or something else) -- G4's decomposition narrows it to "TZ is the outlier,
  QZ/5Z agree much better," but does not explain QZ/5Z's own residual +536.6 cm^-1 CBS shift, which
  remains open.

## 9. Deliverables

- `be2_avas_cbs.py` -- AVAS active-space selection + CASCI/NEVPT2 points, diagnostics, CBS/well
  driver (writes `results/be2_avas_cbs/{points.csv,wells.csv,diagnose.txt}`, tracked).
- `tests/test_be2_avas_cbs_spec.py` -- gates G1-G4.
- `specs/SPEC_be2_cbs.md` R3 updated to point here as the materialized follow-up (fix attempted,
  active-space confound removed, CBS attribution still not confirmed for a different reason).
- BACKLOG.md entry recording the honest finding.
