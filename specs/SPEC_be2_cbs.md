# SPEC: Be2 well depth toward experiment -- CBS extrapolation + core-valence correlation

**Status:** CLOSED -- gates G1-G4 PASS (2026-07-12). `be2_cbs.py` merged. **Finding: the backlog's
original tolerance is NOT met** -- CBS(TZ/QZ) + core-valence NEVPT2 on the validated CAS(4,8)
reference gives a genuine bound well at the right general location (unlike the small-basis
baseline, which has none), but underbinds by roughly half and overshoots Re by ~0.15 A. **MEASURED:
Re=2.58 A, De=469 cm^-1** (baseline was 625 cm^-1 from experiment; this composition is 461 cm^-1
from experiment -- real progress, not closure). Gate G4 records the residual and its attributed
cause instead of forcing the pass.

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

BACKLOG.md claim: core-valence correlation + a cc-pVXZ->CBS extrapolation moves the well depth of
Be2 from the frozen-core CAS(4,8)/cc-pVDZ FCI baseline (`study_be2.py`, ~305 cm^-1) toward the
experimental D_e = 929.7 cm^-1 at R_e = 2.45 A (Merritt, Bondybey & Heaven, *Science* 2009). The
claim is false if the CBS+CV composition does not move measurably toward experiment, or if it
cannot even reproduce a well at the physically correct bond length.

## 2. Background and honest framing

- **Prior art / reference.** Merritt, Bondybey & Heaven, *Science* 323, 1671 (2009): D_e =
  929.7(20) cm^-1, R_e = 2.4498(9) A, resolving a decades-long theory/experiment controversy (HF
  predicts Be2 essentially unbound; the bond is created almost entirely by 2s-2p near-degenerate
  correlation). `study_be2.py` already reproduces the *qualitative* HF-unbound / correlated-bound
  contrast in this repo; this spec asks whether a *quantitative* CBS+core-valence composition of
  already-validated primitives closes the gap.
- **What we can claim if gates pass.** CASCI(4,8) (the validated static-correlation reference) +
  NEVPT2 dynamic correlation (unfrozen core, so core-valence/core-core correlation is included
  automatically), CBS-extrapolated cc-pVTZ/QZ -> QZ, gives a well at a physically sensible bond
  length that is closer to experiment than the frozen-core small-basis baseline.
- **What we cannot claim (stated up front).**
  1. **Quantitative agreement.** The measured result underbinds by ~half and misplaces R_e by
     ~0.15 A (G4) -- this is a *reproduction attempt that falls short*, recorded honestly, not a
     929.7 cm^-1 match.
  2. **CASSCF was tried and rejected, not skipped for convenience.** Orbital-optimized CASSCF(4,8)
     on Be2 is numerically unstable in this basis range -- `mc.converged = False` at multiple
     (basis, R) points during development, and the *converged* solutions moved the well depth
     non-monotonically and non-physically across DZ->TZ->QZ (792 -> 251 cm^-1 in one basis step,
     with `RankWarning`-flagged fits on top). CASCI on fixed canonical HF orbitals (the
     `study_be2.py`-validated choice) converges everywhere it's used here; the price is no orbital
     relaxation, which is part of why G4 underbinds.
  3. **The original backlog baseline number is itself an artifact.** Reproducing it (`fci_energy`
     on the exact `study_be2.py` grid) confirms ~304.8 cm^-1 -- but the discrete minimum sits at
     R=4.5 A, far past any real Be-Be bond. The frozen-core CAS(4,8)/cc-pVDZ curve has **no real
     well near 2.45 A at all**; "~305 cm^-1" was never a quantitative starting point to close a
     100 cm^-1 gap from.
  4. Minimal basis for the CBS pair (TZ/QZ only, the standard minimal Helgaker pair); DZ is
     excluded from the extrapolation (it does not show a well near R_e at all -- see G1).
  5. Isolated-dimer, fixed geometry scan (no vibrational/rotational correction to D_e).

## 3. Approach

`casci_nevpt2_point(R, basis)`: CASCI(4,8) on canonical RHF orbitals (unfrozen core), then
`pyscf.mrpt.NEVPT(mc).kernel()` for the dynamic correlation. Across `R` this gives a curve at each
of cc-pVDZ/TZ/QZ. `cbs_extrapolate_correlation` applies the standard two-point Helgaker form
`E_corr(X) = E_CBS + B/X^3` (X=3,4 for TZ,QZ) to the NEVPT2 correlation energy at each R; the
CASCI reference energy is taken at QZ (near-saturated: DZ->TZ->QZ CASCI-alone shifts by ~2 mHa
then ~0.7 mHa at R=2.45, an order of magnitude below the correlation-energy basis dependence).
`quadratic_well` fits a parabola through 3 bracketing R points to extract (Re, De); the well-depth
reference (dissociated limit) is a separate point at R=8.0 A (validated flat to within ~20-40
cm^-1 of R=6.0 A across all three bases -- the flatness check is G1).

## 4. Public interface

```
be2_cbs.casci_nevpt2_point(R, basis, cas_electrons=4, cas_orbitals=8) -> Be2Point
be2_cbs.cbs_extrapolate_correlation(x_lo, e_lo, x_hi, e_hi) -> float
be2_cbs.cbs_point(R, lo_basis="ccpvtz", hi_basis="ccpvqz", ...) -> float
be2_cbs.quadratic_well(Rs, Es, r_asymptote, asymptote_energy) -> (Re, De)
```

(`hybrid_quantum_solver.dmrg_reference.fci_energy` reused for the baseline reproduction in G3/G4.)

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_be2_cbs_spec.py` (test-first). PySCF only (no block2); a deliberately small
R-grid (R in {2.4, 2.45, 2.6} for the well window + {6.0, 8.0} for the asymptote, at DZ/TZ/QZ as
needed by each gate) keeps the gate to ~1-2 min -- QZ dominates at ~13 s/point.

- **G1 -- the small-basis curve has no real well near Re; TZ/QZ do; the asymptote is converged
  enough to use.** `cc-pVDZ` NEVPT2 energy at R=2.45 is *higher* (less bound) than at R=8.0 -- no
  minimum near the physical bond length, confirming the "~305 cm^-1" baseline is a far-R artifact,
  not a near-Re well -- while `cc-pVTZ` and `cc-pVQZ` are both *lower* at R=2.45 than at R=8.0
  (they do show a real well there, unlike DZ). At both TZ and QZ, R=6.0 and R=8.0 agree to
  `< 150 cm^-1` (MEASURED 110/129 cm^-1 -- an order of magnitude below the well depth, flat
  enough to define De from, though not fully converged to the true R->infinity limit).
- **G2 -- correlation energy grows monotonically with basis (the CBS direction is sane).**
  `|E_corr(TZ)| < |E_corr(QZ)|` at R=2.45 for the CASCI+NEVPT2 correlation energy (basis
  enlargement recovers more dynamic correlation, the expected direction; it is NOT asserted to be
  smooth/monotonic across DZ too -- DZ is excluded from the CBS fit per §2).
- **G3 -- the composition moves toward experiment (definition of "moves ... toward", the literal
  backlog claim).** `|De_cbs - 929.7| < |De_baseline - 929.7|`, where `De_baseline` is the
  frozen-core CAS(4,8)/cc-pVDZ FCI well depth reproduced from `study_be2.py`'s own grid (~304.8
  cm^-1, `fci_energy`) and `De_cbs` is the CBS(TZ/QZ)+CV well depth from a 3-point quadratic fit
  around R=2.4-2.6.
- **G4 -- the honest boundary (definition of "done", NOT the original 100 cm^-1/0.1 A backlog
  gate, which does not hold).** The measured CBS+CV well is pinned as a regression:
  `Re` in **(2.50, 2.70) A** and `De_cbs` in **(400, 550) cm^-1** -- i.e. explicitly OUTSIDE the
  original `|D_e-930|<100`/`|Re-2.45|<0.1` bounds. This gate fails loudly if a future change
  silently regresses the number *or* silently "fixes" it back into the original tolerance without
  a recorded reason (either direction is worth knowing about).

> Definition of done: **G4**, which is the recorded finding -- the original backlog gate is
> falsified, not satisfied, and that falsification (with cause) is what this spec closes on.

## 6. Implementation plan (test-first)

1. `tests/test_be2_cbs_spec.py` encodes G1-G4 against `be2_cbs.py` (initially failing / absent).
2. `be2_cbs.py`: `casci_nevpt2_point`, `cbs_extrapolate_correlation`, `cbs_point`,
   `quadratic_well`, reusing `study_be2.py`'s CASCI(4,8) active-space convention and
   `dmrg_reference.fci_energy` for the baseline.
3. Iterate to green; if G3/G4's actual numbers differ from the investigation above, update the
   pinned ranges and the prose here to match reality (not the reverse).

## 7. Out of scope

- CASSCF orbital relaxation (tried, numerically unstable here -- §2.2).
- Larger/Rydberg-augmented active spaces, higher-order MRPT (NEVPT3), CCSD(T)-F12, or MRCI+Q --
  any of which the literature needed to actually close the last ~450 cm^-1 / 0.15 A. A genuine
  follow-up, not attempted here.
- Vibrational zero-point / rotational corrections to D_e.
- A three-point (DZ/TZ/QZ) CBS fit -- DZ is excluded per G1's finding.

## 8. Caveats and risks

- **R1 (materialized, not hypothetical):** CASSCF reoptimization was the first approach tried and
  was numerically unstable (non-monotonic well depth across basis, unconverged CASSCF at several
  points) -- switched to fixed-orbital CASCI, which is why no orbital relaxation is captured.
- **R2:** the 3-point quadratic well fit (gate-cheap) is a local approximation of the true
  minimum; cross-checked in `be2_cbs.py.__main__` against a 13-point/quartic-fit curve and agrees
  to within ~4 cm^-1 / 0.02 A (informal, not gated -- the gate uses the cheap 3-point version).
  Honest limitation: reproduces a known result, and reproduces it incompletely.

## 9. Deliverables

- `be2_cbs.py` -- CASCI+NEVPT2 points, CBS extrapolation, quadratic well fit, `__main__` driver
  (writes `data/be2_cbs_curve.csv`).
- `tests/test_be2_cbs_spec.py` -- gates G1-G4.
- BACKLOG.md entry moved to Done with the honest numbers (not the original target) recorded.
