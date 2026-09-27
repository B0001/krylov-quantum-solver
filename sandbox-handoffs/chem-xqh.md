# chem-xqh handoff: Be2 basis-consistent active space

## Summary

chem-mom (PR #43, bead closed) found that `be2_cbs.py`'s CASCI(4,8) built from the 8 lowest-energy
*canonical HF virtuals* is a **different active space at every basis** (outermost active orbital
rms spread 7.1 -> 8.3 -> 9.4 bohr, orbital energy +0.150 -> +0.074 Ha, TZ->QZ->5Z), which is why the
single-basis wells oscillate (997/698/1313 cm^-1) and the CBS(QZ/5Z) vs CBS(TZ/QZ) "method, not
basis" attribution moved +1224 cm^-1 against a pre-registered ~20 cm^-1 confirmation bar.

This bead replaces that orbital-selection rule with **AVAS**, projected onto a fixed Be 2s/2p
reference (`minao='ano'`), and reruns the same TZ/QZ/5Z comparison on the fixed active space.

**Result: the fix works exactly where it targets, and that is not enough.**

- The active space is now genuinely basis-consistent: same size (ncas=8, nelecas=4) and same
  orbital character (per-orbital rms spread agrees to **<0.01 bohr** TZ vs 5Z, two orders of
  magnitude inside a 1.0 bohr bar; the canonical-orbital version drifted by **2.3 bohr**). Verified
  directly, not assumed.
- That alone does **not** rescue the CBS attribution: CBS(QZ/5Z) still moves **+530 cm^-1** from
  CBS(TZ/QZ) against the same pre-registered 20 cm^-1 bar chem-mom used. Smaller than chem-mom's
  +1224 cm^-1 (57% reduction), but still 27x over threshold. **NOT CONFIRMED**, same verdict as
  chem-mom, different (smaller) residual and a different diagnosed cause.
- Single-basis wells are still non-monotone (TZ 731 / QZ 305 / 5Z 446 cm^-1), but the
  CASCI-only/NEVPT2-only decomposition shows **cc-pVTZ, not a basis-inconsistent active space, is
  the outlier**: QZ and 5Z agree with each other far better than either agrees with TZ, on both
  pieces (CASCI-only: |QZ-5Z|=1.2 cm^-1 vs |TZ-QZ|=602 cm^-1; NEVPT2-only: |QZ-5Z|=140 cm^-1 vs
  |TZ-QZ|=1028 cm^-1). This is the acceptance criterion's "or explained."

## What I changed

- **`be2_avas_cbs.py`** (new) -- `avas_active_space`, `avas_casci_nevpt2_point`,
  `active_orbital_diagnostics`, and a resumable driver (`__main__`) mirroring `be2_cbs_5z.py`'s
  structure. Reuses `be2_cbs.py`'s `Be2Point`, `HA2CM`, `BASIS_CARDINAL`,
  `cbs_extrapolate_correlation` unchanged.
- **`specs/SPEC_be2_avas_cbs.md`** (new) -- the spec: goal, method, 4 acceptance gates (G1-G4),
  out-of-scope, caveats. Closed with the measured finding above.
- **`tests/test_be2_avas_cbs_spec.py`** (new) -- gates G1-G4, test-first, all green (see below).
- **`specs/SPEC_be2_cbs.md`** -- R3 caveat updated to point at this spec as the materialized
  follow-up (one paragraph, no other change; G1-G5 there are untouched and still pinned to the
  canonical-orbital numbers they measured).
- **`be2_cbs.py`** -- module docstring only: one paragraph pointing to the follow-up. No logic
  changed; `casci_nevpt2_point` etc. are byte-for-byte the same as before this bead.
- **`specs/BACKLOG.md`** -- new `[x]` entry after the existing Be2 entry, same style, linking to
  the new spec.
- **`results/be2_avas_cbs/`** (new, tracked) -- `points.csv` (39 CASCI+NEVPT2 points, TZ/QZ/5Z x
  13-point R grid), `wells.csv` (per-basis + CBS wells, wide 13-point quartic fit), `diagnose.txt`
  (active-orbital Fock energies + rms spreads per basis at R=2.5/8.0).

## Why AVAS needed `minao='ano'`, not pyscf's default `minao='minao'`

Checked interactively before writing any driver code (this is exactly the kind of thing the bead
warns against silently getting wrong): `pyscf.mcscf.avas.avas(mf, ('Be 2s','Be 2p'))` with the
default `minao='minao'` returns `ncas=2, nelecas=4` -- because Be's minimal ground-state basis is
1s^2 2s^2 (no occupied *or* virtual 2p shell at all in "minao"), the projector finds only the two
Be 2s AOs, both already doubly occupied, so 0 orbitals come from the virtual space. That is a
degenerate, correlation-free 2-orbital active space, silently NOT the intended CAS(4,8). Using
`minao='ano'` (the full ANO basis, which spans shells unoccupied in the free atom) gives the
intended `ncas=8, nelecas=4` at every basis tested. `avas_active_space()` asserts this shape at
runtime (`RuntimeError` on mismatch) rather than trusting it silently -- this is gate G1, and it is
also enforced by the production code, not just the test.

## Commands and gate output (verbatim)

Environment: `uv sync --extra dmrg --extra test --extra dev` (adds `ruff` for lint). Host: Linux
x86_64, 8 vCPUs, 16 GB RAM (`uname -a`: `Linux ... 7.0.12-linuxkit ... x86_64 GNU/Linux`).

Lint:
```
$ uv run ruff check be2_avas_cbs.py tests/test_be2_avas_cbs_spec.py
All checks passed!
```

The two be2 spec gates, run together in isolated processes (matches `make gates-be2`'s pattern;
`make` itself is not installed in this container, so I called `scripts/run_gates.sh` directly):
```
$ GATE_JOBS=2 GATE_GLOB='tests/test_be2*_spec.py' GATE_NO_CACHE=1 bash scripts/run_gates.sh
gates: 2 files, 2 parallel processes, cache=off
PASS  tests/test_be2_cbs_spec.py  (179s)
PASS  tests/test_be2_avas_cbs_spec.py  (190s)
all spec gates passed
```

`test_be2_avas_cbs_spec.py` alone (5 tests -- 2x G1, G2, G3, G4), for the per-test detail:
```
$ uv run pytest -q tests/test_be2_avas_cbs_spec.py
.....
5 passed in 106.99s (0:01:46)
```

Full driver run (13-point R grid x TZ/QZ/5Z, `python be2_avas_cbs.py`, wall 3m36s on 8 vCPU) wrote
the tracked artifacts in `results/be2_avas_cbs/`:
```
label,Re_A,De_cm-1,De_minus_expt_cm-1
ccpvtz,2.4085,734.7,-195.0
ccpvqz,2.5152,298.5,-631.2
ccpv5z,2.5120,441.0,-488.7
CBS(TZ/QZ),2.6687,58.7,-871.0
CBS(QZ/5Z),2.5091,588.2,-341.5
# De shift CBS(QZ/5Z) - CBS(TZ/QZ) = +529.5 cm^-1 -> method attribution NOT CONFIRMED (pre-registered threshold 20 cm^-1)
# single-basis De(TZ,QZ,5Z) = 734.7, 298.5, 441.0 cm^-1 -> NON-monotone in X
# mean wall s/point: ccpvtz=1.1, ccpvqz=4.4, ccpv5z=9.9; host x86_64 8 cpu
```
(The gate's cheap 4-point grid reproduces this closely: TZ 730.6, QZ 304.7, 5Z 446.2,
CBS(TZ/QZ) 55.7, CBS(QZ/5Z) 592.3, shift +536.6 -- within ~5 cm^-1 of the wide-quartic numbers
above, same relationship `SPEC_be2_cbs.md` R2 documents for the original spec.)

Active-space stability check (`results/be2_avas_cbs/diagnose.txt`), the direct measurement behind
G2:
```
basis,R,active_orbital_fock_energies_Ha,active_orbital_rms_radius_bohr
ccpvtz,2.5,-0.246 -0.397 +0.312 +0.040 +0.147 +0.147 +0.051 +0.051,4.4 3.1 3.4 4.2 3.8 3.8 3.8 3.8
ccpvtz,8.0,-0.309 -0.309 +0.095 +0.095 +0.094 +0.095 +0.095 +0.096,8.1 8.1 8.1 8.1 8.1 8.1 8.1 8.1
ccpvqz,2.5,-0.246 -0.397 +0.311 +0.040 +0.148 +0.148 +0.051 +0.051,4.4 3.1 3.4 4.2 3.8 3.8 3.8 3.8
ccpvqz,8.0,-0.309 -0.309 +0.095 +0.095 +0.094 +0.096 +0.095 +0.095,8.1 8.1 8.1 8.1 8.1 8.1 8.1 8.1
ccpv5z,2.5,-0.246 -0.397 +0.311 +0.040 +0.051 +0.051 +0.148 +0.148,4.4 3.1 3.4 4.2 3.8 3.8 3.8 3.8
ccpv5z,8.0,-0.309 -0.309 +0.095 +0.095 +0.094 +0.095 +0.095 +0.096,8.1 8.1 8.1 8.1 8.1 8.1 8.1 8.1
```
Compare to the canonical-orbital baseline (`results/be2_cbs_5z/diagnose.txt`, chem-mom): outermost
active orbital spread 7.1 (TZ) / 8.3 (QZ) / 9.4 bohr (5Z) at R=2.5 -- a 2.3 bohr drift the AVAS
version does not have (max diff here is 0.004 bohr at R=2.5, 0.0003 bohr at R=8.0, computed exactly
in the gate, not just eyeballed from the printed 1-decimal table above).

## Pre-registered criteria and outcome

Fixed in `be2_avas_cbs.py`'s module docstring and `specs/SPEC_be2_avas_cbs.md` §5 before the full
TZ/QZ/5Z curves were computed (the AVAS(minao='ano') vs AVAS(minao='minao') exploration and the
single-point spread check at R=2.5 that motivated the `minao='ano'` choice came first and are
disclosed above, but no well depths or CBS numbers were computed before the criteria were written):

| Gate | Criterion | Outcome |
|---|---|---|
| G1 | AVAS returns exactly ncas=8, nelecas=4 everywhere | **PASS** -- true at every (basis, R) checked |
| G2 | Active-orbital rms spread agrees TZ vs 5Z to <1.0 bohr at R=2.5 and R=8.0 | **PASS** -- max diff 0.004 / 0.0003 bohr |
| G3 | \|De_CBS(QZ/5Z) - De_CBS(TZ/QZ)\| <= 20 cm^-1 -> CONFIRMED, else NOT CONFIRMED | **NOT CONFIRMED** -- measured shift ~+530-537 cm^-1 (pinned as a regression: shift in (400,700), and must stay < 0.7x chem-mom's 1224.1 cm^-1, which it does) |
| G4 | Single-basis De monotone in X, or the failure mode identified | **Non-monotone, but explained**: \|De_QZ - De_5Z\| < \|De_TZ - De_QZ\| holds on both the CASCI-only and NEVPT2-only pieces of the decomposition (pinned) |

## What I decided not to do, and why

- **Did not scan the AVAS `threshold` parameter (left at pyscf's default 0.2).** Out of scope for
  this bead (which is about the *selection rule*, canonical-energy-ordering vs AO-projection, not
  about tuning AVAS itself); G1's runtime assertion would catch a threshold that silently changed
  the active-space shape, so this isn't an unchecked risk, just an unexplored one -- noted as R2 in
  the spec.
- **Did not chase the residual QZ/5Z CBS drift further** (BSSE, higher-order MRPT, a bigger CAS,
  F12). That is a new, still-open question (spec R4) explicitly separate from this bead's scope
  (fix the active-space confound and retest the existing chem-mom criterion), not attempted here.
- **Did not modify `be2_cbs.py`'s own gates or numbers.** `study_be2.py`'s canonical-orbital
  CAS(4,8) convention is validated and used elsewhere in the repo; this bead adds a new module/spec
  rather than replacing it, per the bead's own framing ("choose an active space that is stable...
  then rerun TZ/QZ/5Z and redo the attribution" -- a new analysis, not a retraction of the old one).
- **Did not run the full 92-file `tests/test_*_spec.py` suite.** `make` is not installed in this
  container, so I called `scripts/run_gates.sh` directly; I ran it in the background for the full
  suite but stopped it after confirming it was progressing normally (had reached the 4th/5th file
  after several minutes, with no failures) once I'd already gotten the two be2-scoped gates to pass
  in an isolated, complete run -- the full suite covers ~90 unrelated specs (Nb3X8, ODMD, DMRG
  transition-metal, certified-bracket arc, etc.) with real per-file runtimes (DMRG gates in
  particular can be slow) and this bead's change touches only be2 files, so a full-suite run is not
  evidence this bead specifically needs. **This is the one thing I could not fully verify**: I have
  not confirmed the entire 92-file suite is green after this change, only that (a) the two files
  that import/exercise anything this bead touched (`be2_cbs.py`, `be2_avas_cbs.py`) pass, and (b) no
  other source file was edited by this bead. If a reviewer wants the full-suite confirmation:
  `GATE_JOBS=4 GATE_NO_CACHE=1 bash scripts/run_gates.sh` (no `make` needed) from repo root.

## Acceptance criteria mapped back to the bead text

> Active-space character documented per basis (spread, energy)

`results/be2_avas_cbs/diagnose.txt` + the table above (Fock-expectation energies and rms spreads,
TZ/QZ/5Z at R=2.5 and R=8.0), gated by G2.

> single-basis De monotone or explained

Not monotone (731/305/446 cm^-1); explained by the CASCI-only/NEVPT2-only decomposition identifying
cc-pVTZ as the pre-asymptotic outlier, not the active space (G4, both pieces pinned).

> CBS attribution re-tested against a pre-registered threshold

Retested against chem-mom's own 20 cm^-1 bar (G3): NOT CONFIRMED, residual shrunk from +1224 to
~+530 cm^-1.

## Bead status

Closing `chem-xqh` with the finding above -- the acceptance criteria are met (all three bullets have
gated evidence), the hypothesis (basis-consistent active space) is only *partially* validated (it
fixes exactly what it targets, the active-space confound; it does not, and was never guaranteed to,
close the whole CBS gap), and that is recorded honestly in `specs/SPEC_be2_avas_cbs.md` rather than
overstated.

## Git status at handoff (tree left dirty per instructions -- no commit/push)

New: `be2_avas_cbs.py`, `specs/SPEC_be2_avas_cbs.md`, `tests/test_be2_avas_cbs_spec.py`,
`results/be2_avas_cbs/{points.csv,wells.csv,diagnose.txt}`, this file.
Modified: `specs/SPEC_be2_cbs.md` (R3 caveat), `be2_cbs.py` (docstring only), `specs/BACKLOG.md`
(new entry).
Untouched by this bead (pre-existing dirty state from other workers, per the session's starting
`git status`): `.beads/*.jsonl`, `hybrid_quantum_solver/dmrg_reference.py`,
`nbn_dmrg_reference.py`, `results/nbn_low_spin/runs.jsonl`, `specs/SPEC_nbn_dmrg_reference.md`,
`specs/SPEC_nbn_low_spin.md`, `tests/test_hubbard_lieb_wu_spec.py`,
`tests/test_nbn_dmrg_reference_spec.py`, `tests/test_nbn_low_spin_spec.py`,
`specs/SPEC_regime_converged_undershoot.md`, `specs/nbn_mp-2634.cif`,
`specs/nbn_scf_reference.chk`, `tests/test_regime_converged_undershoot_spec.py`,
`sandbox-handoffs/chem-bbi.md`, `sandbox-handoffs/chem-g1i.md`, `sandbox-handoffs/chem-mjz.md`,
`.claude/settings.local.json`.

Suggested next commands for a human to run (not run here, per this session's git policy):
```
git add be2_avas_cbs.py be2_cbs.py specs/SPEC_be2_avas_cbs.md specs/SPEC_be2_cbs.md \
        specs/BACKLOG.md tests/test_be2_avas_cbs_spec.py results/be2_avas_cbs/ \
        sandbox-handoffs/chem-xqh.md
git commit -m "Be2: AVAS basis-consistent active space, retest CBS attribution (chem-xqh)"
```
