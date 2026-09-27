# chem-owr handoff — scan-vs-chk discrepancy diagnosed: not a code bug, an SCF-instability coin flip

## Bead (verbatim, `bd show chem-owr`)

> While correcting SPEC_nbn_low_spin.md (chem-g1i), re-ran nbn_low_spin.py scan at d=3.7255 A --
> the real Nb-N separation the committed data/nbn_scf.chk actually encodes (per chem-bbi's fix in
> SPEC_nbn_dmrg_reference.md). Expected the scan pipeline's fresh SCF-from-synthetic-CIF
> reconstruction to reproduce the committed chk's energy at this matching distance. It does not:
>
>   E(10,4) scan  = -110.041189 Ha (exact FCI, synthetic single-cell CIF reconstruction)
>   E(10,4) chk   = -110.046028 Ha (exact FCI, real vendored data/nbn_scf.chk)
>   discrepancy   = 4.84 mHa
>
> Ruled out: the scan pipeline's own MIC/unwrapping bug ... so this is a different, undiagnosed
> mismatch between the two code paths (candidates not yet ruled out: basis-set/ECP assignment
> differences, SCF spin-scan differences, or a genuine remaining geometry difference).

No formal acceptance criteria were recorded on the bead (checked `bd show chem-owr --json`) — it
is an investigative bug report. I read it as: find and record the root cause, since three prior
sessions (chem-bbi, chem-g1i) had each ruled out one candidate and explicitly punted the rest here.

## Diagnosis

**Not a geometry, basis, ECP, or MIC bug.** Verified directly:

- The real vendored CIF's own raw (no-MIC) Nb-N distance is 3.7255149 Å; the scan's `write_cif`
  target was 3.7255 Å exactly — a 1.5e-5 Å difference, negligible at this energy scale.
- `_get_smart_basis`/`ecp` assignment is a pure function of atomic number (Nb=41 → ECP+basis,
  N=7 → all-electron) — identical on both code paths by construction, verified by reading
  `hybrid_quantum_solver/chemistry_gateway.py:_get_smart_basis` and `benchmark_nbn.py:124`.
- Running `benchmark_nbn.ground_state_mf` on the **real, vendored CIF directly**
  (`specs/nbn_mp-2634.cif`, not the synthetic reconstruction) still fails to reproduce the
  committed chk's energy — so the synthetic-CIF reconstruction path isn't even implicated. This is
  the key experiment: same file, same code, same distance, different answer.

**Root cause: the committed (10,4)/S=3 UHF solution sits essentially exactly at a `pyscf` internal
stability boundary, and `benchmark_nbn._tight_scf`'s Newton/SOSCF stability-following loop lands on
one of (at least) two genuinely, separately locally-stable UHF minima depending on multi-threaded
BLAS floating-point non-associativity — not on any geometry difference.**

- The two minima are 0.548 mHa apart in raw SCF energy and 4.84 mHa apart after
  CASCI(14,14)/FCI(10,4) (CASCI amplifies the gap — plausible if the differently-characterized
  solutions shuffle which physical orbitals fall inside the 14-orbital active window; not
  independently confirmed, not needed for the diagnosis).
- Both minima independently pass an explicit `mf.newton().stability(return_status=True)` check
  (`stat_i=True`, MOs unchanged) when checked directly — i.e. they are both real, stable
  stationary points, not one being a wrongly-flagged saddle.
- Which one `_tight_scf` finds is non-deterministic under default multi-threading (this container:
  8 cores) and deterministic under `OMP_NUM_THREADS=1`.

## Commands run and their exact output

**Baseline: real cif's own raw distance matches the target, ruling out geometry** (instant):
```
$ uv run python3 -c "
from ase.io import read as ase_read
import nbn_low_spin as nl
nl.write_cif('/tmp/scan_3p7255.cif', 3.7255)
a1 = ase_read('/tmp/scan_3p7255.cif'); a2 = ase_read('specs/nbn_mp-2634.cif')
print('scan dist (no mic):', a1.get_distance(0,1))
print('real dist (no mic):', a2.get_distance(0,1))
"
scan dist (no mic): 3.725500000009935
real dist (no mic): 3.725514948764577
```

**Decisive experiment 1 — fresh SCF on the REAL cif (not the reconstruction), default threading**
(~5s per run, this container, 8 cores, x86_64):
```
$ uv run python3 -c "
import benchmark_nbn as bn
bn.CHKFILE = '.dmrg_tmp/real_refit.chk'   # fresh, not restored
mf = bn.ground_state_mf('specs/nbn_mp-2634.cif')
print(mf.e_tot)
"
run1: -110.02813374559007   (found 'internal instability', followed it down)
run2: -110.02758538272639   (found 'stable' -- matches the committed chk to 9 decimals)
run3: -110.02813374568413   (instability, lower)
run4: -110.02813374500829   (instability, lower)
```
Reference: committed `specs/nbn_scf_reference.chk` e_tot = **-110.0275853827266** (matches run2).

CASCI(14,14)/FCI(10,4) on run1's orbitals: **-110.0411912172265** — matches the bead's reported
scan discrepancy (-110.041189) to 5 decimals. CASCI(14,14)/FCI(10,4) on run2's orbitals:
**-110.04602841723879** — matches the committed reference (-110.046028) to 6 decimals. Same code,
same file, same distance — the *only* thing that differs between runs is which SCF solution
`_tight_scf`'s stability loop lands on.

**Decisive experiment 2 — same command, `OMP_NUM_THREADS=1`** (removes BLAS-threading
non-associativity):
```
$ OMP_NUM_THREADS=1 uv run python3 -c "... same as above ..."
run5: -110.02813374576796
run6: -110.02813374576796   (bit-identical to run5)
```
Both single-threaded trials deterministically find the instability and land on the lower solution
— never on the committed one.

**Confirming both solutions are genuinely, separately stable** (not a false instability report on
one of them):
```
$ uv run python3 -c "
from pyscf import lib, scf
mol = lib.chkfile.load_mol('specs/nbn_scf_reference.chk')
res = lib.chkfile.load('specs/nbn_scf_reference.chk', 'scf')
mf = scf.UHF(mol); mf.mo_coeff, mf.mo_occ = res['mo_coeff'], res['mo_occ']
mfn = mf.newton(); mfn.mo_coeff, mfn.mo_occ = mf.mo_coeff, mf.mo_occ; mfn.conv_tol = 1e-10
mo_new, _, stat_i, _ = mfn.stability(return_status=True)
print(stat_i, allclose(mo_new, mfn.mo_coeff))
"
-> wavefunction is stable in the internal stability analysis
   True True
```
(the lower solution's own stability check, printed inline during the `_tight_scf` loop above, also
says "stable" once it's reached — both are real minima.)

Wall time: every run above is a single UHF/CASCI(14,14) construction, 5-9s each; none crossed the
1-minute reporting threshold.

## What changed

- **`specs/SPEC_nbn_low_spin.md`** — §0 (the "second anomaly" paragraph), §7 (Out of scope), §8
  (Caveats), and the Table 3 caption updated from "undiagnosed" to the diagnosis above, with the
  exact reproducing commands and numbers. No table values changed (nothing above was wrong given
  the chk it used — the *chk's own status* is now in question, filed separately). No code files
  touched.
- **No `.py` files changed.** This bead's answer did not require a code fix — the "bug" was never
  in `nbn_low_spin.py`'s reconstruction; ruling that out *was* the diagnosis.
- **Filed `chem-y6w`** (P2): the committed `data/nbn_scf.chk` / `specs/nbn_scf_reference.chk` is
  itself the *higher*-energy, non-deterministically-reached one of two locally-stable solutions —
  a deterministic (single-threaded) rerun of the current code reliably finds a lower one instead.
  Every number built on the committed chk (chem-bbi's G3/G3b, chem-g1i's Tables 1-2, the S1-S3 gap,
  G4/G5) is anchored to it; none of them are computed wrong given that chk, but whether that chk
  should keep being "the" reference is now an open, filed question — explicitly out of chem-owr's
  scope (which only asked why the scan pipeline disagreed with it).

## Pre-registered criteria vs outcome

The bead had no pre-registered acceptance criteria (checked `bd show --json`; description only).
I treated "diagnose the discrepancy, ruling in/out the candidates the bead itself listed" as the
implicit ask, since that is exactly what the three prior sessions' handoffs (chem-bbi, chem-g1i)
each explicitly deferred to this bead. All three listed candidates were addressed:

| Candidate (bead's own list) | Outcome |
|---|---|
| Basis-set/ECP assignment differences | Ruled out — pure function of atomic number, identical on both paths by inspection of the code |
| SCF spin-scan differences | Ruled out — both paths pick ground spin S=3 (2S+1=7); the divergence happens *within* the S=3 sector's own `_tight_scf`, not at spin selection |
| Genuine remaining geometry difference (lattice constants, MIC) | Ruled out — real cif's own raw distance (3.7255149 Å) matches the reconstruction's target (3.7255000 Å) to 1.5e-5 Å; using the real cif directly (bypassing the reconstruction entirely) still reproduces the discrepancy |
| (not on the bead's list, found here) SCF-solution multiplicity + BLAS-threading non-determinism in the stability-following loop | **Confirmed as the actual root cause** — reproduced 6 times (4 default-thread, 2 single-thread), with an independent stability check on both solutions |

## What I decided not to do, and why

- **Did not re-vendor a "corrected" chk or touch `nbn_dmrg_reference.py` / `benchmark_nbn.py`.**
  Fixing the non-determinism (e.g. pinning `OMP_NUM_THREADS`, or deciding which of the two minima
  is "the" reference) ripples into two already-closed beads' worth of gated numbers
  (`test_nbn_dmrg_reference_spec.py` G3/G3b, `test_nbn_low_spin_spec.py` G4/G5) and is a real
  decision for a human or a dedicated bead, not something to fold into a P3 diagnosis bead. Filed
  as **chem-y6w** (P2, reflects that it calls into question already-published reference numbers)
  instead.
- **Did not pin down *why* CASCI(14,14) amplifies a 0.548 mHa SCF gap to 4.84 mHa** (candidate:
  the two solutions' differing orbital character shuffles which physical orbitals fall inside the
  active window, versus the same window with a rotated interior). Not needed to answer "is this a
  code bug" (it isn't); flagged as unexplored, not fixed here.
- **Did not investigate whether pyscf version, BLAS backend (OpenBLAS vs MKL), or a specific
  thread count threshold controls the flip rate.** 8-thread default vs 1-thread was sufficient to
  demonstrate the mechanism; a full characterization (2/4/8/16 threads, repeat counts) would be
  more precision than this bead needs.
- **Did not run the ~7-8 minute DMRG spec gates** (`test_nbn_dmrg_reference_spec.py`,
  `test_nbn_low_spin_spec.py`). No `.py` file was changed this session (only the `.md` spec), and
  neither test file references the changed prose — `grep -n "chem-owr" tests/test_nbn_low_spin_spec.py`
  returned nothing. Re-running them would not exercise anything different from what chem-bbi/
  chem-g1i already verified green.
- **Did not commit, push, or `bd dolt push`** — per standing session instructions.

## What I could not verify

- Whether the same non-determinism reproduces on a different host/BLAS backend than this
  container's (x86_64, 8 cores, whatever OpenBLAS/MKL build `uv sync` pulled in) — only tested
  here.
- The exact numerical size of the near-zero Hessian eigenvalue driving the flip (would need to
  instrument pyscf's stability solver internals; inferred its existence from the flip behavior,
  did not measure it directly).

## Gate / lint status

```
$ uv run --extra dev ruff check specs/ nbn_low_spin.py benchmark_nbn.py nbn_dmrg_reference.py
All checks passed!
```
No test files were touched; `tests/test_nbn_low_spin_spec.py` and
`tests/test_nbn_dmrg_reference_spec.py` were last verified green by chem-g1i/chem-bbi (see their
handoffs) and are unaffected by this session's spec-only edit.

## Git status — NOT committed, NOT pushed (per standing instructions)

```
 M specs/SPEC_nbn_low_spin.md     <- chem-owr (this session: §0/§7/§8 + Table 3 caption diagnosis)
 M .beads/interactions.jsonl      <- bd bookkeeping (claim/notes on chem-owr, creation of chem-y6w)
 M .beads/issues.jsonl            <- same
```
Plus pre-existing dirty files from other sessions (chem-bbi/chem-g1i/chem-mjz and others) — see
`sandbox-handoffs/chem-bbi.md` and `sandbox-handoffs/chem-g1i.md` for those; none touched here.

**Suggested commands for a human to run** (chem-owr-scoped file only):
```bash
git add specs/SPEC_nbn_low_spin.md
git commit -m "nbn: diagnose scan-vs-chk discrepancy as SCF-instability threading nondeterminism (chem-owr)

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```
The `.beads/*.jsonl` files pick up the chem-owr claim/close and the chem-y6w creation; stage
alongside if the sync convention here expects it.

## Bead status

Closing **chem-owr** — the discrepancy is diagnosed with reproduced, decisive evidence: it is not
a code bug in `nbn_low_spin.py`'s geometry/basis reconstruction (all three candidates the bead
itself listed are ruled out), but genuine SCF-solution multiplicity in the committed (10,4)/S=3
sector combined with BLAS-threading non-determinism in `pyscf`'s stability-following loop. The
consequence this raises for the committed reference's own status is filed as **chem-y6w**, not
resolved here — appropriately out of scope for a bead whose question was "why do these two numbers
disagree," now answered.
