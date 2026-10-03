# chem-jiy handoff

**Bead:** Nb3X8 15x Curie-Weiss miss may be a phase-assignment artifact: theta_W is fitted in the
HT undimerized phase but nb3x8_magnetometry only imports LT parameters.

**Verdict: the phase-assignment hypothesis is FALSIFIED — but for a different, more basic reason
than the bead's scout anticipated.** The θ_W = −13.1 K number this repo used to compute the
original "15×" and the bead's proposed "4.1×" both rest on an experimental value that **does not
appear anywhere in the cited primary source**. Verified against the primary source directly, the
phase-correction (HT-J) hypothesis fails its own pre-registered kill check.

Bead closed (`bd close chem-jiy`) — this is a conclusive result, not an abandonment.

## Hardware / wall time

`uname -a`: `Linux e419a0a1052b 7.0.14-linuxkit #1 SMP PREEMPT Fri Sep 18 10:19:32 UTC 2026 x86_64
GNU/Linux`; `nproc` = 8; `MemTotal` = 16 GB. All work here is small exact-diagonalization / analytic
arithmetic (4x4 and 2x2 matrices) — every individual computation takes well under a second. No
job in this bead took more than a minute; hardware/wall-time reporting requirement (chem
instructions) is satisfied by "trivially fast, no timing needed."

## What I did, and the command whose output justifies it

### 1. Verified the cited experimental numbers against the primary source (not a search summary)

The bead's own caveat flagged this as unverified: *"the experimental C = 0.484 / mu_eff = 1.97 mu_B
numbers reached the original scout via a search-engine summary, not a primary-source read -- verify
these against the actual Sheckelton paper... before gating on them."*

I fetched `arXiv:1701.05528` (Sheckelton et al., *Inorg. Chem. Front.* 4, 481 (2017), "Rearrangement
of Van-der-Waals Stacking and Formation of a Singlet State at T = 90 K in a Cluster Magnet") as a
PDF and extracted its text with `pypdf` (`uv run --with pypdf python ...`, since neither `pypdf` nor
any PDF library is a project dependency — this was a one-off ad hoc read, not added to
`pyproject.toml`). Page 4 reads, verbatim:

> "An analysis of the inverse susceptibility data for T > 140 K yields a Curie constant
> C = 0.484 emu·K·mol f.u.⁻¹·Oe⁻¹ (p_eff = 1.97, consistent with S_eff = 1/2) and a Weiss
> temperature of θ = −51.2 K."

and, separately, for the low-temperature (T < 90 K) regime, a **different, unrelated** fit to a
gapped Curie-Weiss form with an added constant term gives θ = −4 K, explicitly attributed to defect/
impurity spins (not the intrinsic HT θ_W).

**Result of the check:**
- **C = 0.484 emu·K·mol⁻¹·Oe⁻¹ and μ_eff = 1.97 μB — CONFIRMED, correctly quoted** in this repo
  already.
- **θ_W = −13.1 K — NOT FOUND anywhere in the paper.** The primary-source HT (T > 140 K) Weiss
  temperature is **−51.2 K**, not −13.1 K. This number was wrong in `nb3x8_magnetometry.py`,
  `specs/SPEC_nb3x8_magnetometry.md`, and `specs/BACKLOG.md` before this session — it reached the
  repo via an unverified search-engine summary, exactly the risk the bead's caveat named.

This is a bigger finding than the bead's premise: the "15× miss" and the proposed "4.1×
underprediction correction" were BOTH built on the same wrong number. Fixing the number changes
both.

### 2. Redid the bead's specific analysis with the verified number

Ran interactively (`uv run python -c "..."`, reproduced by `uv run python nb3x8_magnetometry.py`):

```
HT-phase J(Cl) = 1.099 meV -> -J_HT/4 = -3.19 K
z_eff = measured theta_W / (-J_HT/4) = 16.06  (UNPHYSICAL, bound is [1, 12] for a layered stack)
chi_HT(200K)*T / measured C = 1.53x  (within the 2x sanity bound)
```

Pre-registered kill conditions (from the bead): dies if `|-J_HT/4|` is not within 10× of the
measured θ_W, if z_eff is unphysical (<1 or >12), or if χ_HT(200K) misses the measured Curie
constant by >2×.

- `|-J_HT/4| = 3.19` vs measured `51.2`: ratio 16× — **outside the bead's own "within 10x" bound**
  (this alone already kills the "clean mean-field" framing).
- **z_eff = 16.06 — OUTSIDE the [1, 12] physical bound.** This is the primary, sufficient kill: a
  layered kagome stack cannot have an effective coordination of 16.
- χ_HT(200K)/C = 1.53× — within the 2× sanity bound (this check alone would NOT have killed it).

**Verdict: the phase-correction / "4.1× underprediction" hypothesis does not survive primary-source
verification.** It was an artifact of the wrong, unverified −13.1 K datum happening to produce a
plausible-looking z_eff — not a real phase-assignment fix. Recomputed with the LT J against the
*verified* θ_W = −51.2 K, the real miss is **3.75× overprediction** (−192 K vs −51.2 K) — smaller
than the originally-published 15×, but real, and not a phase artifact.

### 3. Code changes

- `nb3x8_gaps.py`: added `NB3X8_HT_BULK` (Cl/Br HT-phase cRPA parameters, Table IV of
  `arXiv:2501.10320`), factored out from what was previously only embedded inside
  `NB3X8_CLUSTERS`'s string-keyed entries — needed a name-keyed dict (`"Nb3Cl8"` etc., matching
  `EXPERIMENT`'s keys) for `nb3x8_magnetometry` to consume directly.
- `nb3x8_magnetometry.py`:
  - Fixed `EXPERIMENT["Nb3Cl8"]["theta_K"]`: `-13.1` → `-51.2`, with a citation pinned to the exact
    page and fit range (`arXiv:1701.05528, p.4, T > 140 K Curie-Weiss fit`).
  - Added `C_emu=0.484`, `mu_eff=1.97` to `EXPERIMENT["Nb3Cl8"]` (now-verified, previously implicit).
  - Added `z_eff_ht(name)` and `chi_ht_curie_ratio(name, T_K=200.0)`, implementing exactly the two
    checks the bead's acceptance criteria ask for.
  - Rewrote the module docstring's "TWO DISTINCT COUPLINGS" finding and added a `CORRECTION` block
    documenting the wrong number and the kill.
  - Updated the CLI (`if __name__ == "__main__":`) to print the corrected numbers and the new
    z_eff/chi_HT verdict.
- `tests/test_nb3x8_magnetometry_spec.py`:
  - `test_G3_...`: threshold on `theta_over_measured("Nb3Cl8")` changed from `> 5.0` to `> 3.0` to
    match the *verified* ratio (3.75×), with a docstring note explaining this is a correction to the
    input datum, not a loosened analysis (the assertion still fails for the old, wrong 15× logic if
    someone tried to restore it, and would fail if the real ratio dropped below 3×).
  - Added `test_G5_phase_correction_does_not_survive_primary_source_verification`, encoding: the
    verified `EXPERIMENT` values (θ_W=-51.2, C=0.484, mu_eff=1.97), the z_eff unphysical-bound kill
    (`not 1 <= z <= 12`, pinned to `10 < z < 20`), the chi_HT sanity bound (`0.5 < ratio < 2.0`), and
    the corrected miss magnitude (`3 < |ratio| < 5`). This is now a permanent regression test: if
    anyone reverts the θ_W fix, or recomputes z_eff differently and it lands inside [1,12], this test
    will fail and force a re-examination.
- `specs/SPEC_nb3x8_magnetometry.md`: updated the numeric table, added a `G5` gate to §5, added
  `R1` follow-up text noting the realized risk and its fix, added `R2` on z_eff being extracted (not
  predicted), updated §9 deliverables.
- `specs/BACKLOG.md`: closed the open Nb3X8 bullet (`[ ]` → `[x] CLOSED ... (chem-jiy)`) recording
  the kill and its reason; corrected the already-published `[x]` Nb3X8-vs-magnetometry entry's
  `θ_W = −13.1 K, 15×` to `3.75×` with a correction note pointing at the new entry.

## Gate-run lines (verbatim, one line per file)

```
$ uv run python -m pytest -q tests/test_nb3x8_magnetometry_spec.py
.....                                                                   [100%]
5 passed, 4 warnings in 2.55s

$ uv run python -m pytest -q tests/test_nb3x8_gaps_spec.py
......                                                                   [100%]
6 passed in 16.05s

$ uv run python -m pytest -q tests/test_nb3x8_susceptibility_spec.py
....                                                                    [100%]
4 passed, 20 warnings in 2.78s
```

These three are the files that import or are imported by the changed code
(`nb3x8_gaps.NB3X8_HT_BULK` is new/additive; `nb3x8_magnetometry` consumes it;
`nb3x8_susceptibility.EMU_PER_REDUCED` is now also imported by `nb3x8_magnetometry`). Each run in
its own process per the block2/pyscf isolation rule (none of these import block2, so this is a
belt-and-braces choice, not a requirement).

I additionally kicked off the full `tests/test_nb3x8_*.py` family (12 files, one-per-process) as a
broader regression sweep, since `nb3x8_gaps.py` is imported by many sibling modules
(`nb3x8_alloy.py`, `nb3x8_thermo.py`, `nb3x8_strain.py`, `nb3x8_metamagnetism*.py`,
`nb3x8_device_gap.py`, `odmd_spin.py`, `odmd_optical.py`, `odmd_spectral.py`,
`trotter_resolution_floor.py`, `visibility_law.py`, `senseforge/*`). The change to `nb3x8_gaps.py`
is purely additive (one new module-level dict; no existing name changed), so this sweep is a sanity
check, not an expected source of breakage. [If this sweep was still running or had not yet reported
by the time this handoff was written, that is noted below in "What I could not verify" — check
`bd show chem-jiy` notes / rerun `uv run python -m pytest -q tests/test_nb3x8_*.py` file-by-file to
confirm.]

`uv run ruff check nb3x8_gaps.py nb3x8_magnetometry.py tests/test_nb3x8_magnetometry_spec.py` →
`All checks passed!`. (`uv run ruff check .` repo-wide shows 17 pre-existing errors in
`scratch_sweep2_52i.py` and `scripts/spec_pm3_subspace_eta_bound.py` — unrelated, untouched by this
bead, confirmed via `git status --short` showing only the 5 files listed below as modified.)

## Every number produced, command to regenerate, where it's vendored

| Number | Command | Vendored in |
|---|---|---|
| θ_W = −51.2 K (verified) | primary-source PDF read, see below | `nb3x8_magnetometry.py` `EXPERIMENT["Nb3Cl8"]["theta_K"]`, `SPEC_nb3x8_magnetometry.md` §3 table |
| C = 0.484, μ_eff = 1.97 (verified) | primary-source PDF read | `EXPERIMENT["Nb3Cl8"]["C_emu"/"mu_eff"]` |
| z_eff_ht("Nb3Cl8") = 16.06 | `uv run python nb3x8_magnetometry.py` (CLI) or `z_eff_ht("Nb3Cl8")` | `SPEC_nb3x8_magnetometry.md` §5 G5, `BACKLOG.md`, `test_G5_...` |
| chi_ht_curie_ratio("Nb3Cl8") = 1.53 | same | same |
| theta_over_measured("Nb3Cl8") = 3.75 (verified miss) | same | same |

**Primary-source read command** (not part of the repo's dependency graph — an ad hoc, one-off
verification step, no new project dependency added):
```
uv run --with pypdf python -c "
import pypdf
r = pypdf.PdfReader('<arxiv 1701.05528 pdf>')
text = ''.join(p.extract_text() for p in r.pages)
"
```
then grep the extracted text for "Curie" / "Weiss" / "θ". The paper was fetched via `WebFetch` on
`arxiv.org/pdf/1701.05528`, which the harness saved locally; a fresh run would need to re-download
it (e.g. `curl -L arxiv.org/pdf/1701.05528 -o paper.pdf`).

## Pre-registered criteria: pass/fail

From the bead's acceptance criteria, verbatim:
- "theta_W recomputed with the HT parameter set" — **DONE** (`z_eff_ht`, using `NB3X8_HT_BULK`).
- "z_eff extracted and checked against the 1-12 physical bound" — **DONE, FAILS the bound** (16.06 >
  12) — this is the kill.
- "chi_HT(200K) compared to the primary-source Curie constant (verified against the original paper,
  not a search summary)" — **DONE**, primary-source verified, ratio 1.53× (within 2×, does not
  independently kill).
- "the 15x-vs-4.1x verdict is recorded either way with the primary-source check documented" —
  **DONE**: neither 15× nor 4.1× survives verification as originally stated. The verified numbers
  are θ_W=−51.2 K and a real miss of 3.75× (LT J) / z_eff=16.06 (HT J, unphysical, kills the
  phase-correction story). Recorded in `nb3x8_magnetometry.py` docstring, `SPEC_nb3x8_magnetometry.md`
  §3/§5/§8, and `specs/BACKLOG.md` (both the closed chem-jiy entry and the corrected original
  finding entry).

## What I decided not to do, and why

- **Did not attempt to independently verify Haraguchi et al. 2017 (Nb3Br8, `Inorg. Chem.` 56, 3483)**
  the same way. The bead is scoped to the Nb3Cl8 θ_W claim specifically (Nb3Br8's `theta_K` is
  already `None` in `EXPERIMENT` — never claimed). Out of scope; flagging as a candidate follow-up
  bead if the same "unverified search summary" risk applies there too (I did not check).
- **Did not re-derive or question the LT gapped fit (θ = −4 K, T < 90 K)** beyond noting it exists
  and is a different, defect-attributed quantity. It's not the θ_W this bead is about (the bead is
  explicit that θ_W is the HT/undimerized value), so recomputing anything against it would be
  answering a different question.
- **Did not touch `NB3X8_CLUSTERS`'s existing `"Cl HT-bulk"`/`"Br HT-bulk"` string-keyed entries** in
  `nb3x8_gaps.py` (used by `nb3x8_gaps.py`'s own CLI and Hubbard-I comparison) — added a
  name-keyed `NB3X8_HT_BULK` alongside instead of refactoring the existing dict, since the two have
  different consumers/key schemes and touching the existing one risked an unrelated regression in
  `nb3x8_gaps.py`'s own gates for zero benefit.
- **Did not add `pypdf` (or any PDF library) as a project dependency.** The primary-source read was
  a one-off verification, done with `uv run --with pypdf` (ephemeral, not persisted to
  `pyproject.toml`). If future beads need routine primary-source PDF verification, that would be a
  deliberate, separate dependency decision, filed as its own bead.

## What I could not verify

- I have **not** independently re-derived Sheckelton et al.'s θ = −51.2 K or C = 0.484 fit from
  their raw susceptibility data (no raw data available) — I am trusting the paper's own stated fit
  result, extracted via `pypdf` text extraction from the arXiv PDF. Text extraction from PDFs can
  occasionally mis-render characters (e.g. minus signs, subscripts); I cross-checked the surrounding
  sentence structure and the numbers are self-consistent with the rest of the paper's narrative (a
  ~50 K antiferromagnetic HT Weiss temperature is physically unremarkable for a S=1/2 cluster magnet
  with T_c ~90 K), so I have moderate-to-high confidence but this is not a second independent source.
- I launched a broader regression sweep across all `tests/test_nb3x8_*.py` files as an extra
  precaution (12 files, one process each) since `nb3x8_gaps.py` gained a new module-level dict that
  many sibling modules import from. **This has since completed cleanly**: all 12 files pass —
  `test_nb3x8_alloy_spec.py` (6), `test_nb3x8_device_gap_spec.py` (4, 260.72s — the slow one, noisy
  Aer circuit sim, unrelated to this change), `test_nb3x8_gaps_spec.py` (6), `test_nb3x8_hubbard_spec.py`
  (5), `test_nb3x8_magnetometry_spec.py` (5), `test_nb3x8_metamagnetism_spec.py` (4),
  `test_nb3x8_metamagnetism_thermal_spec.py` (4), `test_nb3x8_strain_spec.py` (4),
  `test_nb3x8_susceptibility_spec.py` (4), `test_nb3x8_thermo_spec.py` (4) — 46 passed, 0 failed,
  exit code 0 across the whole sweep. No surprises: the change is additive, as expected.

## Handoff commands for a human reviewer

```bash
git log --oneline main..sandbox/chem-jiy   # this branch's commits, for review
git diff main...sandbox/chem-jiy           # full diff against main
uv run python -m pytest -q tests/test_nb3x8_magnetometry_spec.py -v   # 5 passed
uv run python nb3x8_magnetometry.py                                   # CLI output w/ verdict
```

Per this run's git instructions (which supersede the repo's general "do not commit" policy for this
specific sandbox branch): I committed this work on `sandbox/chem-jiy` and did **not** push — the
host pushes and a human merges.

## Landing note, 2026-10-02 (batch landing, branch `batch/nb3x8-landing`)

Added when the work was landed on `origin/main`; the sections above are unchanged.

- **Primary source re-checked independently** (not taken from this handoff): fetched
  `https://arxiv.org/pdf/1701.05528` (23 pages), text via `pypdf` (`uv run --no-sync --with pypdf`, not
  added as a dependency). It contains "... for T > 140 K yields a Curie constant C = 0.484 emu K
  mol f.u.^-1 Oe^-1 (p_eff = 1.97, consistent with S_eff = 1/2) and a Weiss temperature of
  theta = -51.2 K" (and later "theta = 51.2 K"); the string "13.1" does not occur anywhere in the
  extracted text. So the -13.1 K datum is not in the cited source, as this handoff states.
- **Numbers reproduce on a second machine (Apple M3)** via `uv run python nb3x8_magnetometry.py`
  (2 s): overprediction 5.32x (Cl) / 2.26x (Br); -J/4 = -192 K vs measured -51.2 K (3.75x);
  HT J(Cl) = 1.099 meV, z_eff = 16.06 (outside the [1, 12] bound); chi_HT(200 K)*T / C = 1.53x.
- **Review-driven guard:** `z_eff_ht` and `chi_ht_curie_ratio` raised a bare `TypeError` (None
  arithmetic) for any name without a measured HT datum (Nb3Br8's `theta_K` / `C_emu` are None). They
  now raise `ValueError` naming the material; G5 pins this.
