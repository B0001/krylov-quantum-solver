# chem-cfl — Nb3X8 Tc overprediction: candidate 1 (cooperative/first-order structural transition)

## Outcome

Ran the pre-registered `F_phase(T)` free-energy crossing check
(`nb3x8_tc_attribution.UNTESTED_CANDIDATES["cooperative_structural_transition"]`, pre-registered by
chem-phh) using elastic/lattice data sourced and primary-text-verified from a second paper this
sandbox could actually fetch over the network: Mortazavi, Zhuang & Rabczuk, "A first-principles
study on the physical properties of two-dimensional Nb3Cl8, Nb3Br8 and Nb3I8," `arXiv:2412.03360v1`.

**Verdict: candidate 1 SURVIVES** (both halides) — see "What was decided not to claim" below for
the honest limits of that verdict. Closed `chem-cfl` with `bd close`.

**Important secondary finding, not the bead's main deliverable but worth flagging:** chem-phh and
chem-a04 both concluded "no working internet access" in this sandbox — but that conclusion was based
only on `WebSearch`, which still returns a canned refusal every time I tried it. `WebFetch` on
direct URLs (arXiv abstract pages, the `export.arxiv.org/api/query` API, and PDF URLs — whose binary
content it saves locally for the `Read` tool to then extract as full multi-page text) **does** work.
This is how `arXiv:2412.03360` was sourced. I did not act on this for chem-b0y (in-plane kagome
exchange) — that's explicitly out of this bead's scope — but I noted it in both specs and
`UNTESTED_CANDIDATES["in_plane_kagome_exchange"]["blocked_on"]` so a future session doesn't repeat
the same dead end.

## What I changed

1. **Brought forward from the unmerged `sandbox/chem-phh` branch** (never merged to `main`, but
   `bd show chem-phh` confirms the bead itself is closed — beads state syncs independently of git
   branches): `nb3x8_tc_attribution.py`, `specs/SPEC_nb3x8_tc_attribution.md`,
   `tests/test_nb3x8_tc_attribution_spec.py`, and chem-phh's edits to `nb3x8_gaps.py`,
   `specs/SPEC_nb3x8_magnetometry.md`, `specs/BACKLOG.md` — via
   `git show sandbox/chem-phh:<path> > <path>`, without switching branches (per this session's git
   policy). Verified the brought-forward baseline reproduces cleanly before touching anything
   (`22 passed`, see below).
2. **`nb3x8_tc_attribution.py`** — added, reusing only already-validated primitives:
   - `mean_energy(U0, t, Us, T)` — exposes `<E>(T)` from the same N=2 Boltzmann trace
     `nb3x8_thermo.entropy`/`heat_capacity` already run internally (no new physics; `entropy`
     computes `e_avg` internally and discards it, this just returns it).
   - `NB3X8_EXFOLIATION` — vendored exfoliation energies (γ, J/m²) and monolayer hexagonal lattice
     constants (a, Å) from arXiv:2412.03360, with full citation.
   - `exfoliation_energy_per_fu(name)` — converts γ to meV per Nb3X8 formula unit
     (area = (√3/2)·a², assuming 1 formula unit per hexagonal monolayer cell).
   - `required_e_lattice(name)` — solves `E_lattice = Tc·ΔS_mag(Tc) − ΔE_mag(Tc)` directly (linear
     in `E_lattice`, no root-finding) using `NB3X8_HT_BULK`/`NB3X8_LT_BULK`.
   - `candidate1_verdict(name)` — the kill/survive dict, mirroring `candidate3_verdict`'s pattern.
   - Removed `"cooperative_structural_transition"` from `UNTESTED_CANDIDATES` (it is tested now);
     updated `"in_plane_kagome_exchange"`'s `blocked_on` text with the WebFetch finding.
   - Updated the module docstring and `__main__` CLI block to report both candidates.
3. **`tests/test_nb3x8_tc_attribution_spec.py`** — added gates G5–G7 (below); renumbered the old G5
   (honest-boundary gate) to G4 and narrowed it to candidate 2 only, since candidate 1 is no longer
   in `UNTESTED_CANDIDATES`.
4. **`specs/SPEC_nb3x8_tc_attribution.md`** — added §3b (candidate 1 approach, data, result);
   updated the public interface, gates G4–G7, out-of-scope, caveats (new R4), and deliverables
   sections; updated the status line.
5. **`specs/BACKLOG.md`** — added a new `[x] CLOSED 2026-10-15 (chem-cfl)` entry recording the
   sourced reference, the check, the exact numbers, and the caveat, directly under chem-phh's entry.
6. **`specs/SPEC_nb3x8_magnetometry.md` §7** — added a short cross-reference update (chem-phh had
   already left a pointer here for consistency; I extended it) noting candidate 1 is now tested.
7. **`nb3x8_gaps.py`** docstring — one-line update to the existing chem-phh note, same pattern.

No changes to any validated physics primitive (`nb3x8_susceptibility.py`, `nb3x8_thermo.py`,
`nb3x8_gaps.py`'s tables/functions, `nb3x8_magnetometry.py`) beyond the docstring note in (7).

## The sourced data — verify against primary text (chem-jiy lesson)

`WebFetch` on `https://arxiv.org/pdf/2412.03360` saved the PDF locally; the `Read` tool then
extracted the full multi-page text (not an abstract or summary). Verbatim numbers used, read
directly from that text:

- Exfoliation energy: Nb3Cl8 = 0.24 J/m², Nb3Br8 = 0.27 J/m² (Nb3I8 = 0.28 J/m², not used).
- Monolayer hexagonal lattice constant `a`: Nb3Cl8 = 6.783 Å, Nb3Br8 = 7.121 Å.

These are independent of `arXiv:2501.10320` (the only source `nb3x8_gaps.py` otherwise vendors):
this is a DFT structural/elastic study (VASP, PBE-GGA+D3) with no electronic or magnetic parameters,
and it reports only the single known (LT, dimerized) ground-state structure — it does **not**
resolve an HT phase. That is the central limitation of using it (see R4 below).

## The check and the exact numbers

Formula (derived, not fit): `F_phase(T) = E_mag,phase(T) - T*S_mag,phase(T) + E_lattice,phase`,
LT carrying `-E_lattice` relative to HT. Setting `F_HT(Tc) = F_LT(Tc)` at the MEASURED Tc and solving
for `E_lattice` (linear, exact):

```
E_lattice = Tc_measured * [S_mag,HT(Tc) - S_mag,LT(Tc)] - [E_mag,HT(Tc) - E_mag,LT(Tc)]
```

Reproduced by `uv run python nb3x8_tc_attribution.py`:

```
Candidate 1 (cooperative/first-order structural transition, F_phase(T) crossing):
  Nb3Cl8: needs E_lattice = 9.94 meV to match measured Tc (1.7% of the 596.9 meV exfoliation-energy bound, arXiv:2412.03360) -> survives
  Nb3Br8: needs E_lattice = 41.61 meV to match measured Tc (5.6% of the 740.1 meV exfoliation-energy bound, arXiv:2412.03360) -> survives
```

Both `E_lattice` values are positive (correct sign: LT must be net lattice-favored at T=0, matching
the observed LT ground state) and small fractions of the exfoliation-energy upper bound — comfortably
inside, not tuned right up against it (robust even to a 2× tightening of the bound; see R4(b)).

Checked the "unique crossing" assumption numerically (`test_G7_*`): scanning `E_lattice_needed(T)`
over `[0.01, 3000]` meV, it is non-decreasing everywhere for both halides — so a fixed `E_lattice`
cannot cross `F_HT − F_LT` more than once, and the crossing temperature is unambiguous.

## What was decided not to claim (honest limits, stated up front — not rounded into the verdict)

1. **The exfoliation-energy bound is a generous upper bound, not a measured HT/LT energy
   difference.** `arXiv:2412.03360` does not resolve an HT phase — no paper directly measures the
   quantity this check actually needs (the energy cost of the specific HT→LT distortion). Full
   layer separation is a defensible *ceiling* on any partial structural relaxation, nothing more.
2. **The formula-unit-per-cell assumption is a convention, not a re-verified atom count.** I used
   the standard 1-formula-unit hexagonal monolayer cell. Even halving the bound (2 formula
   units/cell: 298/370 meV) leaves it far above the needed `E_lattice` (9.9/41.6 meV) — the verdict
   is robust to this specific ambiguity, but it remains an assumption I did not independently
   re-derive from the paper's POSCAR atom counts in this pass.
3. **`NB3X8_HT_BULK`/`NB3X8_LT_BULK` are single-dimer electronic parameters with no lattice
   feedback.** The check asks "is a lattice energy of this size plausible," not "does this specific
   model predict Tc=90 K/382 K ab initio." Surviving the check rules out only energetic
   implausibility — it is not confirmation that the cooperative structural transition is the
   correct mechanism.
4. Candidate 2 (in-plane kagome exchange) is **still untested** — neither vendored source reports
   an in-plane exchange constant. Left pre-registered in `UNTESTED_CANDIDATES`, now with the
   WebFetch note instead of the stale WebSearch-only "no internet" claim.

## Commands run and their exact output

```
$ uv sync --extra dmrg --extra test --extra dev        # environment setup, done first

$ uv run pytest tests/test_nb3x8_tc_attribution_spec.py -v
... 14 passed, 22 warnings in 1.54s

$ uv run pytest tests/test_nb3x8_gaps_spec.py tests/test_nb3x8_magnetometry_spec.py tests/test_nb3x8_tc_attribution_spec.py -q
... 28 passed, 23 warnings in 20.10s

$ uv run ruff check .
All checks passed!

$ uv run python nb3x8_tc_attribution.py
(full CLI report above; candidate 3 KILLED for both halides as before, candidate 1 survives for
both halides, candidate 2 reported BLOCKED)
```

**What I did NOT run:** the full `make gates` (115+ spec files, each in its own process) or
`make test`. I ran only the three directly-relevant gate files (this module's own, plus its two
upstream dependencies, `test_nb3x8_gaps_spec.py` and `test_nb3x8_magnetometry_spec.py`, which pin
the tables and predictor this module reuses) plus a repo-wide `ruff check .`. I did not touch any
other module, so a regression elsewhere is very unlikely, and chem-phh made the same resource-
conscious call for the identical reason (shared 8-CPU/16GB compute envelope). If a reviewer wants
the full gate sweep run, `make gates` is the command.

## Acceptance criteria (verbatim from `bd show chem-cfl`) — pass/fail

- "E_lattice (or an equivalent elastic/lattice energy scale) for Nb3Cl8 and Nb3Br8 sourced from a
  cited primary reference (not a search-engine summary — verify against the primary text per the
  chem-jiy lesson)": **met** — `arXiv:2412.03360v1`, fetched via `WebFetch`, read in full via `Read`
  on the saved PDF, exact numbers quoted above and in `nb3x8_tc_attribution.NB3X8_EXFOLIATION`.
- "the pre-registered F_phase(T) crossing check implemented and run": **met** —
  `required_e_lattice`/`candidate1_verdict` in `nb3x8_tc_attribution.py`, gated by
  `tests/test_nb3x8_tc_attribution_spec.py::test_G5_*`/`test_G6_*`/`test_G7_*`.
- "kill condition evaluated and recorded either way in specs/SPEC_nb3x8_tc_attribution.md and
  specs/BACKLOG.md": **met** — both updated; verdict is SURVIVES, recorded with the caveats above,
  not rounded into an overclaim.

## Files touched (for the diff review)

- Brought forward (new in this branch, authored by chem-phh): `nb3x8_tc_attribution.py`,
  `specs/SPEC_nb3x8_tc_attribution.md`, `tests/test_nb3x8_tc_attribution_spec.py`.
- Modified (chem-phh's prior edits, now extended by this session):
  `nb3x8_gaps.py`, `specs/SPEC_nb3x8_magnetometry.md`, `specs/BACKLOG.md`.
- `.beads/issues.jsonl` — passive export, updated by `bd update --claim` / `bd close`.
- Not touched: `.claude/settings.local.json` (pre-existing untracked local file, left alone).

## Suggested next commands (git policy: I commit, I do not push)

```
git add nb3x8_tc_attribution.py specs/SPEC_nb3x8_tc_attribution.md tests/test_nb3x8_tc_attribution_spec.py
git commit   # chem-phh prerequisite files, brought forward verbatim

git add nb3x8_gaps.py specs/SPEC_nb3x8_magnetometry.md specs/BACKLOG.md
git commit   # chem-phh's prior edits to shared files, also brought forward verbatim

git add nb3x8_tc_attribution.py tests/test_nb3x8_tc_attribution_spec.py specs/SPEC_nb3x8_tc_attribution.md specs/SPEC_nb3x8_magnetometry.md specs/BACKLOG.md nb3x8_gaps.py
git commit   # chem-cfl's own candidate-1 implementation, gates, and spec/backlog updates

git add sandbox-handoffs/chem-cfl.md .beads/issues.jsonl
git commit   # handoff + beads claim/close state
```
