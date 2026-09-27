# chem-4y9 — H-chain EOS: reproduce the published 10-point R grid at n=10 (FCI check)

**Status: CLOSED.** Acceptance criteria met.

## What was asked

`specs/BACKLOG.md`'s "one fixed R was the easy point" item, part (a): the pipeline had only ever
validated the H-chain geometry/basis convention at one bond length, R=1.8 Bohr. Run `fci_energy`
at n=10 across the published Simons Collaboration 10-point R grid (Motta et al., *PRX* 7, 031059,
2017) and check agreement with the published FCI totals to <1e-4 Ha at all ten R, validating the
geometry convention across the whole Mott crossover, not just the one point every other H-chain
spec in this repo assumes. Part (b) of the same backlog item was already closed by chem-1uu.

## What I changed

- `specs/SPEC_hchain_R_grid.md` — new spec (goal, reference, gates G1/G2, caveats).
- `specs/hchain_R_grid_motta_fci.csv` — vendored reference table, hand-transcribed from Motta et
  al.'s own Table II ("Potential energy curve of H10 with the minimal (STO-6G) basis", FCI column,
  energy per atom), with the total (×n=10) computed alongside for direct comparison.
- `tests/test_hchain_R_grid_spec.py` — gate test, G1 (grid agreement) + G2 (canonical/localized
  invariance across all 10 R, generalizing chem-1uu's 3-point check).
- `specs/BACKLOG.md` — updated the "one fixed R was the easy point" entry to record part (a) as
  done alongside part (b) (chem-1uu).

No production code changed: `benchmark_hchain_tdl.integrals(n, R_bohr=...)` and
`hybrid_quantum_solver.dmrg_reference.fci_energy` already existed (the `--R` flag was added by
chem-1uu) and were exactly the primitives this gate needed.

## How the reference numbers were obtained (this is the part worth checking)

No internet route gave a machine-readable table for this benchmark on the first few tries:
- `WebSearch` returned canned "I can't search the web" text for every query in this sandbox —
  it is non-functional here, not just rate-limited. Do not trust its output if you see this again.
- `WebFetch`'s *summarizing* pass on the paper hallucinated a plainly-wrong table once (total
  energies around −73 Ha for H10 — physically impossible; H10/STO-6G total energy is ≈−5.4 Ha).
  **Do not trust WebFetch's prose summary of a numeric table** — it is a small model paraphrasing,
  not a parser.
- What worked: `WebFetch` on `https://arxiv.org/pdf/1705.01608v2` saved the actual PDF bytes to
  `~/.claude/projects/-workspace/*/tool-results/*.pdf`, which I then read directly with
  `uv run --with pypdf python` (ephemeral dependency via `uv run --with`, nothing added to
  `pyproject.toml`/`uv.lock`) and located the correct table by text search rather than trusting a
  summary. arXiv ID 1705.01608 itself was found via the `export.arxiv.org` Atom search API
  (`ti:"many-electron problem" AND au:Motta`) after a guessed ID (1701.00054) turned out to be an
  unrelated paper — guessed arXiv IDs are not reliable, verify via the search API or don't guess.

Table II of Motta et al. (*PRX* 7, 031059, 2017; arXiv:1705.01608) gives FCI as its own column, at
all 10 R values, energy **per atom**, 6 decimals. The paper states "DMET[5], MRCI and MRCI+Q
energies coincide with FCI to within 10⁻⁶" — i.e. multiple independent method families agree with
FCI at this small system size, so this is solid ground truth, not one method's output. Internal
sanity check before trusting the transcription: the paper's prose separately states the exact
equilibrium point R_e=1.786 Bohr, E0=−0.542457 Ha/atom; Table II's R=1.8 row gives FCI=−0.542439,
consistent with being one small step off the true minimum. That cross-check is what gave me
confidence the OCR/extraction was reading the right column, not an adjacent one.

## Command whose output justifies the closing claim

```
cd /workspace && uv run python -c "
from benchmark_hchain_tdl import integrals
from hybrid_quantum_solver.dmrg_reference import fci_energy
R_grid = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6]
motta_fci_per_atom = {1.0:-0.382439,1.2:-0.476638,1.4:-0.520509,1.6:-0.538436,1.8:-0.542439,
                      2.0:-0.538963,2.4:-0.522794,2.8:-0.505024,3.2:-0.491038,3.6:-0.481870}
for R in R_grid:
    h1, eri, ne, ec, e_hf = integrals(10, localize=False, R_bohr=R)
    e_fci = fci_energy(h1, eri, ne, ec)
    print(R, e_fci, e_fci - 10*motta_fci_per_atom[R])
"
```

Output (our FCI total, diff vs published × 10, Ha):

| R (Bohr) | our FCI total | Motta total (×10) | diff |
|---|---|---|---|
| 1.0 | −3.82438879 | −3.824390 | +1.21e-6 |
| 1.2 | −4.76637746 | −4.766380 | +2.54e-6 |
| 1.4 | −5.20509432 | −5.205090 | −4.32e-6 |
| 1.6 | −5.38436053 | −5.384360 | −5.25e-7 |
| 1.8 | −5.42438538 | −5.424390 | +4.62e-6 |
| 2.0 | −5.38962606 | −5.389630 | +3.94e-6 |
| 2.4 | −5.22793643 | −5.227940 | +3.57e-6 |
| 2.8 | −5.05024496 | −5.050240 | −4.96e-6 |
| 3.2 | −4.91038284 | −4.910380 | −2.84e-6 |
| 3.6 | −4.81870080 | −4.818700 | −7.99e-7 |

Max |diff| = 5.25e-6 Ha, 20x inside the 1e-4 Ha gate at every point, and the sign of the diff
flips unsystematically across R (consistent with it being the paper's own 6-decimal print
rounding, not a real disagreement).

## Gate run (verbatim)

```
$ cd /workspace && uv sync --extra dmrg --extra test --extra dev
Resolved 106 packages ... Installed 7 packages (block2, pytest, ruff, ...)

$ uv run python -m pytest -q tests/test_hchain_R_grid_spec.py -v
tests/test_hchain_R_grid_spec.py::test_G1_fci_matches_published_grid PASSED
tests/test_hchain_R_grid_spec.py::test_G2_localized_invariant_across_grid PASSED
2 passed in 74.32s

$ GATE_JOBS=1 GATE_GLOB="tests/test_hchain_R_grid_spec.py" bash scripts/run_gates.sh
gates: 1 files, 1 parallel processes, cache=on
PASS  tests/test_hchain_R_grid_spec.py  (72s)
all spec gates passed

$ uv run ruff check .
All checks passed!
```

I did **not** re-run the full `make gates`/`make test` suite (dozens of files, many DMRG, each its
own process — expensive, and this bead touches nothing outside the two new files + BACKLOG.md
prose). The new gate file is pure PySCF (no block2 import), confirmed by grep of
`hybrid_quantum_solver/dmrg_reference.py` showing `import block2` is deferred inside function
bodies, not at module load — so it needs no process isolation and cannot be the thing that
segfaults an isolated run.

## Pre-registered criteria (from the bead) and outcome

- FCI energy at n=10 matches the published reference to <1e-4 Ha at all 10 grid points — **PASS**
  (max 5.25e-6 Ha, ~20x margin).
- Numbers vendored into a tracked file a gate reads — **PASS**
  (`specs/hchain_R_grid_motta_fci.csv`, read by `tests/test_hchain_R_grid_spec.py`).

I additionally pre-registered and passed a G2 (canonical/localized FCI invariance at all 10 R,
not just the 3 chem-1uu checked) as a cheap extra falsifier — not required by the bead, but free
given the same `integrals()` calls, and it generalizes an existing pattern
(`test_hchain_tdl_spec.py::test_G6`) that had previously only been checked at 3 of the 10 points.

## What I decided not to do, and why

- **Did not try to locate the paper's actual data repository** (ref. [11] in the arXiv preprint
  says "Github repository (link to be finalized)" — the placeholder was apparently never resolved
  in this preprint version, and I did not find a working link to the finalized PRX-published
  version's SI). The in-paper Table II was sufficient and is independently cross-checked (see
  above), so chasing the primary data repo would have been effort without changing the verdict.
- **Did not extend to n≠10 or to DMRG at these R values.** The bead is explicitly n=10/FCI-only;
  `SPEC_hchain_largen2.md` §12 already covers DMRG-vs-R mechanism at n=20 for 3 of these R values
  (chem-1uu), and extending the full TDL curve to all 10 R would be a much larger compute effort
  — out of scope, noted as such in the new spec's §7.
- **Did not re-run the whole `make gates` suite.** Justified above; the change surface is two new
  files plus prose in BACKLOG.md, and the new gate itself was run both directly and through
  `scripts/run_gates.sh`'s isolation path.

## What I could not independently verify

- The transcription of Table II from the arXiv PDF was done by extracting text with `pypdf` and
  matching columns to the header row by position — I did not have a second, independent source to
  diff it against (no PRX-hosted SI, no resolved GitHub data repo). The internal consistency check
  against the paper's own prose-stated equilibrium value (R_e, E0) is the only independent check I
  had, and it passed. This is recorded as caveat R1 in `specs/SPEC_hchain_R_grid.md` §8 rather than
  hidden.

## Git state

Tree is dirty and NOT committed, per this session's git policy (do not commit/push). Files touched
by this bead:
```
 M specs/BACKLOG.md
?? specs/SPEC_hchain_R_grid.md
?? specs/hchain_R_grid_motta_fci.csv
?? tests/test_hchain_R_grid_spec.py
```
Everything else showing as modified/untracked in `git status` predates this session (other beads'
uncommitted work) and was not touched.

Suggested commands for a human to commit just this bead's work:
```
git add specs/SPEC_hchain_R_grid.md specs/hchain_R_grid_motta_fci.csv \
        tests/test_hchain_R_grid_spec.py specs/BACKLOG.md
git commit -m "H-chain EOS: validate FCI geometry convention across 10-point Motta R grid (chem-4y9)"
```

## Environment

Linux 7.0.12-linuxkit x86_64, 8 vCPUs (`nproc`). Gate wall time: 72-74s (pure PySCF FCI, n=10,
10 R points, canonical + localized = 30 FCI solves total). No block2/DMRG run for this bead.
