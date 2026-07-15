# SPEC: nb_validate_sweep — the ragged-schema CSV union, checked not just asserted in a docstring

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

> **A more consequential finding surfaced while probing this spec, fixed alongside it:** this
> repo's `pyproject.toml` had no `[tool.pytest.ini_options]` section, so it didn't qualify as a
> pytest config file — meaning pytest's rootdir search silently walked UP past this repo into its
> parent directory, found an UNRELATED sibling project's `pyproject.toml` there (one with a
> qualifying `[tool.pytest.ini_options]`), and loaded THAT project's `conftest.py` — which force-
> stubs `sys.modules["polars"]` globally (`_stub("polars", col=MagicMock())`, for that other
> project's own test isolation) before any test file even runs. Every `pytest` invocation from
> within this repo — including `make gates`, which shells out to exactly the same `python -m
> pytest <file> -q` pattern — was silently vulnerable to any gate that calls real `polars`
> functionality beyond what the stub's single mocked attribute (`col`) covers; this went unnoticed
> only because no gate before `test_nb_validate_sweep_schema_spec.py` happened to construct a real
> `polars.DataFrame`. Fixed by adding a minimal `[tool.pytest.ini_options]` (with `testpaths =
> ["tests"]`, itself a correctness improvement) to `pyproject.toml`, anchoring rootdir at this repo
> and confirmed to stop the sibling `conftest.py` from loading at all. See the BACKLOG entry and PR
> for the full regression evidence (54 pre-existing gates re-run clean after the fix).

## 1. Goal

`nb_validate_sweep.py` wires `validate_and_cost` into `run_nbn_sqd_sweep.py`'s spin-sector
machinery, writing one CSV row per physically reachable spin sector so "each Nb3 sector carries a
validated energy and an FT cost tag in your telemetry." Never gated — only exercised by its own
`__main__` demo, printed and eyeballed. Its `_to_frame` helper carries a precise, checkable claim in
its own docstring: rows may be RAGGED (an "OK" sector emits ~19 fields, a guarded failure emits only
`{mol_spin, multiplicity, status}`), and `_to_frame` "union[s] the keys and null-fill[s] the gaps
... while preserving first-seen column order." That is exactly the kind of claim that is easy to get
subtly wrong (a `set()` instead of `dict.fromkeys()` would silently scramble column order) and had
never been checked. This spec gates the sweep end-to-end and the ragged-schema claim specifically.
False if the union/null-fill/ordering claim doesn't hold, or if the sweep's own two spin sectors
don't reproduce independently-computed CASCI references.

## 2. Background and honest framing

- `nb_validate_sweep.py` reuses `run_nbn_sqd_sweep.py` and `validate_and_cost.py` unmodified — no
  new physics, only falsifiers around the sweep/CSV-flattening layer.
- **What you can claim if the gates pass:** on LiH CAS(2,2) (the module's own `__main__` example —
  small and safe; see R1 for why this matters), the sweep produces exactly the two physically
  reachable spin sectors (singlet, triplet), both `status == "OK"`, with `reference_energy` matching
  independently-computed CASCI to `< 1e-6` Ha for each sector; every OK row's `ft_status ==
  "no_openfermion"` in the standard `chem` env (no `openfermion`), the graceful-degradation path
  through the CSV-flattening layer, not just `validate_and_cost`'s own return dict
  (`SPEC_validate_and_cost_composition.md` already checked the latter); `_to_frame`'s ragged-schema
  handling is exactly what its docstring claims — union of keys, null-filled gaps, first-seen column
  order — verified on synthetic mixed-schema rows matching the sweep's own OK/FAILED shapes; and the
  sweep actually WRITES a valid, re-readable CSV to disk (an end-to-end I/O check, not just the
  in-memory row list).
- **What you cannot claim:** anything about a real CIF-loaded structure or the FT-cost numbers
  themselves (both out of scope — `load_geometry` needs a real CIF file this repo doesn't bundle,
  per `SPEC_nb3x8_gaps.md`'s own honest scope note, and FT-cost needs the `chem-ft` env, per
  `SPEC_validate_and_cost_composition.md`'s R1); that the per-sector guarded-failure path (`try/except`
  around `build_scf`/`integrals_for_spin`/`validate_and_cost`) is exercised by a REAL failing SCF
  here — both tested sectors converge cleanly, so that path is covered only via `_to_frame`'s
  synthetic ragged-row test (G3), not a genuine SCF non-convergence.
- **Reference:** independently-computed CASCI energies (dense diagonalization, singlet via RHF,
  triplet via UHF) for the accuracy claims; direct inspection of the returned rows and the written
  CSV's columns/values for the schema and I/O claims.

## 3. Approach

Reuse `valid_spin_sectors`, `validate_sweep`, `_to_frame` unmodified. G1/G2/G4: call
`validate_sweep("Li 0 0 0; H 0 0 1.6", "sto-3g", {}, cas_electrons=2, cas_orbitals=2,
output_csv=<pytest tmp_path>)` — the module's own `__main__` system, writing to a pytest-managed
temp path (never the repo's own directory) so no stray CSV is left behind. G3: construct synthetic
mixed-schema row dicts directly (an "OK"-shaped row with the sweep's own ~19 keys, a "FAILED"-shaped
row with only 3) and call `_to_frame` on them directly — no need to force a real SCF failure.

## 4. Public interface

No new library code — this spec adds only test-file assertions around `nb_validate_sweep.py`'s
existing `validate_sweep`/`_to_frame`, reused unchanged.

## 5. Acceptance criteria (validation gates)

- **G1 — the sweep produces the correct sectors with correct energies.** On LiH CAS(2,2):
  `valid_spin_sectors(2, 2) == [0, 2]` (singlet, triplet); the sweep returns exactly 2 rows, both
  `status == "OK"`; `reference_energy` for `mol_spin=0` matches an independently-computed RHF/CASCI
  singlet (`-7.8621288...` Ha) to `< 1e-6`, and for `mol_spin=2` matches an independently-computed
  UHF/CASCI triplet (`-7.7663313...` Ha) to `< 1e-6`.
- **G2 — graceful FT degradation through the CSV-flattening layer, not just `validate_and_cost`'s
  own dict.** Every OK row has `ft_status == "no_openfermion"` and every `ft_*` numeric field is
  `None` in the standard `chem` env (confirmed via `validate_and_cost._HAVE_FT is False`).
- **G3 — THE FINDING (definition of done): the ragged-schema union is exactly what the docstring
  claims.** `_to_frame` on a synthetic `[OK-shaped row, FAILED-shaped row]` list produces a
  DataFrame whose columns are the UNION of both rows' keys, in FIRST-SEEN order (all of the
  OK row's ~19 keys first, in their original order, then any FAILED-only keys appended); the
  FAILED row's missing fields are null (not absent, not crashed, not misaligned with the OK row's
  values).
- **G4 — the sweep actually writes a valid, re-readable CSV.** After `validate_sweep(...,
  output_csv=<tmp path>)`, the file exists, has exactly 2 data rows when re-read with
  `polars.read_csv`, and its columns match `_to_frame(rows).columns` exactly — the on-disk
  artifact the sweep's own docstring promises ("each Nb3 sector carries ... in your telemetry"),
  not just the in-memory return value.

> Definition of done: **G3**. G1/G2/G4 confirm the sweep works end-to-end on its own known-safe
> example; G3 is the precise, previously-unchecked docstring claim about the CSV schema itself.

## 6. Implementation plan (test-first)

1. Write `tests/test_nb_validate_sweep_schema_spec.py` encoding G1-G4 (RED in the sense these checks
   are new, even though `nb_validate_sweep.py`'s functions are not).
2. No changes to `nb_validate_sweep.py` — every gate calls its existing public functions directly.
3. Targeted pytest to green; ruff clean.

## 7. Out of scope

- Real CIF-loaded structures (`load_geometry`/`from_cif`) — no CIF bundled in this repo, per
  `SPEC_nb3x8_gaps.md`'s own scope note.
- FT-cost numbers themselves (needs `chem-ft`, out of scope per `SPEC_validate_and_cost_composition.md`).
- A genuinely failing SCF sector — the per-sector `try/except` guard is covered only at the
  `_to_frame` schema level (G3), not by forcing a real non-convergence.
- Systems beyond LiH CAS(2,2) — a natural follow-up, not attempted here. **Explicitly does not**
  attempt a larger CAS that would exceed `validate_and_cost`'s default `qubit_dense_max_orb=7` —
  `SPEC_validate_and_cost_composition.md`'s R2 already recorded the OOM risk of pushing taper past
  its tractable range; this spec stays at the module's own small, known-safe example.

## 8. Caveats and risks

- **R1 — staying at LiH CAS(2,2) is a deliberate safety choice, not an oversight.** A larger CAS run
  through this sweep (as the module's real intended use, Nb3 sectors) would call `validate_and_cost`
  per sector, which internally taps `taper_hamiltonian` — the same exponential-in-qubit-count cost
  the immediately-prior spec's probing (`SPEC_validate_and_cost_composition.md` R2) OOM-killed a
  process over. This spec never raises `cas_orbitals` past the module's own tiny demo value.
- Honest limitation: one system, two sectors, no real CIF, no FT-cost numbers.

## 9. Deliverables

- `tests/test_nb_validate_sweep_schema_spec.py` — gates G1-G4 (no library code changes).
- Results summary (with the R1 caveat) in the PR description / BACKLOG entry.
