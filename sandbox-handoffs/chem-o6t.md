# chem-o6t handoff — ruff extend-exclude vs CLAUDE.md quarantine claim

## Status: fixed on branch `sandbox/chem-o6t`; bead left OPEN for the operator to close.

Worked by hand from a cloud session. Per the operator's instructions the claim was already held,
and this session ran no `claims.py` and no `bd` write (`update`, `close`, `dolt push`). Only
`bd show chem-o6t` was run, to read the spec.

## Fix

The acceptance criteria allowed either direction. I chose **config follows docs**: I added
`hybrid_quantum_solver/quantum_sampler.py` to `[tool.ruff] extend-exclude` in `pyproject.toml`
and updated the comment above it. Every doc that mentions the file treats it as quarantined
alongside `orchestrate_hybrid_pipeline.py`: CLAUDE.md, README.md:158, prd.md:103,
`hybrid_quantum_solver/__init__.py:10` and `.claude/skills/spec-invent/SKILL.md:24`
("do not touch"). Correcting CLAUDE.md instead would have left ruff free to demand edits to a
file that five places say is frozen.

No docs changed. CLAUDE.md is now accurate as written.

## Verification

- Reproduced first. Before the change, `uv run ruff check --show-settings
  hybrid_quantum_solver/quantum_sampler.py` listed only `orchestrate_hybrid_pipeline.py` in
  `extend_exclude`, and `ruff check . --show-files` listed `quantum_sampler.py`. After the change,
  both files appear in `extend_exclude`, and `--show-files` lists neither.
- `uv run ruff check hybrid_quantum_solver/ tests/` → `All checks passed!`
- `uv run --extra test python -m pytest tests/ --ignore-glob='tests/test_*_spec.py'
  --ignore=tests/test_dmrg_reference.py -q` (the first half of `make test`) → **55 passed, 4
  skipped** in 27.5 s. This includes `tests/test_reference_energies.py`, the one test file that
  references the legacy fixtures.
- **Not run:** the DMRG/spec half of `make test` and `make gates`. This change touches only ruff
  config, which no test reads, and those gates run for hours.

## Caveats / follow-ups (not fixed here, out of scope)

- **`make lint` (`ruff check .`) is red on `main`, independent of this bead.** There are 17
  errors, all in committed scratch files (`scratch_*_52i.py`,
  `scratch_before_after_hbv.py`) and `scripts/spec_pm3_subspace_eta_bound.py`. None of them are
  in the files this bead touched. The bead's claim that "`ruff check .` is green" held on
  2026-09-28 but no longer does. Worth a bead, which I did not file because `bd create` is a
  write this session was told not to make.
- `extend-exclude` only affects file discovery. `ruff check hybrid_quantum_solver/quantum_sampler.py`
  given as an explicit path still lints it, because `force-exclude` is off. That is ruff's normal
  behavior and matches how `orchestrate_hybrid_pipeline.py` has always been handled.
- `bd` (1.3.0) rewrites `.beads/config.yaml`/`.beads/.gitignore` and drops `.beads.gate.lock`
  at the repo root on first use. None of that is committed on this branch.
