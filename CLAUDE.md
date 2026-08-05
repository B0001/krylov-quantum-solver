# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A hybrid quantum–classical pipeline for molecular ground-state energies. The validated core path
is **real-time quantum Krylov subspace diagonalization**: PySCF extracts active-space integrals →
a vetted Jordan–Wigner transform builds a qubit Hamiltonian → the Krylov solver estimates the
ground state, validated against FCI/CASCI and DMRG.

This repo was **rebuilt after a scientific audit** found the original physics core to be entirely
broken (see `REFACTOR_PLAN.md`). The defining culture here is *falsifiable honesty*: every claimed
result must have a reference (FCI, DMRG, experiment, or analytic limit) and a test that could kill
it. When adding capability, follow this — do not introduce numbers without a ground-truth check.

## Environment & commands

## Environment Rules

- **Package Manager:** `uv` is the exclusive dependency manager. `pip`, `poetry`, and `conda` are strictly prohibited.
- **Dependency Installation:** Use `uv add <package>` for new dependencies and `uv remove <package>` for removal. Never invoke `pip install`.
- **Script Execution:** Execute all Python scripts and commands via `uv run` (e.g., `uv run python script.py`). Never invoke a bare `python` or `python3` command.

```bash
make gates     # run every tests/test_*_spec.py, each in its OWN process (see segfault note)
make test      # full suite (block2/DMRG tests isolated into a second process)
make lint      # ruff check .
```

Run a single test:
```bash
conda run -n chem python -m pytest tests/test_krylov_convergence.py -v
conda run -n chem python -m pytest tests/test_krylov_convergence.py::test_name -v
```

`bash run_in_chem.sh` does the full validation walkthrough (versions → deps → tests → benchmarks).

### Critical: block2/DMRG process isolation

`block2` (the DMRG backend) initializes its own OpenMP runtime and **segfaults if it loads into a
process that already imported `pyscf` or `qiskit-aer`**. This is why the Makefile and `run_in_chem.sh`
run the DMRG/spec tests (`tests/test_dmrg_reference.py` + every `tests/test_*_spec.py`) in a
*separate* pytest process from the rest of the suite, and why `make gates` runs each spec gate in
its own process. Preserve this isolation when adding tests or build steps — the `test_*_spec.py`
glob keeps it correct automatically as new spec gates are added.

### Dependency pin

`scipy` is pinned `>=1.8,<1.16`: scipy 1.16–1.17 crash the qiskit-nature import chain (a vendored
`array_api_compat` torch-detection bug). Do not bump it past 1.15.x. DataFrames standardized on
**polars** (not pandas).

## Spec-driven development (SDD)

New capabilities are added through the loop in `specs/README.md`, not ad hoc:

1. Pick/add a hypothesis in `specs/BACKLOG.md` — a *claim + a cheap check that could disprove it*.
2. Copy `specs/SPEC_TEMPLATE.md` → `specs/SPEC_<slug>.md` (goal, interface, acceptance gates,
   out-of-scope, honest caveats).
3. Write `tests/test_<slug>_spec.py` encoding the gates (initially failing).
4. Implement the **minimum** code to pass — reuse validated primitives over new code.
5. `make gates` until green. **If a gate proves unsatisfiable, revise the spec and record why** —
   that mismatch is the finding, not a failure (see `SPEC_hchain_tdl.md` for a worked example, and
   the `[-]` killed entries in `BACKLOG.md`).

Conventions: specs are `specs/SPEC_<slug>.md`, gate tests are `tests/test_<slug>_spec.py`. Every
result is checked against a reference; no reference → not yet a spec.

## Architecture

### `hybrid_quantum_solver/` — the package (validated live path)

The public API is exported from `__init__.py`. The pipeline flows:

- `chemistry_gateway.py` — PySCF integral extraction (CIF/geometry → CASCI active space).
- `molecular_hamiltonian.py` — vetted Jordan–Wigner qubit Hamiltonian (`SparsePauliOp`) + the
  Hartree–Fock reference state; nuclear/core constant tracked as an energy offset.
- `quantum_krylov_solver.py` — real-time quantum Krylov: builds |φₖ⟩ = e^(−ikΔtH)|φ_HF⟩, forms
  Hermitian H and S, solves `Hc = ESc` by thresholded canonical orthogonalization. Variationally
  bounded below by the true ground state (a key correctness invariant the tests rely on).
- `trotter_krylov.py` — Trotter-circuit Krylov + qiskit-aer device-noise expectation path.
- `hardware_krylov.py` — on-hardware Krylov (Hij/Sij via ancilla Hadamard tests).
- `noise.py` — shot-noise model (added Hermitian-symmetrically), Aer `NoiseModel`, ZNE.
- `pipeline.py` — end-to-end `run_geometry` / `run_from_integrals` → `PipelineResult`.
- `dmrg_reference.py` — classical reference: `reference_energy(method="auto")` uses DMRG (block2)
  when installed, falls back to exact FCI. The integral convention here is shared with
  `build_hamiltonian_from_integrals` and pinned by the test suite — keep them in sync.

**Quarantined / do not use:** `orchestrate_hybrid_pipeline.py` and `quantum_sampler.py` are the old
broken core, retained only as regression fixtures. They are intentionally not exported and excluded
from ruff. Do not import, "fix", or build on them.

### Top-level scripts

These are studies and harnesses, not library code. Roughly:

- **Near-term benchmarks/validation:** `benchmark_krylov.py` (FCI table), `benchmark_n2.py` (N₂
  dissociation vs CASCI), `benchmark_resources.py` (circuit cost), `benchmark_dmrg*.py`,
  `benchmark_nbn.py`, `benchmark_hchain_tdl.py`, `cross_check.py` (4 independent methods must agree),
  `validate_and_cost.py` (integrals → costed cross-validated verdict).
- **Fault-tolerant stack:** `qubitization_blueprint.py`, `qpe_walk_readout.py`, `iterative_qpe.py`,
  `lambda_ladder.py` / `df_factorization.py` / `ft_resource_estimator.py` (Hamiltonian
  factorization → qubitization 1-norm λ → T-gate budget), `taper_qubits.py` (Z2 tapering).
- **Variational/other:** `adapt_vqe.py`, `krylov_subspace_solver.py`, SQD sweeps
  (`run_nbn_sqd_sweep.py`, `smoke_test_sqd_plumbing.py`).
- **Materials path:** `run_hybrid_solver.py --input_file …cif --active_space 8,8`. Caveat: the CIF
  path builds a *finite molecular cluster* with no periodic boundary conditions — it is **not** a
  periodic-solid calculation. Transition-metal systems are a research target, not a validated result.

### `tests/` — three validation gates plus spec gates

`test_reference_energies.py` (qubit Hamiltonian == FCI), `test_krylov_convergence.py` (converges to
FCI, respects the variational floor), `test_noise_resilience.py` (shot noise bounded, improves with
shots), plus `test_*_spec.py` SDD gates and `test_dmrg_reference.py`.

## Long-running jobs

Hour-long DMRG/benchmark jobs die on session teardown if launched from here — have the user run
those in their own terminal. block2 DMRG scratch goes to `.dmrg_tmp/` (git-ignored).


<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:6cd5cc61 -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

**Architecture in one line:** issues live in a local Dolt DB; sync uses `refs/dolt/data` on your git remote; `.beads/issues.jsonl` is a passive export. See https://github.com/gastownhall/beads/blob/main/docs/SYNC_CONCEPTS.md for details and anti-patterns.

## Agent Context Profiles

The managed Beads block is task-tracking guidance, not permission to override repository, user, or orchestrator instructions.

- **Conservative (default)**: Use `bd` for task tracking. Do not run git commits, git pushes, or Dolt remote sync unless explicitly asked. At handoff, report changed files, validation, and suggested next commands.
- **Minimal**: Keep tool instruction files as pointers to `bd prime`; use the same conservative git policy unless active instructions say otherwise.
- **Team-maintainer**: Only when the repository explicitly opts in, agents may close beads, run quality gates, commit, and push as part of session close. A current "do not commit" or "do not push" instruction still wins.

## Session Completion

This protocol applies when ending a Beads implementation workflow. It is subordinate to explicit user, repository, and orchestrator instructions.

1. **File issues for remaining work** - Create beads for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **Handle git/sync by active profile**:
   ```bash
   # Conservative/minimal/default: report status and proposed commands; wait for approval.
   git status

   # Team-maintainer opt-in only, unless current instructions forbid it:
   git pull --rebase
   git push
   git status
   ```
5. **Hand off** - Summarize changes, validation, issue status, and any blocked sync/commit/push step

**Critical rules:**
- Explicit user or orchestrator instructions override this Beads block.
- Do not commit or push without clear authority from the active profile or the current user request.
- If a required sync or push is blocked, stop and report the exact command and error.
<!-- END BEADS INTEGRATION -->
