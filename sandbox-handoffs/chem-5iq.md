# chem-5iq handoff

## What this bead was

chem-tn4 fixed conda/pip instructions in `run_in_chem.sh`, `.claude/skills/spec-invent/SKILL.md`,
`prd.md`, and the four certkit robot gates (+ `certkit_bridge.py`) -- but that work landed on
branch `sandbox/chem-tn4` (commit `74a1b07`), which is **not yet merged into the base this branch
was cut from**, so those files still show conda in this working tree. That's expected and out of
scope here; I did not touch them. chem-5iq is the repo-wide grep for everything chem-tn4's named
scope did *not* cover.

## What I changed, and why

### Conda/bare-pip instructions converted to uv (acceptance criterion 1: "every
README/setup-script instruction that claims to set up or run this repo's main (non-chem-ft)
environment uses uv")

- `README.md:109-110` -- `conda create ... && conda activate` + `pip install -r requirements.txt`
  → `uv sync --extra dmrg --extra test --extra dev`.
- `README.md:126` -- the materials-path quickstart's bare `python run_hybrid_solver.py ...` →
  `uv run python run_hybrid_solver.py ...`. Not separately named in the bead but directly adjacent
  to the install-instructions fix and the same category of problem (CLAUDE.md: "Never invoke a
  bare `python` or `python3` command").
- `submit_hipergator_job.sh:38-43` -- replaced `module load conda/24.1.2` + `conda activate chem`
  with: ensure `uv` is on PATH (install via the official installer if missing), `uv sync --extra
  dmrg --extra test --extra dev`, then `source .venv/bin/activate`. I kept the final
  `source .venv/bin/activate` (rather than switching every downstream call to `uv run`) because
  the script it invokes, `run_sprint_benchmarks.sh`, calls bare `python` internally (line 81) and
  that file is not in this bead's named scope -- activating the venv makes that still resolve
  correctly without editing a file outside scope.
  - Caveat I could not verify: I have no HiPerGator access, so I could not actually submit this
    job. `bash -n submit_hipergator_job.sh` is clean (syntax only).
- `logs/run_hchain_headline.sh`, `logs/run_nbn_schedules.sh` -- `conda run -n chem python -u` →
  `uv run python -u` (5 occurrences total, verbatim `sed` substitution, diff reviewed).
- `run_grand_sweep_on_cloud.sh:62` -- `conda run -n chem python run_grand_sweep.py "${API_URL}"` →
  `uv run python run_grand_sweep.py "${API_URL}"`.

### Scattered `pip install` asides (acceptance criterion 2: "reviewed and either converted to
uv-native wording or left with a stated reason")

Converted (the dependency is already declared in `pyproject.toml`, so the install hint should name
the repo's actual install path instead of a generic `pip install <pkg>`):

- `benchmark_dmrg.py`, `benchmark_dmrg_large.py`, `benchmark_nbn.py`, `benchmark_hchain_tdl.py`,
  `hybrid_quantum_solver/dmrg_reference.py` (4 call sites) -- `pip install block2` →
  `uv sync --extra dmrg` (block2 is the `dmrg` extra).
- `hybrid_quantum_solver/quantum_krylov_solver.py:192` -- `pip install cupy-cuda12x` →
  `uv sync --extra gpu` (cupy-cuda12x is the `gpu` extra). Note CLAUDE.md tells *this sandbox* not
  to install the gpu extra (no GPU target here) -- that's a sandbox-environment instruction, not a
  statement that the extra shouldn't exist or be referenced; the error message is for a user who
  actually has an NVIDIA GPU.
- `run_nbn_sqd_sweep.py:35` -- `pip install qiskit-addon-sqd (inside the chem conda env)` →
  "It is a core dependency -- run `uv sync`" (qiskit-addon-sqd is already in `pyproject.toml`'s
  top-level `dependencies`, not even an extra -- the old message was stale/misleading on two
  counts, conda AND implying it's optional).
- `run_hardware_odmd.py:160` -- same fix, same reasoning, for `qiskit-ibm-runtime` (also a
  top-level dependency already).

Left as-is, with the reason stated here (and, for the one in-repo case, in the file itself):

- `ft_resource_estimator.py:22-23` -- `pip install openfermion pytest`. Confirmed by hard import
  (`from openfermion.resource_estimates import df, sf` at module level, no try/except) that this
  script genuinely cannot run in the uv-managed main env -- `openfermion` is not in
  `pyproject.toml`. This is exactly the `chem-ft` conda env the bead itself called out as
  "intentionally separate ... not in pyproject.toml's dependencies." I updated the docstring to
  point at `bash setup_chem_ft.sh` (the actual repo mechanism) instead of a bare, incomplete `pip
  install openfermion pytest` line, but did not convert it to uv. I verified
  `uv run python -c "import ft_resource_estimator"` fails with `ModuleNotFoundError: No module
  named 'openfermion'` in this sandbox, confirming the env split is real, not just asserted.
  - I also ran the feasibility check the bead asked for before any conversion: `uv pip compile
    pyproject.toml <(echo openfermion) --extra dmrg --extra test` resolves cleanly (openfermion
    1.8.1 + scipy 1.15.3, satisfying this repo's `scipy<1.16` pin) -- so openfermion's dependency
    chain *would* resolve under uv if someone wanted to fold it in later. That's a separate
    architectural decision (declaring an `ft`/`openfermion` extra and retiring
    `setup_chem_ft.sh`'s conda clone), which I did not make unilaterally -- it changes the
    documented env-isolation architecture and deserves its own bead if the owner wants to pursue
    it. Did not file that bead myself since it's speculative, not a found defect.
- `setup_chem_ft.sh` itself -- untouched. Same chem-ft rationale; this is the one script the bead
  explicitly flagged as possibly-intentionally-separate and asked me to confirm rather than
  convert. Confirmed via the feasibility check above and left alone.
- `pytest_invariants/README.md:26-28` -- `pip install -e .` / `pip install pytest-invariants`.
  This is a generically-packaged pytest plugin whose README describes both in-repo use and
  eventual standalone PyPI publication (the file's own text: "This plugin was extracted from a
  hybrid quantum-chemistry solver ... or, once published: pip install pytest-invariants"). This is
  the acceptance criterion's own example case ("generic pip-installable library note, not a repo
  setup step"). I added one clarifying paragraph above the existing instructions noting that
  within this repo `uv sync` already installs it, without touching the generic pip install
  guidance that's correct for external reuse.
- `specs/tasks/01-certchem-m1.md`, `05-invariant-plugin.md`, `12-odmd-anywhere.md` -- `pip install
  -e . && pytest`, `pip install pytest-invariants`, `pip install` (ODMD-Anywhere). Read
  `specs/tasks/00-MASTER-LIST.md`: this directory is a historical portfolio task-breakdown
  (20 mini-specs, several entirely outside chemistry -- e.g. #13 Leyning Coach / #14 Tanakh
  Corpus are Biblical-Hebrew cantillation projects). #5 and #12's `pip install` lines are
  literally the definition-of-done for *hypothetical future standalone PyPI packages*
  (pytest-invariants, ODMD-Anywhere) -- not instructions for setting up this repo. #1's line is a
  historical milestone (`tag v0.1.0`) for `certchem`, which is now folded into this repo's own
  `pyproject.toml` package list; rewriting a historical done-criterion to retroactively say "uv"
  would misrepresent what was actually run at the time. Left all three unchanged. This directory
  is not part of CLAUDE.md's live SDD loop (that's `specs/SPEC_*.md` + `tests/test_*_spec.py`;
  `specs/tasks/*.md` is a separate, older planning artifact).

## Filed, not fixed (new bead, out of this bead's scope)

`chem-5fl` -- ~16 top-level scripts' docstrings still say `Run:  python <script>.py` (bare
python), contradicting CLAUDE.md's "never invoke a bare python/python3 command" rule. Found while
auditing but the original bead named only conda/pip wording, not bare-python invocation, so I
filed it separately rather than doing it here. List of files is in the bead body.

## Verification run

```
$ uv sync --extra dmrg --extra test --extra dev
Resolved 118 packages ... Installed 8 packages
$ uv run ruff check .
All checks passed!
$ bash -n submit_hipergator_job.sh logs/run_hchain_headline.sh logs/run_nbn_schedules.sh run_grand_sweep_on_cloud.sh
(all four: syntax OK, no output)
$ uv run pytest tests/ --ignore-glob='tests/test_*_spec.py' --ignore=tests/test_dmrg_reference.py -q
56 passed, 4 skipped, 25 warnings in 51.89s
$ uv run pytest tests/test_dmrg_reference.py -q      # separate process, per block2/pyscf isolation rule
2 passed, 1 skipped in 2.79s
$ uv run python -c "import benchmark_dmrg, benchmark_dmrg_large, benchmark_nbn, benchmark_hchain_tdl, run_hardware_odmd; import run_nbn_sqd_sweep"
(imports clean)
$ uv run python -c "import ft_resource_estimator"
ModuleNotFoundError: No module named 'openfermion'     # expected -- confirms the chem-ft env split is real
```

Did not run `make gates` (full spec-gate suite, one-process-per-file, can run long per CLAUDE.md's
"hour-long jobs" warning) -- this bead changed only string literals in error/help messages and
documentation, touched zero logic, and the full non-DMRG suite plus the isolated DMRG reference
test both pass clean. If a reviewer wants the full gate run for extra assurance, run
`scripts/run_gates.sh` (or `make gates`); I'd expect it green but have not verified it myself.

Hardware for the above: inside the sandbox container (`uname -a`: Linux, x86_64 per
`uv 0.12.19 (x86_64-unknown-linux-musl)`); every command above ran in well under a minute, so no
wall-clock/resource report is owed per the "report hardware and wall time" rule (that's for
anything over a minute).

## Pre-registered acceptance criteria, checked against what's above

1. "Every README/setup-script instruction that claims to set up or run this repo's main
   (non-chem-ft) environment uses uv, not conda/bare pip." -- PASS for every file in the bead's
   named list (README.md, submit_hipergator_job.sh, logs/run_hchain_headline.sh,
   logs/run_nbn_schedules.sh, run_grand_sweep_on_cloud.sh). setup_chem_ft.sh is explicitly
   excluded by the "(non-chem-ft)" qualifier and was left alone on purpose.
2. "Package-install asides for genuinely optional third-party deps ... are reviewed and either
   converted to uv-native wording or left with a stated reason." -- PASS: every file in the bead's
   named list was reviewed; conversions and stated-reason decisions are itemized above.

## What I decided not to do, and why

- Did not touch `run_in_chem.sh`, `.claude/skills/spec-invent/SKILL.md`, or `prd.md` -- these are
  chem-tn4's named scope, already fixed on `sandbox/chem-tn4` (not yet merged here); re-touching
  them here would conflict with that branch's own commit.
- Did not fold openfermion/`ft_resource_estimator.py`/`setup_chem_ft.sh` into uv, despite
  confirming it's technically resolvable -- that's an architecture decision outside a docs-wording
  audit bead, see above.
- Did not edit `requirements.txt` -- nothing in this bead's scope references it anymore after the
  README fix; it's still read by `run_in_chem.sh` (chem-tn4's unmerged branch) and mentioned in
  `REFACTOR_PLAN.md`. Removing/rewriting it is a separate decision for whoever reconciles the two
  branches.
- Filed (did not fix) the bare-`python` docstring issue as `chem-5fl` -- out of this bead's named
  scope.

## What I could not verify

- `submit_hipergator_job.sh` actually working on a real HiPerGator node -- no cluster access from
  this sandbox. Verified syntax only (`bash -n`).
- `run_grand_sweep_on_cloud.sh` against real GCP/terragrunt infra -- no cloud credentials here.
  Verified syntax only.
- The full `make gates` spec-gate suite -- not run (see Verification run section for why, and
  what I ran instead).

## Suggested commands for the reviewer

```bash
git status
git diff -- README.md submit_hipergator_job.sh logs/ run_grand_sweep_on_cloud.sh \
  benchmark_dmrg.py benchmark_dmrg_large.py benchmark_nbn.py benchmark_hchain_tdl.py \
  hybrid_quantum_solver/dmrg_reference.py hybrid_quantum_solver/quantum_krylov_solver.py \
  ft_resource_estimator.py run_nbn_sqd_sweep.py run_hardware_odmd.py pytest_invariants/README.md
git add README.md submit_hipergator_job.sh logs/run_hchain_headline.sh logs/run_nbn_schedules.sh \
  run_grand_sweep_on_cloud.sh benchmark_dmrg.py benchmark_dmrg_large.py benchmark_nbn.py \
  benchmark_hchain_tdl.py hybrid_quantum_solver/dmrg_reference.py \
  hybrid_quantum_solver/quantum_krylov_solver.py ft_resource_estimator.py run_nbn_sqd_sweep.py \
  run_hardware_odmd.py pytest_invariants/README.md
git commit -m "docs: convert remaining conda/bare-pip instructions to uv (chem-5iq)"
```
