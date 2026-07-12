# Spec-driven development targets. See specs/README.md.
#
# Override the env runner if you don't use conda (e.g. `make gates RUN=`):
#   ENV ?= chem   ->  commands run via `conda run -n chem`
ENV ?= chem
RUN ?= conda run -n $(ENV)

GATE_JOBS ?= $(shell nproc 2>/dev/null || echo 4)

.PHONY: gates gates-seq gates-fresh gates-clean-cache test lint help

help:
	@echo "make gates             - all spec gates, $(GATE_JOBS) parallel processes, cached"
	@echo "make gates-seq         - all spec gates, sequential (GATE_JOBS=1), cached"
	@echo "make gates-fresh       - all spec gates, parallel, ignore cache"
	@echo "make gates-<component> - gates matching tests/test_<component>*_spec.py (e.g. gates-certchem)"
	@echo "make gates-clean-cache - drop all cached gate passes"
	@echo "make test              - full test suite (block2/DMRG tests isolated)"
	@echo "make lint              - ruff check"
	@echo "Overrides: GATE_JOBS=N  RUN='' (no conda)"

# Acceptance gates. Each spec's gate file runs in its OWN process so block2's OpenMP runtime
# never loads into an interpreter that already imported pyscf/qiskit-aer (which segfaults).
# Parallelism is across processes, so that isolation is preserved. Passes are cached against
# a conservative global source digest (any .py change invalidates everything — safe, never
# wrongly skips). See scripts/run_gates.sh.
gates:
	@GATE_RUN="$(RUN)" GATE_JOBS=$(GATE_JOBS) bash scripts/run_gates.sh

gates-seq:
	@GATE_RUN="$(RUN)" GATE_JOBS=1 bash scripts/run_gates.sh

gates-fresh:
	@GATE_RUN="$(RUN)" GATE_JOBS=$(GATE_JOBS) GATE_NO_CACHE=1 bash scripts/run_gates.sh

# Component gates: `make gates-certchem` runs tests/test_certchem*_spec.py, etc.
gates-%:
	@GATE_RUN="$(RUN)" GATE_JOBS=$(GATE_JOBS) GATE_GLOB='tests/test_$**_spec.py' bash scripts/run_gates.sh

gates-clean-cache:
	@rm -rf .gates-cache && echo "gate cache cleared"

# Full suite. The block2/DMRG tests (test_dmrg_reference + every spec gate) must NOT share a
# process with pyscf/qiskit-aer (block2's OpenMP runtime segfaults), so they run separately.
# Excluding test_*_spec.py by glob keeps this correct as new spec gates are added.
test:
	$(RUN) python -m pytest tests/ \
		--ignore-glob='tests/test_*_spec.py' --ignore=tests/test_dmrg_reference.py -q
	$(RUN) python -m pytest tests/test_dmrg_reference.py tests/test_*_spec.py -q

lint:
	$(RUN) ruff check .
