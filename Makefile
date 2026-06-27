# Spec-driven development targets. See specs/README.md.
#
# Override the env runner if you don't use conda (e.g. `make gates RUN=`):
#   ENV ?= chem   ->  commands run via `conda run -n chem`
ENV ?= chem
RUN ?= conda run -n $(ENV)

.PHONY: gates test lint help

help:
	@echo "make gates  - run every tests/test_*_spec.py, each in its OWN process"
	@echo "make test   - full test suite (block2/DMRG tests isolated)"
	@echo "make lint   - ruff check"

# Acceptance gates. Each spec's gate file runs in a separate process so block2's OpenMP runtime
# never loads into an interpreter that already imported pyscf/qiskit-aer (which segfaults).
gates:
	@set -e; \
	shopt -s nullglob 2>/dev/null || true; \
	files=$$(ls tests/test_*_spec.py 2>/dev/null); \
	if [ -z "$$files" ]; then echo "no spec gate tests (tests/test_*_spec.py)"; exit 0; fi; \
	for f in $$files; do \
		echo "=== gates: $$f ==="; \
		$(RUN) python -m pytest $$f -q || exit 1; \
	done; \
	echo "all spec gates passed"

# Full suite: block2/DMRG tests in their own process (see run_in_chem.sh).
test:
	$(RUN) python -m pytest tests/ \
		--ignore=tests/test_dmrg_reference.py --ignore=tests/test_hchain_tdl_spec.py -q
	$(RUN) python -m pytest tests/test_dmrg_reference.py tests/test_hchain_tdl_spec.py -q

lint:
	$(RUN) ruff check .
