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

# Full suite. The block2/DMRG tests (test_dmrg_reference + every spec gate) must NOT share a
# process with pyscf/qiskit-aer (block2's OpenMP runtime segfaults), so they run separately.
# Excluding test_*_spec.py by glob keeps this correct as new spec gates are added.
test:
	$(RUN) python -m pytest tests/ \
		--ignore-glob='tests/test_*_spec.py' --ignore=tests/test_dmrg_reference.py -q
	$(RUN) python -m pytest tests/test_dmrg_reference.py tests/test_*_spec.py -q

lint:
	$(RUN) ruff check .
