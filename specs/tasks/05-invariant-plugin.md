# Task Breakdown 5 — #20 Invariant-Test Framework (pytest plugin)
Goal: extract the floor-check pattern into `pytest-invariants`; dogfood on the solver's own CI. Cheapest real invention in the portfolio.

1. **API design doc (1 page)** — decorators: `@invariant.lower_bound(fn_or_value)`, `@invariant.contains(reference, key=...)`, `@invariant.monotone(param, direction)`, `@invariant.conserved(quantity, tol)`. Decide: decorate test functions (v1, simple) vs instrument production functions (v2, defer).
   ✓ Doc reviewed against the solver's actual checks — every existing hand-rolled check expressible. (S)
2. **Plugin skeleton** — cookiecutter pytest plugin; entry point registers an `invariants` marker; `pytest --invariant-report` flag stub.
   ✓ Installs; `pytest -p pytest_invariants` loads clean. (S)
3. **Implement `lower_bound` + `contains`** — the two the solver needs. Failure message format: expected relation, observed values, provenance of the bound.
   ✓ Red/green demo tests; failure message readable without reading the plugin source. (M)
4. **Dogfood refactor** — replace the solver's hand-rolled containment/floor assertions in the golden suite with plugin decorators. Zero behavior change intended.
   ✓ CI green before and after; diff is deletions + decorators. This diff IS the README's centerpiece. (M)
5. **Invariant-coverage report** — `--invariant-report` emits which invariants ran, over which tests, pass counts.
   ✓ Report for the solver repo generated in CI artifacts. (S)
6. **Write README from the dogfood diff** — motivation section = the old codebase's −hundreds-of-Hartree story (anonymized or not, your call).
   ✓ Quickstart reproducible; PyPI-ready metadata. (S)
7. **Publish `v0.1`** — PyPI + repo; solver depends on it from the next release.
   ✓ `pip install pytest-invariants` works cold. (S)
