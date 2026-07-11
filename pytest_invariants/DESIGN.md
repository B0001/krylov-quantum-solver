# pytest-invariants — API design note (v1)

**Goal.** Extract this repo's hand-rolled "reference-checked" test assertions — the variational
*floor* check and the *bracket-contains-a-reference* check — into reusable decorators, so a test
declares *what must hold* about its result and the failure explains itself. v1 decorates **test
functions** (simple, no production instrumentation). Instrumenting production functions is v2 and
deliberately deferred.

## Model

A decorated test **returns an observation** (a scalar, a `(lower, upper)` tuple, or any object).
Each stacked `@invariant.<name>(...)` pulls what it needs from that observation via an optional
`key` and asserts one relation, raising `InvariantViolation` (an `AssertionError` subclass) on
failure. Decorators compose onto a single wrapper, so the test body runs once and the wrapper
returns `None` to pytest (no `PytestReturnNotNoneWarning`).

```python
@invariant.lower_bound(FCI_H2, provenance="FCI/STO-3G", key=lambda b: b.best_estimate_hartree)
@invariant.contains(FCI_H2, provenance="FCI/STO-3G", key=lambda b: (b.lower_hartree, b.upper_hartree))
def test_h2():
    return certified_energy("H 0 0 0; H 0 0 0.735", "sto-3g", (2, 2)).bracket
```

## v1 decorators

| decorator | asserts | replaces |
|-----------|---------|----------|
| `lower_bound(fn_or_value, *, key, tol, provenance)` | `observed >= bound` | `certchem.floor_guard` |
| `contains(reference, *, key, provenance)` | `lower <= reference <= upper` | golden bracket check |

`fn_or_value` / `reference` may be a literal **or** a zero-arg callable (resolved at test time,
e.g. an FCI lookup). `provenance` is a human string naming where the bound came from; it is echoed
verbatim in the failure message, so a failure reads without opening any source.

## Expressibility audit — every existing hand-rolled check maps

Checked against the solver's actual assertions (`certchem/core.py`, `certchem/scorecard` reuse in
`chemcheck`, and `tests/test_certchem_core_spec.py`):

| Hand-rolled check (where) | Expressed as |
|---|---|
| `estimate >= floor` — `certchem.floor_guard` (the `-hundreds of Hartree` guard) | `lower_bound(floor, provenance="variational floor")` |
| `floor_guard(-999, -1.1)` unit raises — `test_floor_guard_unit_raises_below_floor` | `lower_bound` red case |
| `b.lower <= fci <= b.upper` — `test_golden_bracket_contains_fci_within_tol` | `contains(fci, key=lambda b: (b.lower_hartree, b.upper_hartree))` |
| `abs(estimate - fci) <= tol` — same golden test (accuracy band) | `contains(estimate, reference-window)`: `contains(fci, key=lambda e: (e-tol, e+tol))` |
| `chemcheck.scorecard` floor detector (reuses `floor_guard`) | `lower_bound` (same as above) |
| scorecard accuracy bands (`err_mha <= PASS/MARGINAL`) | `contains` over the ± tolerance window |
| `b.width >= 0.0` — golden test | `lower_bound(0.0, key=lambda b: b.width, units="Ha")` |

**Honest boundaries.**
- `result.certificate.floor_check == "pass"` is a *status equality*, not a numeric invariant — it
  is the recorded *consequence* of the `lower_bound` check, so it is out of scope by design.
- The scorecard's one-sided `headroom <= 1.0` is an *upper* bound. v1 ships only the two decorators
  the solver's golden suite needs; a mirror `upper_bound` is a one-line addition (or expressible
  today as `lower_bound(-limit, key=lambda x: -x)`), recorded here rather than built speculatively.
- `monotone(param, direction)` and `conserved(quantity, tol)` (from the task's wishlist) have **no**
  current hand-rolled counterpart in the golden suite, so they are v2 — not implemented, to avoid
  untested surface area (falsifiable honesty: no capability without a check that exercises it).

## Plugin surface (v1)

Registers the `invariants` marker (auto-applied to decorated tests, so `-m invariants` selects
them) and `--invariant-report` (which invariants ran, over which tests, with pass counts). Loads
via `pytest -p pytest_invariants` or the `pytest11` entry point. Because the decorators enforce
their own checks, a test's pass/fail never depends on the plugin being active — the plugin only
adds discovery and reporting.
