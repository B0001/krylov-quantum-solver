# pytest-invariants

Reference-checked test **invariants** as decorators. A test declares *what must physically hold*
about its result — that an energy estimate never sinks below a variational floor, that a certified
bracket brackets a known reference — and a violation fails the test with a message that explains
itself: the expected relation, the observed values, and the **provenance** of the bound.

## Motivation — the "-hundreds of Hartree" story

This plugin was extracted from a hybrid quantum-chemistry solver that had been **rebuilt after a
scientific audit**. The audit's headline finding: with noise injected, the old physics core
returned ground-state energies of **-5, -39, -191, -813 Hartree** for a molecule whose true minimum
is **-1.85 Ha**. A variational / Krylov method *can never* go below the true ground state — those
numbers were the loud signature of a broken estimator, yet every unit test passed, because the
tests asserted the broken code's *own output* and never compared against a reference.

The fix was a culture rule: **every result has a reference (FCI, DMRG, experiment, analytic limit)
and a test that could kill it.** The single most valuable such test is the *variational floor*:
`estimate >= floor`. That one check, hand-rolled as `certchem.floor_guard`, is exactly what
`@invariant.lower_bound` packages — so any project can catch a below-the-floor regression the
instant it appears, with a failure message a reviewer can read at a glance.

## Install

```bash
pip install -e .          # from this repo; registers the pytest11 entry point
# or, once published:
pip install pytest-invariants
```

The plugin auto-loads via its `pytest11` entry point. To load explicitly:

```bash
pytest -p pytest_invariants
```

## Quickstart (reproducible)

```python
# test_energy.py
from pytest_invariants import invariant

FCI_H2 = -1.137270  # reference ground-state energy, Ha (provenance below travels into failures)

@invariant.lower_bound(FCI_H2, provenance="FCI/STO-3G")
def test_estimate_respects_the_floor():
    estimate = -1.10                 # your solver's estimate
    return estimate                  # the test RETURNS the observation to check

@invariant.contains(FCI_H2, provenance="FCI/STO-3G")
def test_bracket_brackets_fci():
    lower, upper = -1.20, -1.05      # your certified bracket
    return (lower, upper)
```

```bash
pytest test_energy.py --invariant-report
```

A decorated test **returns the observation** (a scalar, a `(lower, upper)` tuple, or an object plus
a `key=` extractor). Both bound and reference may be a literal or a zero-arg callable (e.g. an FCI
lookup computed at test time).

### What a violation looks like

A sub-floor estimate (`return -813.0`) fails with:

```
Invariant `lower_bound` violated in test_estimate_respects_the_floor:
  expected:  observed >= bound
  observed:  -813 Ha
  bound:     -1.13727 Ha  (provenance: FCI/STO-3G)
  margin:    observed - bound = -811.86273 Ha  (must be >= 0)
```

A bracket that misses the reference (`return (-2.00, -1.90)`) fails with:

```
Invariant `contains` violated in test_bracket_brackets_fci:
  expected:  lower <= reference <= upper
  reference: -1.13727 Ha  (provenance: FCI/STO-3G)
  interval:  [-2 Ha, -1.9 Ha]
  reference is ABOVE the interval by 0.76273 Ha
```

## Invariant-coverage report

`--invariant-report` prints which invariants ran, over which tests, with pass counts:

```
========================== invariant coverage report ===========================
invariant         checks  passed   tests
----------------------------------------
contains               1       1       1
lower_bound            1       1       1
----------------------------------------
TOTAL                  2       2

per-test detail:
  [PASS] lower_bound    test_energy.py::test_estimate_respects_the_floor
  [PASS] contains       test_energy.py::test_bracket_brackets_fci
```

Decorated tests are auto-tagged with the `invariants` marker, so `pytest -m invariants` runs only
the reference-checked ones.

## API

| decorator | asserts | typical use |
|-----------|---------|-------------|
| `@invariant.lower_bound(fn_or_value, *, key=None, tol=0.0, provenance="…", units="Ha")` | `observed >= bound` | variational floor |
| `@invariant.contains(reference, *, key=None, provenance="…", units="Ha")` | `lower <= reference <= upper` | certified bracket contains FCI |

`key` extracts the checked value from the returned observation (identity by default). Stacking both
on one test checks both against the same observation. See [`DESIGN.md`](DESIGN.md) for the full
design and the audit showing every hand-rolled check in the source solver is expressible.

## Status

v0.1 — the two decorators the solver's golden suite needs. `upper_bound`, `monotone`, and
`conserved` are on the roadmap (see `DESIGN.md`); they are intentionally not shipped until a real
test exercises them.
