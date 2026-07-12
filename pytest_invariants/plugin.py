"""The pytest plugin: registers the ``invariants`` marker, the ``--invariant-report`` flag,
and the invariant-coverage report.

Loaded either explicitly (``pytest -p pytest_invariants``) or automatically via the
``pytest11`` entry point declared in ``pyproject.toml``. The decorators in :mod:`.api` enforce
their invariants on their own (inside the wrapper), so a test's pass/fail does not depend on
this plugin being active; the plugin adds discovery (the marker), the coverage report, and
``-m invariants`` selection.
"""

from __future__ import annotations

from collections import Counter
from typing import Any

import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    """Register ``--invariant-report`` (which invariants ran, over which tests, pass counts)."""
    group = parser.getgroup("invariants", "pytest-invariants: reference-checked test invariants")
    group.addoption(
        "--invariant-report",
        action="store_true",
        default=False,
        help="After the run, print which invariants ran, over which tests, with pass counts.",
    )


def pytest_configure(config: pytest.Config) -> None:
    """Register the ``invariants`` marker and attach the run-collector to the session config."""
    config.addinivalue_line(
        "markers",
        "invariants: test carries one or more @invariant.* reference checks (auto-applied).",
    )
    # list of (nodeid, invariant_name, passed) rows, one per invariant per decorated test.
    config._invariant_records = []  # type: ignore[attr-defined]


def _invariants_of(item: pytest.Item) -> list:
    """Return the InvariantSpec list a test carries, or [] if it is not invariant-decorated."""
    func = getattr(item, "obj", None)
    return list(getattr(func, "__invariants__", []) or [])


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Auto-tag every invariant-decorated test with the ``invariants`` marker so ``-m
    invariants`` selects them without the author adding the marker by hand."""
    for item in items:
        if _invariants_of(item):
            item.add_marker("invariants")


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item: pytest.Item, call: pytest.CallInfo) -> Any:
    """Record, on the 'call' phase, each invariant the test carried and whether it passed."""
    outcome = yield
    report = outcome.get_result()
    if report.when != "call":
        return
    specs = _invariants_of(item)
    if not specs:
        return
    records = getattr(item.config, "_invariant_records", None)
    if records is None:
        return
    passed = report.outcome == "passed"
    for spec in specs:
        records.append((report.nodeid, spec.name, passed))


def pytest_terminal_summary(
    terminalreporter: Any, exitstatus: int, config: pytest.Config
) -> None:
    """Print the invariant-coverage report when ``--invariant-report`` was passed."""
    if not config.getoption("--invariant-report"):
        return
    records = getattr(config, "_invariant_records", [])
    tw = terminalreporter
    tw.write_sep("=", "invariant coverage report")
    if not records:
        tw.write_line("No @invariant.* decorated tests ran.")
        return

    checks: Counter = Counter()
    passed_checks: Counter = Counter()
    tests_by_name: dict[str, set] = {}
    for nodeid, name, passed in records:
        checks[name] += 1
        if passed:
            passed_checks[name] += 1
        tests_by_name.setdefault(name, set()).add(nodeid)

    tw.write_line(f"{'invariant':<16}{'checks':>8}{'passed':>8}{'tests':>8}")
    tw.write_line("-" * 40)
    for name in sorted(checks):
        tw.write_line(
            f"{name:<16}{checks[name]:>8}{passed_checks[name]:>8}"
            f"{len(tests_by_name[name]):>8}"
        )
    tw.write_line("-" * 40)
    tw.write_line(f"{'TOTAL':<16}{sum(checks.values()):>8}{sum(passed_checks.values()):>8}")
    tw.write_line("")
    tw.write_line("per-test detail:")
    for nodeid, name, passed in records:
        mark = "PASS" if passed else "FAIL"
        tw.write_line(f"  [{mark}] {name:<14} {nodeid}")
