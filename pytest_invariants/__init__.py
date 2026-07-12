"""pytest-invariants — reference-checked test invariants as decorators.

Extracted from this repo's hand-rolled variational-floor / bracket-contains-reference checks
(``certchem.floor_guard`` and the golden-suite ``lower <= FCI <= upper`` assertions). A test
returns its observation and declares what must hold about it:

    from pytest_invariants import invariant

    @invariant.lower_bound(FCI_H2, provenance="FCI/STO-3G", key=lambda b: b.best_estimate)
    @invariant.contains(FCI_H2, provenance="FCI/STO-3G", key=lambda b: (b.lower, b.upper))
    def test_h2_bracket():
        return certified_energy("H 0 0 0; H 0 0 0.735", "sto-3g", (2, 2)).bracket

The pytest plugin (:mod:`pytest_invariants.plugin`) registers the ``invariants`` marker and a
``--invariant-report`` coverage report.
"""

from __future__ import annotations

from .api import InvariantSpec, InvariantViolation, invariant

# Re-export the plugin hooks into the package namespace so THIS module *is* the pytest plugin:
# both `pytest -p pytest_invariants` and the pytest11 entry point resolve to this one module
# object, so pluggy registers the hooks exactly once (no double-report when both are in play).
from .plugin import (  # noqa: F401  (re-exported for pluggy hook discovery)
    pytest_addoption,
    pytest_collection_modifyitems,
    pytest_configure,
    pytest_runtest_makereport,
    pytest_terminal_summary,
)

__all__ = ["invariant", "InvariantViolation", "InvariantSpec"]

__version__ = "0.1.0"
