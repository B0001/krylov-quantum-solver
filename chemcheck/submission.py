"""Submission validation (task 2) against architecture/interfaces/chemcheck-submission.schema.json.

Pointer-precise errors: on failure we raise ``SubmissionError`` whose message names the JSON
path to the offending field, so a vendor sees exactly what to fix.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import jsonschema
from jsonschema import Draft202012Validator
from referencing import Registry, Resource

_SCHEMA_DIR = Path(__file__).resolve().parents[1] / "architecture" / "interfaces"
_SUBMISSION_SCHEMA = "chemcheck-submission.schema.json"


class SubmissionError(ValueError):
    """A submission failed schema validation; ``.path`` is the JSON pointer to the bad field."""

    def __init__(self, message: str, *, path: str) -> None:
        super().__init__(message)
        self.path = path


@lru_cache(maxsize=1)
def _validator() -> Draft202012Validator:
    # The submission schema $refs compiler-manifest.schema.json (a sibling); register both so
    # the local $ref resolves without network access.
    resources = []
    for name in (_SUBMISSION_SCHEMA, "compiler-manifest.schema.json"):
        doc = json.loads((_SCHEMA_DIR / name).read_text())
        resources.append((name, Resource.from_contents(doc)))
    registry = Registry().with_resources(resources)
    schema = json.loads((_SCHEMA_DIR / _SUBMISSION_SCHEMA).read_text())
    return Draft202012Validator(schema, registry=registry)


def validate_submission(submission: dict[str, Any]) -> None:
    """Validate a submission dict. Raises :class:`SubmissionError` with a JSON-pointer path."""
    errors = sorted(_validator().iter_errors(submission), key=lambda e: list(e.absolute_path))
    if not errors:
        return
    err = errors[0]
    pointer = "/" + "/".join(str(p) for p in err.absolute_path) if err.absolute_path else "/"
    raise SubmissionError(f"{pointer}: {err.message}", path=pointer)


def is_mode_b(submission: dict[str, Any]) -> bool:
    """True when the submission carries device ``runs`` (Mode B), not just a spec sheet."""
    return bool(submission.get("runs"))


# Re-export for callers that want the raw jsonschema exception type too.
ValidationError = jsonschema.ValidationError
