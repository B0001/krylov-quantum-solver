#!/usr/bin/env python3
"""
senseforge.headers -- the ADR-0003 mandatory scope note on every SenseForge artifact.

architecture/adr/0003-finite-cluster-scope-for-materials.md, decision 1: "Every artifact derived
from the CIF path carries an automatic header: 'cluster-model prediction; not validated for the
periodic solid.'" Non-optional, no flag to disable (task 4). SenseForge's Nb3X8 modeling is not
literally CIF-derived (see senseforge/hamiltonian.py's deviation note) but is exactly the finite
isolated-dimer cluster model ADR-0003 was written to scope -- the header applies with equal force.
"""
from __future__ import annotations

ADR_0003_NOTE = "cluster-model prediction; not validated for the periodic solid."


def csv_header_lines(config_header: dict, cluster: str) -> list:
    """Comment-prefixed (``#``) lines to prepend to a CSV artifact."""
    return [
        f"# {ADR_0003_NOTE}",
        f"# cluster: {cluster}",
        f"# config: {_format_kv(config_header)}",
    ]


def markdown_header_block(config_header: dict, cluster: str) -> str:
    """A markdown blockquote to open a ``.md`` artifact (candidates.md, design cards) with."""
    lines = [
        f"> **{ADR_0003_NOTE}**",
        f"> Cluster: {cluster}",
        f"> Config: `{_format_kv(config_header)}`",
        "",
    ]
    return "\n".join(lines)


def _format_kv(d: dict) -> str:
    return ", ".join(f"{k}={v}" for k, v in sorted(d.items()))


def has_adr0003_note(text: str) -> bool:
    """True if ``text`` contains the mandatory scope note verbatim -- what the gate greps for."""
    return ADR_0003_NOTE in text
