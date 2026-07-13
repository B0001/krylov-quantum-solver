#!/usr/bin/env python3
"""
senseforge.candidates -- FoM ranking, candidates.md, and per-candidate design cards.

FoM = |S| / bracket_width at the operating point (PRD sec 5): "a huge slope with a huge error
bar is worthless." DEVIATION, recorded: every gap bracket in this repo's Nb3X8 dimer model is
EXACT (zero-width -- see senseforge/hamiltonian.py's deviation note on why certified_gaps.py's
approximate Krylov machinery does not apply to this exactly-diagonalizable system). The literal
formula divides by zero there. Rather than hide that behind a tiny epsilon, ``figure_of_merit``
falls back to FoM = |S| when width == 0, and every emitted row/card is honestly flagged
``bracket=exact`` so a reader sees this is a degenerate (fully-certain) case of the general
formula, not a silently-inflated score.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence

from senseforge.headers import ADR_0003_NOTE, markdown_header_block
from senseforge.sensitivity import Sensitivity


@dataclass(frozen=True)
class Candidate:
    halide: str
    axis: str                # "strain" | "field"
    x: float                 # operating point (eps, dimensionless; or B, tesla)
    gap: float                # best-estimate gap (meV) at x
    bracket_width: float       # width of the GAP's own bracket at x (0.0 = exact)
    slope: float               # S = d(gap)/d(axis)
    slope_lower: float
    slope_upper: float
    second_derivative: float
    fom: float
    validation_state: str = "screened"


def figure_of_merit(slope: float, bracket_width: float) -> float:
    """|S| / bracket_width, falling back to |S| at bracket_width == 0 (see module docstring)."""
    if bracket_width <= 0.0:
        return abs(slope)
    return abs(slope) / bracket_width


def build_candidates(halide: str, axis: str, xs: Sequence[float], gaps: Sequence[float],
                     bracket_widths: Sequence[float],
                     sensitivities: Sequence[Sensitivity]) -> List[Candidate]:
    """One :class:`Candidate` per interior grid point (``sensitivities`` already excludes the
    two edge points a central difference cannot cover)."""
    # sensitivities[i] corresponds to xs[i+1] (interior points only, see sensitivity.py).
    x_to_idx = {x: i for i, x in enumerate(xs)}
    out = []
    for s in sensitivities:
        i = x_to_idx[s.x]
        fom = figure_of_merit(s.slope, bracket_widths[i])
        out.append(Candidate(
            halide=halide, axis=axis, x=s.x, gap=gaps[i], bracket_width=bracket_widths[i],
            slope=s.slope, slope_lower=s.slope_lower, slope_upper=s.slope_upper,
            second_derivative=s.second_derivative,
            fom=fom,
        ))
    return out


def rank_candidates(candidates: Sequence[Candidate]) -> List[Candidate]:
    """Descending by FoM -- the ranked screening order the PRD's candidates.md wants."""
    return sorted(candidates, key=lambda c: c.fom, reverse=True)


#: Verdicts from :func:`ranking_verdict`.
RANKING_DEGENERATE = "degenerate"
RANKING_MONOTONE = "monotone"
RANKING_INTERIOR = "interior"


def ranking_verdict(candidates: Sequence[Candidate], *, rel_tol: float = 1e-9) -> str:
    """Does the FoM ranking actually DISCRIMINATE between operating points, or not?

    THE FINDING (SPEC_senseforge.md, gate G9). On the Nb3X8 dimer -- the only model this repo can
    reach -- the ranking never identifies a genuine operating point, for two different reasons:

    * ``"degenerate"`` -- every FoM is equal. The ranking is then pure sort order and carries NO
      information; "rank 1" is not a recommendation. This is the FIELD axis: the Zeeman response
      is EXACTLY linear (gate G2), so d(gap)/dB is the same constant at every B, and all 19
      operating points tie at |S| = 0.1158 meV/T. The pre-fix design_card_1.md nonetheless
      published "+2 T" as the top operating point -- sort noise presented as a design choice.
    * ``"monotone"`` -- |S| is monotone across the swept grid, so the top-ranked point is the
      WINDOW EDGE, not an optimum: widen the window and the "best" point obediently moves with it
      (verified: [-2%,+2%] ranks +1.75%; [-5%,+5%] ranks +4.75%). This is the STRAIN axis, where
      J(t) is smooth and monotone over the accessible range (|S| spans only 123.3 -> 126.7, a
      2.8% spread).
    * ``"interior"`` -- the top-ranked point is a genuine interior optimum. NOT OBSERVED on this
      model; it is what a discriminating screener would have to produce to be worth anything.

    Callers must surface this: :func:`render_candidates_md` and :func:`render_design_card` stamp
    the verdict on every artifact, so a reader cannot mistake a tie or a window edge for a result.
    """
    if len(candidates) < 2:
        return RANKING_INTERIOR
    foms = [c.fom for c in candidates]
    spread = max(foms) - min(foms)
    if spread <= rel_tol * max(abs(f) for f in foms):
        return RANKING_DEGENERATE

    ranked = rank_candidates(candidates)
    xs = sorted(c.x for c in candidates)
    if ranked[0].x in (xs[0], xs[-1]):
        return RANKING_MONOTONE
    return RANKING_INTERIOR


def _verdict_banner(candidates: Sequence[Candidate]) -> str:
    """The honest one-liner every artifact carries about what the ranking does (and does not) mean."""
    verdict = ranking_verdict(candidates)
    if verdict == RANKING_DEGENERATE:
        return ("> **NO DISCRIMINATION (degenerate ranking):** every operating point in this sweep "
                "has the SAME figure of merit, so the rank order below is sort order, not a "
                "recommendation -- no point is better than any other. (The response is exactly "
                "linear on this axis; see SPEC_senseforge.md.)")
    if verdict == RANKING_MONOTONE:
        return ("> **NO INTERIOR OPTIMUM (monotone ranking):** |S| increases monotonically across "
                "the swept window, so rank 1 is simply the WINDOW EDGE -- widen the window and it "
                "moves with it. This is not a discovered operating point. (See SPEC_senseforge.md.)")
    return "> Ranking discriminates: rank 1 is an interior optimum of |S| within the swept window."


def render_candidates_md(candidates: Sequence[Candidate], config_header: dict, cluster: str) -> str:
    """The ranked table (PRD sec 6: halide, operating point, S, FoM, bracket width, flags)."""
    ranked = rank_candidates(candidates)
    lines = [markdown_header_block(config_header, cluster)]
    lines.append("# SenseForge candidates\n")
    lines.append(_verdict_banner(candidates) + "\n")
    lines.append("| rank | halide | axis | operating point | gap (meV) | S | bracket | FoM | flags |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for rank, c in enumerate(ranked, start=1):
        bracket_note = "exact" if c.bracket_width <= 0.0 else f"{c.bracket_width:.4g}"
        unit = "" if c.axis == "strain" else " T"
        lines.append(
            f"| {rank} | Nb3{c.halide}8 | {c.axis} | {c.x:+.4g}{unit} | {c.gap:.4f} | "
            f"{c.slope:.4g} | {bracket_note} | {c.fom:.4g} | {c.validation_state} |"
        )
    return "\n".join(lines) + "\n"


def render_design_card(rank: int, c: Candidate, config_header: dict, cluster: str,
                       siblings: Sequence[Candidate] = ()) -> str:
    """One design card (PRD sec 6): operating point, sensitivity, FoM, flags, validation_state.

    ``siblings`` is the full candidate set this card was ranked within -- required to state
    honestly whether "rank 1" means anything at all (see :func:`ranking_verdict`). A card is a
    RECOMMENDATION; it must not imply one was made when the ranking was a tie or a window edge.
    """
    unit = "" if c.axis == "strain" else " T"
    bracket_note = "exact (zero-width -- see senseforge/hamiltonian.py)" if c.bracket_width <= 0.0 \
        else f"[{c.slope_lower:.4g}, {c.slope_upper:.4g}]"
    lines = [
        markdown_header_block(config_header, cluster),
        f"# Design card #{rank}: Nb3{c.halide}8 {c.axis} sensor\n",
        _verdict_banner(siblings) + "\n" if siblings else "",
        f"- **Operating point:** {c.axis} = {c.x:+.4g}{unit}",
        f"- **Gap at operating point:** {c.gap:.4f} meV",
        f"- **Sensitivity S = d(gap)/d({c.axis}):** {c.slope:.4g} meV/{'unit strain' if c.axis == 'strain' else 'T'}",
        f"- **Sensitivity bracket:** {bracket_note}",
        f"- **Second derivative (plateau check):** {c.second_derivative:.4g}"
        + (" (flat -- stable operating point)" if abs(c.second_derivative) < abs(c.slope) * 0.1
           else " (curved -- a knife-edge point, see PRD sec 5)"),
        f"- **Figure of merit:** {c.fom:.4g}",
        f"- **Validation state:** {c.validation_state}",
        "",
        f"> {ADR_0003_NOTE}",
    ]
    return "\n".join(lines) + "\n"


def write_candidate_report(candidates: Sequence[Candidate], config_header: dict, cluster: str,
                           output_dir: str) -> Path:
    """Writes candidates.md + design_card_{1,2,3}.md for the top-3; returns candidates.md path."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    ranked = rank_candidates(candidates)

    md_path = out / "candidates.md"
    md_path.write_text(render_candidates_md(candidates, config_header, cluster))

    for rank, c in enumerate(ranked[:3], start=1):
        card_path = out / f"design_card_{rank}.md"
        card_path.write_text(render_design_card(rank, c, config_header, cluster,
                                                siblings=candidates))

    return md_path
