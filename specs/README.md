# specs/ — spec-driven development

Each new capability starts as a **spec**: a falsifiable hypothesis with concrete acceptance gates,
implemented test-first. This keeps invention honest — every claim has a check that could kill it
(the discipline this whole project was rebuilt around; see `../REFACTOR_PLAN.md`).

## The loop

1. Pick a hypothesis from [`BACKLOG.md`](BACKLOG.md) (or add one — claim + a cheap check).
2. Copy [`SPEC_TEMPLATE.md`](SPEC_TEMPLATE.md) to `SPEC_<slug>.md`; fill in goal, interface,
   **acceptance gates**, out-of-scope, honest caveats. Get it reviewed *as a doc*.
3. Write `tests/test_<slug>_spec.py` encoding the gates (initially failing).
4. Implement the **minimum** code to pass them (reuse validated primitives over new code).
5. `make gates` until green. If a gate proves unsatisfiable, **revise the spec and record why** —
   that mismatch is the finding, not a failure.

## Running the gates

```
make gates     # every tests/test_*_spec.py, each in its own process (block2 OpenMP isolation)
make test      # full suite (block2 tests isolated)
make lint      # ruff
```

## Conventions

- Specs: `specs/SPEC_<slug>.md`. Gate tests: `tests/test_<slug>_spec.py`.
- Every result is checked against a ground-truth **reference** (FCI, DMRG, experiment, analytic
  limit). No reference → not yet a spec.
- State scope limits and "reproduction vs novelty" **up front**, not after the numbers.

Worked example: [`SPEC_hchain_tdl.md`](SPEC_hchain_tdl.md) (Hₙ thermodynamic limit) — including a
G1 that was found *unsatisfiable* during implementation and honestly revised.
