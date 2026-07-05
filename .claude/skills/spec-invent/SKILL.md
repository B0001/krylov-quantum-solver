---
name: spec-invent
description: Add a new falsifiable, test-gated capability to this hybrid quantum-classical chemistry repo (or extend the ODMD suite / study the Nb3X8 clusters) using its spec-driven "falsifiable honesty" loop. Use when the user asks to invent, add, or extend a method/study here, when they say "spec-invent", or when work in this repo would produce a numeric claim that needs a reference and a test that could kill it. Produces a SPEC_<slug>.md, a failing tests/test_<slug>_spec.py gate, minimal code reusing validated primitives, a recorded boundary/finding, and an honest commit.
---

# spec-invent — the falsifiable-honesty invention loop

This repo was rebuilt after an audit found its original core entirely broken (fabricated
"Krylov" bases, asymmetric noise on H, a hashing "integrity" check that certified collapsed
runs). The culture that replaced it: **every claimed number has a ground-truth reference (FCI,
DMRG, experiment, or an analytic limit) and a test that could kill it; when a gate proves
unsatisfiable, you revise the spec and record why — the boundary IS the finding, not a failure.**
This skill runs a new capability through that loop. Do not shortcut it; the shortcut is how the
repo broke the first time.

## Non-negotiables (read `CLAUDE.md` first)

- **Every shell command uses `conda run -n chem`.** (Fault-tolerant/openfermion work uses the
  separate `chem-ft` env — see the memory notes.) Never assume the env is active.
- **block2/DMRG process isolation:** block2 segfaults if it loads into a process that already
  imported pyscf or qiskit-aer. `make gates` runs each `tests/test_*_spec.py` in its own process
  and keeps DMRG tests separate — preserve this; the `test_*_spec.py` glob keeps it correct.
  DMRG-only test files must not import qiskit.
- **Do not touch the quarantined files** (`orchestrate_hybrid_pipeline.py`, `quantum_sampler.py`)
  — regression fixtures only, never import/fix/build on them.
- **scipy is pinned `>=1.8,<1.16`** (1.16+ crashes the qiskit-nature import chain); DataFrames are
  **polars**, not pandas.
- **Reuse validated primitives over new code.** The public API is in
  `hybrid_quantum_solver/__init__.py`; the ODMD suite (`odmd.py`, `odmd_spectral.py`,
  `odmd_optical.py`, `odmd_spin.py`, `device_odmd.py`, `trotter_odmd.py`, `visibility_law.py`,
  `odmd_uq.py`, `temple_bounds.py`, `msd.py`, …) composes it. `docs/ODMD_SUITE.md` is the map.

## The loop (do these in order)

1. **Probe first — get real numbers before freezing any gate.** Write a throwaway script in the
   scratchpad (NOT the repo) that measures the thing you're about to claim, against the exact
   reference (dense diagonalization / `fixed_filling_energy` / `mh.ground_state_energy()` /
   analytic closed form). Every gate threshold in this repo is a *measured* number with margin,
   never a guess. If the probe kills the idea or moves the boundary, that is a success — record it
   and adjust. Long DMRG probes: run detached in the user's terminal, not the session.

2. **Write `specs/SPEC_<slug>.md`** from `specs/SPEC_TEMPLATE.md`. Mandatory sections: the single
   falsifiable claim; **honest framing up front** — what you CAN claim if gates pass and what you
   CANNOT (reproduction vs novelty, "no quantum advantage at this scale", finite-cluster/
   minimal-basis, exact-statevector vs shot noise — state these *before* the numbers, not after);
   the reference every result is checked against; the acceptance gates G1..Gn as concrete numeric
   pass conditions with the measured value in a comment; which gate is the **definition of done**;
   out-of-scope; caveats/risks. Name it `SPEC_<slug>.md`, slug in `snake_case`.

3. **Write `tests/test_<slug>_spec.py` encoding the gates — and confirm it is RED** (the module
   doesn't exist yet). Medians over noise seeds for stochastic gates; deterministic seeds so gates
   are reproducible. Keep a gate that pins the **boundary/finding** (where the method stops
   working), not just the happy path — that is what makes it falsifiable.

4. **Implement the minimum code to pass**, composing validated primitives. New top-level method
   modules sit at repo root beside `odmd.py`; genuinely core additions go in
   `hybrid_quantum_solver/` and get exported. Match the surrounding style: ~100-col lines, honest
   docstrings that state scope and the found boundary.

5. **`conda run -n chem python -m pytest tests/test_<slug>_spec.py`** to green, then
   `conda run -n chem ruff check <new files>`. If a gate proves unsatisfiable, **revise the gate
   in the spec with a note on what reality showed** (see `SPEC_hchain_tdl.md` G1 and the `[-]`
   killed entries in `specs/BACKLOG.md` for worked examples) — do not weaken it silently.

6. **Record the finding in `specs/BACKLOG.md`** as a `[x]` entry: the headline number, the
   recorded boundary/finding, the gate location, the spec link, and the honest-scope parenthetical.
   Update `docs/ODMD_SUITE.md` if it extends the suite.

7. **Commit** (only when the user asks, or per their standing instruction). Message leads with the
   claim, lists the gated results with numbers, and states the boundary and honest scope — the same
   discipline as the spec. If you found and fixed a latent defect while probing (it happens — this
   session found three), make it a regression gate and say so in the message.

## What "done" means

`make gates` green including the new file; ruff clean on new files; the spec's definition-of-done
gate passing; the boundary recorded. A capability without a test that could have killed it is not
done — it is exactly what this repo exists to prevent.

## Anti-patterns (these are how the original core broke)

- A number with no reference, or a gate that cannot fail.
- Widening a tolerance to pass instead of recording why the claim is narrower than hoped.
- Presenting a reproduction as novelty, or a simulated/finite-cluster result as a materials claim.
- An error bar that quantifies variance and is quoted as if it caught bias (see `SPEC_odmd_uq` G4).
- "Integrity" theatre (hashes, symmetry stamps) in place of a real reference check.
