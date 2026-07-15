# SPEC: validate_and_cost — one threshold quietly governs two independent regime boundaries

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

---

## 1. Goal

`validate_and_cost.py` composes three already-gated stages (taper, cross-check, FT cost) into "one
call that takes an active space from integrals to a costed, cross-validated verdict." Never gated
itself. This spec checks the composition, not the stages (each already has its own spec): does the
full pipeline actually work end-to-end on a small system, does the regime-transition logic correctly
degrade at a larger one, and — the sharper question — does `qubit_dense_max_orb`, ONE parameter
forwarded through this module, quietly control TWO conceptually independent regime boundaries at
once (this module's own taper-skip gate, and `cross_check`'s internal ADAPT-skip gate), with no
documentation saying so? False if the two boundaries diverge (one fires without the other), or if
the pipeline crashes rather than degrading gracefully when a stage is out of regime.

## 2. Background and honest framing

- `validate_and_cost.py` reuses `taper_hamiltonian`, `cross_check`, and (conditionally)
  `accuracy_gate` unmodified — no new physics, only falsifiers around the composition seam.
- **What you can claim if the gates pass:** on a small CAS (H4 CAS(4,4)), the full three-stage
  pipeline runs (taper reduces qubits, cross-check's four methods agree, FT cost gracefully reports
  `None` in an environment without `openfermion` — the standard `chem` env this repo's `make gates`
  runs in); on a larger CAS (N2, 8 orbitals, exceeding the default `qubit_dense_max_orb=7`), taper
  and ADAPT (inside cross-check) BOTH skip simultaneously, driven by the SAME forwarded threshold —
  a coupling the source code creates but nothing documents; and despite two of five
  sub-components being out of regime, the surviving methods (CASCI, Krylov, SQD) still reach a
  valid, non-vacuous agreement verdict — the pipeline degrades gracefully, not catastrophically.
- **What you cannot claim:** anything about stage 3's actual numbers (never exercised here — this
  repo's standard `chem` env has no `openfermion`, and `chem-ft` isn't part of the default `make
  gates` flow; the graceful-skip path is what's gated, not the FT cost formula itself, which
  `SPEC_precision_cost.md`/`SPEC_scdf_lambda.md` already cover elsewhere); that the shared-threshold
  coupling (G2) is a bug — it may be an intentional simplification (both gates ask "is the qubit
  space too big for a dense operation"), but it is currently unexamined and undocumented, which is
  the actual finding.
- **Reference:** the exact CASCI energy at each tested active space (dense diagonalization, the same
  reference `cross_check` itself already validates against) for the agreement claims; structural
  inspection of the returned report dict (keys present/absent) for the regime-boundary claims.

## 3. Approach

Reuse `validate_and_cost`, `print_report` unmodified, driven entirely through the public API. Small
system: H4 CAS(4,4) (well under the default `qubit_dense_max_orb=7`). Large system: N2 CASCI(8
orbitals, 10 electrons) — 8 orbitals exceeds the default threshold, and `(na, nb) = (5, 5)` keeps
the CI dimension (3136) comfortably under `cross_check`'s own `krylov_max_dim`/`fci_max_dim` caps, so
CASCI/Krylov/SQD all still run while only taper and ADAPT (the two `qubit_dense_max_orb`-gated
stages) drop out. Both `print_report(...)` calls exercised directly (not just the returned dict) to
catch any conditional-branch crash the structural checks alone might miss.

The threshold-coupling gate (G2) is proven by LOWERING `qubit_dense_max_orb` on the already-small H4
system (norb=4), not by raising it on the large N2 system. Raising it on N2 was tried while probing
this spec — it lets taper actually attempt `pauli_decompose` at 16 qubits (norb=8), which is
exponential in qubit count (`SPEC_taper_spectrum.md`'s own R1 measured ~1000s at 8 qubits; 16 is many
orders of magnitude worse) and OOM-killed the process. Proving the same coupling in the safe
direction — shrinking an already-tractable system's cap below its own `norb` — reaches the identical
conclusion (both gates move together) without ever running taper outside its known-tractable range.

## 4. Public interface

No new library code — this spec adds only test-file assertions around `validate_and_cost.py`'s
existing `validate_and_cost`/`print_report`, reused unchanged.

## 5. Acceptance criteria (validation gates)

- **G1 — small CAS: the full three-stage pipeline composes correctly.** On H4 CAS(4,4): taper is
  NOT skipped (`"skipped" not in report["taper"]`) and reduces qubit count
  (`n_qubits_tapered < n_qubits_original`); `report["cross_check"]["agree"]` is `True`; `report["ft_cost"]
  is None` (the standard `chem` env has no `openfermion`, confirmed via `validate_and_cost._HAVE_FT
  is False`). *Measured: taper 8 -> 5 qubits; cross-check max deviation 0.0047 mHa; ft_cost None.*
- **G2 — THE FINDING: one threshold, two regime boundaries, no documentation.** On H4 CAS(4,4)
  (norb=4): at the default `qubit_dense_max_orb=7`, taper runs and ADAPT is attempted (both "on");
  at an explicitly LOWERED `qubit_dense_max_orb=3` on the SAME system, taper IS skipped
  (`"skipped" in report["taper"]`) AND `"ADAPT"` appears in `report["cross_check"]["skipped"]` —
  both driven by the identical forwarded threshold, flipping off together, never independently.
  (Proven in the safe direction — see §3 — not by raising the cap on a large system.)
- **G3 — graceful degradation, not a crash: the surviving methods still reach a valid verdict.**
  On the same N2 case: `report["cross_check"]["reference"] is not None`,
  `report["cross_check"]["agree"] is True`, and `{"CASCI", "Krylov", "SQD"} <=
  set(report["cross_check"]["results"])` — three of four methods still ran and agreed despite taper
  and ADAPT both being out of regime. *Measured: max deviation 0.0012 mHa.*
- **G4 — `print_report` doesn't crash across either regime.** `print_report(...)` called directly
  (not just the structural dict checks) on both the H4 (full pipeline) and N2 (taper+ADAPT skipped,
  ft_cost None) reports, with no exception.

> Definition of done: **G2**. G1/G3/G4 build the surrounding confidence that the pipeline works and
> degrades safely; G2 is the actual coupling nobody had checked.

## 6. Implementation plan (test-first)

1. Write `tests/test_validate_and_cost_composition_spec.py` encoding G1-G4 (RED in the sense these
   checks are new, even though `validate_and_cost.py`'s functions are not).
2. No changes to `validate_and_cost.py` — every gate is driven through its existing public API.
3. Targeted pytest to green; ruff clean.

## 7. Out of scope

- Stage 3 (FT cost) numbers themselves — requires the `chem-ft` environment, outside this repo's
  standard `make gates` flow; a separate spec run in that environment would be the natural follow-up.
- Whether the shared-threshold coupling (G2) SHOULD be split into two independent parameters — this
  spec measures and documents the current behavior, doesn't propose a fix.
- Systems beyond the two tested (H4 small, N2 large) — a third intermediate-size case is a natural
  follow-up, not attempted here.

## 8. Caveats and risks

- **R1 — the FT-cost stage is entirely untested by this spec** (out of scope, R1 above); "the
  pipeline degrades gracefully" is verified only for the taper/cross-check seam, not the
  cross-check/FT-cost seam, which would need the `chem-ft` environment to exercise for real.
- **R2 — a near-miss while probing this spec, caught not shipped:** the first draft of G2 tried to
  prove the threshold coupling by RAISING `qubit_dense_max_orb` on the large N2 system, which let
  `taper_hamiltonian` attempt `pauli_decompose` at 16 qubits and OOM-killed the process (see §3).
  `qubit_dense_max_orb` is not a safe knob to raise past a system's own `norb` on a genuinely large
  active space — a caution worth recording for anyone else composing `validate_and_cost`.
- Honest limitation: two systems, one parameter (`qubit_dense_max_orb`) examined for coupling; other
  parameters (`fci_max_dim`, `krylov_max_dim`, `target_mHa`, `tol_mHa`) are not examined for similar
  hidden couplings here.

## 9. Deliverables

- `tests/test_validate_and_cost_composition_spec.py` — gates G1-G4 (no library code changes).
- Results summary (with the R1 caveat) in the PR description / BACKLOG entry.
