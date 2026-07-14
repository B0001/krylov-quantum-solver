# SPEC: The qubitization walk operator recovers EVERY Hamiltonian eigenvalue, exactly twice

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

---

## 1. Goal

`qubitization_blueprint.py`'s `verify_qubitization` checks one direction of the spectral relation
`eig(W) = e^{±i arccos(E_k/lambda)}`: that every RECOVERED value `lambda*cos(theta)` lands near some
true eigenvalue of H. That is the wrong direction for what actually matters — a block-encoding bug
that silently DROPPED the ground state (or any eigenvalue) from `W`'s spectrum would pass that check
trivially, since it only asks "is every recovered value valid," never "is every true eigenvalue
recovered." This spec gates the direction that matters for QPE (which needs every eigenvalue of
interest, especially the ground state, to actually be phase-encodable), and sharpens it to an exact
counting invariant: every eigenvalue of H (with its natural degeneracy) contributes EXACTLY two
`W`-eigenphases (the `theta, -theta` pair), no more, no fewer. False if any eigenvalue of H is
missing from `W`'s spectrum, or if the multiplicity count is ever odd or unequal to exactly 2 per
eigenvalue slot.

## 2. Background and honest framing

- `qubitization_blueprint.py` is the fault-tolerant backbone `adapt_vqe.py` and `taper_qubits.py`
  already reuse (their JW Hamiltonian, Pauli decomposition) — this spec is the first to gate the
  module's own core claim, the thing everything downstream in the FT stack (`qpe_walk_readout.py`,
  `iterative_qpe.py`) rests on.
- **What you can claim if the gates pass:** the qubitization construction in this repo is a genuine,
  exact block encoding at the level that matters for QPE — no eigenvalue of the active-space
  Hamiltonian is silently unreachable through the walk operator's phases, and the `theta/-theta`
  pairing (a structural consequence of the LCU construction's real block structure) holds without
  exception, not just "on average."
- **What you cannot claim:** anything about T-gate cost or circuit-level qubitization (this module's
  `build_walk_operator` is the dense-matrix verification oracle, explicitly not the gate-synthesized
  circuit); that `lambda` is minimal (a separate question, already the subject of
  `SPEC_scdf_lambda.md`); scalability beyond exact diagonalization (`np.linalg.eigvals` on a
  `A*2^n`-dimensional dense matrix — this spec's systems keep that under 1000).
- **Reference:** the exact eigenvalues of the qubit Hamiltonian `H` (`np.linalg.eigvalsh`, dense) —
  the same reference `verify_qubitization` itself already uses, made bidirectional.

## 3. Approach

Reuse `build_qubit_hamiltonian`, `pauli_decompose`, `build_walk_operator` unmodified. For each
system: diagonalize `H` (`eigH`, WITH its natural degeneracy — not deduplicated, so a repeated
eigenvalue contributes one slot per occurrence) and `W` (`np.linalg.eigvals`, complex, since `W` is
unitary not Hermitian); recover `lambda*cos(theta)` for every `W` eigenphase (also not deduplicated).
For each entry in `eigH`, count how many `W`-recovered values fall within `1e-6` of it — the
multiplicity. Sum of multiplicities across all of `eigH` must equal exactly `2 * dim(H)` (every
eigenvalue slot gets exactly one `theta/-theta` pair, no eigenvalue missing, no extra collision).

## 4. Public interface

No new library code — this spec adds only test-file assertions around
`qubitization_blueprint.py`'s existing public functions (`build_qubit_hamiltonian`,
`pauli_decompose`, `build_walk_operator`), reused unchanged.

## 5. Acceptance criteria (validation gates)

- **G1 — bidirectional recovery (the direction the existing check doesn't test).** Every eigenvalue
  in `eigH` has at least one matching `W`-recovered value within `1e-6` Ha, on H2 CAS(2,2) and LiH
  CAS(2,2). *Measured max distance: H2 2.0e-15, LiH 1.0e-15 — machine precision, not just within
  tolerance.*
- **G2 — THE FINDING (definition of done): exactly two `W`-eigenphases per eigenvalue slot, summed
  over ALL of `eigH` (degeneracy included).** `sum(multiplicity(e) for e in eigH) == 2 * len(eigH)`
  exactly, on both systems. *Measured: H2 32 == 2*16; LiH 32 == 2*16.*
- **G3 — the `theta/-theta` pairing holds without exception.** Every individual eigenvalue's
  multiplicity count is even (never an odd, "unpaired" phase) on both systems.
  *Measured: H2 counts [2,4,4,6,6,6,4,4,...] (all even); LiH similarly all even.*
- **G4 — construction bookkeeping.** `lam` (the reported 1-norm) equals `sum(abs(c) for _, c in
  terms)` independently recomputed, and `W.shape[0] == 2**a * 2**n` exactly (ancilla dimension times
  system dimension) — pins the basic dimensions/normalization, not just the spectral relation.

> Definition of done: **G2**. G1 alone (every eigenvalue has SOME match) is necessary but not
> sufficient — a construction that collapsed every eigenvalue onto a single degenerate phase, or
> duplicated one eigenvalue's phase while dropping another's, could still pass G1.

## 6. Implementation plan (test-first)

1. Write `tests/test_qubitization_spectrum_spec.py` encoding G1-G4 (RED in the sense that the
   bidirectional/counting checks are new, even though `qubitization_blueprint.py`'s functions are
   not).
2. No changes to `qubitization_blueprint.py` — the check is a genuine external verification,
   deliberately not routed through `verify_qubitization`'s own one-directional logic.
3. Targeted pytest to green; ruff clean.

## 7. Out of scope

- Circuit-level (gate-synthesized) qubitization — `build_walk_operator` is the dense-matrix
  verification oracle only.
- T-gate/query cost claims (`lambda/epsilon` scaling) — prose in the module docstring, not gated
  here; a separate spec would need an actual QPE simulation.
- Systems beyond CAS(2,2) — `np.linalg.eigvals` on the full `A*2^n`-dim `W` is the bottleneck
  (256x256 / 512x512 here); a CAS(3,3)-scale system would push `dim(W)` into the low thousands,
  likely still tractable but untested — a natural follow-up, not attempted here.

## 8. Caveats and risks

- **R1 — the exact-multiplicity invariant (G2) is checked on two systems only**; a construction bug
  specific to a different term-count/ancilla-padding regime (e.g. `L` exactly a power of 2, so no
  ancilla padding is needed at all) is not ruled out. Both tested systems have `L` NOT a power of 2
  (`H2`: L=15, `LiH`: L=27), so the zero-padding path is exercised, but the no-padding path is not.
- Honest limitation: dense-matrix verification oracle, not circuit-level; two minimal-basis systems.

## 9. Deliverables

- `tests/test_qubitization_spectrum_spec.py` — gates G1-G4 (no library code changes).
- Results summary (with the R1 caveat) in the PR description / BACKLOG entry.
