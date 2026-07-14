# SPEC: Z2 qubit tapering preserves the FULL sector spectrum, independently verified

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

---

## 1. Goal

`taper_qubits.py` implements Z2 symmetry tapering and its `__main__` informally checks one number:
the tapered ground eigenvalue matches CASCI. That is the weakest possible falsifier — a tapering
bug that scrambled every excited state while leaving the ground state untouched would sail through
it. This spec gates the claim tapering actually makes: the Clifford-rotate-then-drop-a-qubit
procedure is a similarity transform followed by an exact restriction to one stabilizer eigenspace,
so it should preserve the ENTIRE eigenspectrum of that sector, not just its lowest eigenvalue — and
that restriction should be a REAL, non-vacuous one (a wrong sector gives a different spectrum).
False if the full-spectrum match fails anywhere it's tested, or if a wrong-sector projection turns
out indistinguishable from the correct one (which would mean the "gate" proves nothing).

## 2. Background and honest framing

- `taper_qubits.py` already reuses validated primitives (`qubitization_blueprint`'s JW Hamiltonian
  and Pauli decomposition, `adapt_vqe.hf_state`) — no new physics, only a stronger falsifier around
  a mechanism the module already implements.
- **What you can claim if the gates pass:** tapering is exact at the full-spectrum level (not just
  verified on the ground state) on three systems, including one OPEN-SHELL case (H3 radical,
  CAS(3,3), spin=1) `taper_qubits.py`'s own `__main__` never exercises — its `reference()` helper
  hardcodes closed-shell RHF; and the qubit count removed always equals exactly the number of
  independent Z-type symmetry generators found, an internal-consistency check on the reduction
  mechanism itself, not just its output energy.
- **What you cannot claim:** that tapering scales to systems this repo actually cares about at that
  qubit count. Probing this spec found `pauli_decompose` (used both by `taper_qubits.py` internally
  and by this spec's independent verification) costs ~1000s at 8 qubits (H4 CAS(4,4) — the LARGEST
  case `taper_qubits.py`'s own `__main__` already runs) against ~0.02s at 4 qubits and ~0.3s at 6 —
  a roughly 16x-per-qubit blowup consistent with its exponential (4^n Pauli strings, each needing an
  O(2^n)-size trace) construction. This spec's gates therefore stop at 6 qubits by design, and that
  ceiling — not scaled up silently, not hidden — is itself part of the recorded finding.
- **Reference:** an INDEPENDENTLY constructed sector-projected Hamiltonian — built directly from the
  computational-basis parity of each Z-type symmetry generator on the reference state (Z-type Pauli
  strings are diagonal, so this needs no Clifford rotation at all) — NOT `taper_hamiltonian`'s own
  Clifford-walk code path. An honest external check, not a round-trip through the same machinery.

## 3. Approach

For each system: build the qubit Hamiltonian and HF state (`qubitization_blueprint`,
`adapt_vqe.hf_state`), find the Z-type symmetry generators (`taper_qubits.find_symmetries`,
`_gf2_independent`, unmodified). Independently: for each generator `g` (a pure Z-string), read its
sector eigenvalue directly off the reference state, `s = <ref|Z_g|ref>` (a computational-basis
parity, no rotation needed); build the projector `P = prod_i (I + s_i Z_gi)/2`; diagonalize `P`'s
range to get the sector-projected spectrum. Compare (sorted, `atol=1e-7`) against
`np.linalg.eigvalsh(taper_hamiltonian(...)["H_tapered"])`. Separately, repeat the same independent
construction with a DIFFERENT reference state (a different computational basis string) to confirm
the two spectra genuinely differ — proof the check has teeth.

## 4. Public interface

No new library code — this spec adds only test-file helpers (the independent projector
construction) around `taper_qubits.py`'s existing public functions (`find_symmetries`,
`_gf2_independent`, `taper_hamiltonian`), reused unchanged.

## 5. Acceptance criteria (validation gates)

- **G1 — full-spectrum equivalence, independently verified.** On H2 CAS(2,2), LiH CAS(2,2), and H3
  radical CAS(3,3) (nelec=(2,1), spin=1, ROHF): the sorted full spectrum of `H_tapered` matches the
  sorted spectrum of the independently-built sector projection to `atol=1e-7`.
  *Measured max diff: H2 1.8e-15, LiH 6.7e-16, H3 (exact match to machine precision).*
- **G2 — the reduction mechanism is internally consistent.** `n_qubits_original - n_qubits_tapered
  == len(z_indep)` (the independent Z-symmetry count) on all three systems — the qubits removed
  match the symmetries found exactly, not a smaller or inconsistent number.
  *Measured: H2 4-1=3=|z_indep|; LiH 4-2=2=|z_indep|; H3 6-3=3=|z_indep|.*
- **G3 — THE FINDING (definition of done): the check is non-vacuous.** On H2 and LiH, the SAME
  independent-projector construction using a WRONG reference state (a different computational basis
  string, flipping the sector of at least one generator) gives a spectrum that does NOT match
  `H_tapered`'s — proving G1 actually discriminates correct tapering from incorrect, rather than
  trivially matching any sector. *Measured: differs on both tested systems.*
- **G4 — open-shell scope, newly exercised.** All of G1-G3 hold on H3 radical CAS(3,3) (spin=1,
  na != nb) — `taper_qubits.py`'s existing `__main__` only ever runs closed-shell RHF systems; this
  is the first check of tapering on an open-shell reference.

> Definition of done: **G3**. G1 alone (spectrum match) is exactly the kind of "a gate that cannot
> fail" this repo's culture warns against without G3 to prove it CAN fail.

## 6. Implementation plan (test-first)

1. Write `tests/test_taper_spectrum_spec.py` encoding G1-G4 (RED — none of the independent-check
   helpers exist yet; the systems/comparisons are new even though the underlying `taper_qubits.py`
   functions are not).
2. No changes to `taper_qubits.py` — the independent verification lives entirely in the test file
   (deliberately: reusing `taper_hamiltonian`'s own internals to "verify" itself would not be an
   independent check).
3. Targeted pytest to green; ruff clean.

## 7. Out of scope

- Systems beyond 6 qubits — `pauli_decompose`'s exponential cost (§2) makes them impractical for a
  CI gate; a follow-up would need a non-exhaustive Pauli-extraction path (e.g. reading JW terms
  directly from the fermionic operators instead of decomposing the dense matrix) before scaling up.
- Composing tapering with a downstream method (e.g. running `adapt_vqe` on the tapered, smaller
  register) — the tapered pool/operator mapping is a separate, nontrivial design question, not
  attempted here.
- X/Y-type (non-diagonal) symmetry generators — `taper_qubits.py` currently only taper Z-type ones
  (`z_syms` in its own code); whether it should is a separate question.

## 8. Caveats and risks

- **R1 — `pauli_decompose`'s exponential scaling is a real, load-bearing limitation of
  `taper_qubits.py` as it stands**, not a probing artifact: it is used inside `taper_hamiltonian`
  itself (every call to `taper_hamiltonian` pays this cost), so the existing H4 CAS(4,4) case in its
  own `__main__` already takes on the order of 15-20 minutes to run — worth flagging in the BACKLOG
  even though fixing it is out of scope here.
- Honest limitation: 3 systems, all minimal-basis STO-3G, small active spaces; Z-type symmetries
  only.

## 9. Deliverables

- `tests/test_taper_spectrum_spec.py` — gates G1-G4 (no library code changes).
- Results summary (with the R1 timing finding) in the PR description / BACKLOG entry.
