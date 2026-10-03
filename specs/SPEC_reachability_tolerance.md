# SPEC: "reachable" means two different things in the certified arc

**Status:** IMPLEMENTED once gates green. Found while verifying the `krylov_refine` chained-overlap
hypothesis (specs/BACKLOG.md, Open → Certified-bounds arc): the chained bound appeared to *violate*
its reference, and the cause was not the bound.

---

## 1. Goal

Every module in the certified arc defines the HF-**reachable** sector by thresholding the HF
amplitude, `|⟨HF|ψ_k⟩|² > tol`. The threshold is a magic number written independently in each file,
and **it is not the same number**:

| tol | sites |
|---|---|
| `1e-10` | `certified_gaps.py:107`, `hf_overlap_certificate.py:76,99`, `certified_dipole.py:118`, `certified_noise.py:80` |
| `1e-8` | `hf_overlap_subspace.py:41` (`_REACHABLE_TOL`), and the ODMD/MSD/Trotter family |

**The falsifiable claim: there exists a physical geometry at which the two thresholds select
DIFFERENT ground states, so two specs in the same arc certify different targets while appearing to
be compared.** `SPEC_hf_overlap_certificate` (1e-10) and `SPEC_hf_overlap_subspace` (1e-8) are
exactly that pair — the head-to-head "the d=1 certificate is vacuous while the d=2 block rescues it"
is the latter's headline finding.

## 2. Background and honest framing

Witness, square H₄ at side **a = 1.1 Å** (the geometry from `hf_overlap_subspace.py`'s own
`__main__` sweep), STO-3G:

| threshold | selected ground state | HF overlap |
|---|---|---|
| `1e-10` | E = −4.556211 | **2.25e-5** |
| `1e-8`  | E = −4.404505 | **0.667** |

A level with HF amplitude² = 5.07e-10 sits between the two thresholds. It is *reachable* to one
module and *unreachable* to the other, and it is the lowest such level — so it becomes "the ground
state" for one and not the other. The overlap being certified differs by **4 orders of magnitude**.

**What we can claim:** the inconsistency exists, it is consequential at a real geometry already in
the repo's own sweep, and any d=1-vs-d=2 comparison at that geometry compares different targets.

**What we cannot claim — and this spec deliberately does NOT do.** We do **not** unify the constant.
Picking a value decides which of two recorded findings is right, and changes results in at least one
shipped spec. That is a scientific decision with a reference behind it, not a refactor, and it is
recorded as a follow-up rather than smuggled into a tolerance edit. This spec makes the divergence
**visible and gated** so it cannot silently widen.

~~**Not claimed:** that either value is wrong; that any recorded number is wrong. Both are internally
consistent — they answer different questions.~~ **← FALSIFIED 2026-07-31, see §2b.**

## 2b. FALSIFICATION — the physics framing above was wrong

The observation in §2 stands: the two thresholds do select different states. **The interpretation
did not.** "Neither value is wrong, they answer different questions" is false. Three measurements:

**(i) The disputed amplitude is an SCF convergence residue, not an overlap.** Varying only
`PySCFDriver(conv_tol=...)` at the witness geometry, the *eigenvalue is unchanged to 10 digits*
while the amplitude moves **19 orders of magnitude**:

| conv_tol | E₀ (electronic) | p₀ = \|⟨HF\|ψ₀⟩\|² |
|---|---|---|
| 1e-6 | −4.5562107647 | 1.49e-09 |
| 1e-9 *(driver default)* | −4.5562107647 | **5.07e-10** |
| 1e-11 | −4.5562107647 | 1.14e-10 |
| 1e-13 | −4.5562107647 | **1.53e-28** |

A physical overlap does not depend on how tightly the SCF was converged. Gated as **G6**.

**(ii) The mechanism: the state is symmetry-forbidden.** Square H₄ is D2h. The RHF MO irreps are
`[Ag, B2u, B3u, Ag]` with occupation (2,2,0,0), so the **HF determinant is Ag**. Symmetry-resolved
FCI in the tightly-converged basis:

| irrep | E₁ | \|c_HF-det\|² |
|---|---|---|
| B1g | −1.95159401 | **0.000e+00** (exactly) |
| Ag | −1.79988864 | 0.4451 |

The level `tol=1e-10` admits is the **B1g** ground state, whose true HF overlap is **exactly zero**.
The level `tol=1e-8` selects is the **Ag** ground state — the physically correct target. Gated as
**G7**. So at this witness **1e-10 is wrong and 1e-8 is right**, and the certified arc has been
certifying an overlap with a state Hartree–Fock cannot reach.

**(iii) But "just pick 1e-8" is also wrong.** Scanning the same family, the artifact amplitude is
**bistable**, not smooth: it sits at ~1e-29 where the driver finds the symmetric SCF solution, and
jumps to ~0.47 where it instead falls into a *lower* broken-symmetry RHF solution — and at
**a = 1.190 Å it lands at 1.4e-8, above the looser threshold too**. Gated as **G8**. **No fixed
constant separates physics from SCF residue.** The correct fix is a symmetry/sector-aware
reachability test, not a better number.

**Consequence for the original G1.** It pinned `ov10 ≈ 2.25e-5` "so a silent drift is loud". That
number *is* the drift — an SCF convergence residue that would break on any pyscf/qiskit-nature bump
or a change to the driver's default `conv_tol`. The pin is removed; the durable statement
(`ov10 < 1e-3`) is kept.

**Consequence for `SPEC_subspace_floor_resolvability` (PR #22).** Its recorded mechanism — "a
~1e-4-amplitude reachable level near the cluster boundary" — is very likely this same
symmetry-forbidden level seen through a looser SCF. That spec's mechanism statement deserves
re-examination on these grounds; its *conclusion* (the guard is insufficient, oracle mode is the
rigorous path) is untouched, since a spurious level breaks the floor guard just as effectively as a
real one. Filed as a follow-up rather than asserted here.

**Connection to a prior finding.** `SPEC_subspace_floor_resolvability` (PR #22) killed the
"fail-safe" claim for the self-mode floor and identified the mechanism as *"a ~1e-4-amplitude
reachable level near the cluster boundary"*. This is the same phenomenon one rung lower: a
near-threshold level whose classification is tolerance-dependent. That spec found it corrupting a
floor; this one finds it corrupting the definition of the target.

## 3. Approach

Pure dense linear algebra on an 8-qubit system — no Krylov, no fitting. Diagonalize, threshold the
HF amplitudes at both values, compare the selected index. **Reference:** exact `eigh`; the claim is
about which index a threshold picks, so the reference is the amplitude spectrum itself.

Expose one documented constant, `REACHABLE_TOL_CERTIFIED = 1e-10`, naming the certified arc's
prevailing value **without changing any call site**, so the follow-up has something to unify *to*
and the gates have something to compare *against*.

## 4. Public interface

```
hf_overlap_certificate.REACHABLE_TOL_CERTIFIED : float = 1e-10   # documentation + gate anchor only
```

No behaviour change. No call site is edited.

## 5. Acceptance criteria (validation gates)

`tests/test_reachability_tolerance_spec.py` — pure, ~8 qubits, seconds. **G1–G3 are
Linux-reference-only (skipped elsewhere) — see §10.**

- **G1 — the witness (DEFINITION OF DONE).** At square H₄ a = 1.1 Å the two thresholds select
  different lowest-reachable eigenstates, and the resulting HF overlaps differ by > 1000×. Killed if
  they select the same index — then the divergence is cosmetic.
- **G2 — the two shipped modules really do disagree.** `hf_overlap_certificate.exact_reachable_overlap`
  and `hf_overlap_subspace.exact_hf_subspace_overlap(..., 1)` — the d=1 references of the two specs
  in the head-to-head — return values differing by > 1000× at that geometry. This is the consequence,
  measured through the public functions rather than asserted from source.
- **G3 — the offending level is genuinely between the thresholds.** Some eigenstate has
  `1e-10 < |⟨HF|ψ⟩|² < 1e-8`. Killed if none does (the divergence would then have another cause).
- **G4 — the boundary: this is near-threshold, not pervasive.** At a ∈ {1.0, 1.2, 1.3, 1.4} the two
  thresholds agree on the lowest reachable index. Killed if they disagree everywhere — the finding
  would be "the arc is globally inconsistent", a much larger claim than the evidence supports.
- **G5 — the constant is pinned.** `REACHABLE_TOL_CERTIFIED == 1e-10` and matches the literal used
  by `hf_overlap_certificate.exact_reachable_overlap`, so a future edit to one without the other
  fails loudly.

## 6. Implementation plan (test-first)

1. `tests/test_reachability_tolerance_spec.py` encoding G1–G5 (initially failing).
2. Add the documented constant; wire `exact_reachable_overlap` to it (value unchanged).
3. `make gates`.

## 7. Out of scope

- **Unifying the tolerance.** Needs a decision on which value is physically right and a re-run of
  every affected recorded result. Filed as a follow-up hypothesis.
- The ODMD/MSD/Trotter family's 1e-8 sites — a different sector of the codebase with its own
  recorded results.
- The `krylov_refine` chained-overlap bound that surfaced this. Its verified numbers are recorded in
  the backlog; it cannot be gated until the target is unambiguous, which is precisely this finding.

## 8. Caveats and risks

- **R1 — a = 1.1 Å is one geometry.** The claim is existence ("there is a geometry where this
  bites"), not prevalence; G4 bounds it explicitly.
- **R2 — G2 goes through public functions**, so it will keep passing if both modules are later
  changed *together*. That is the intended semantics: the gate protects the *relationship*, not the
  literals. G5 covers the literal.
- The witness depends on the atom ordering in the geometry string (`H 0 0 0; H a 0 0; H a a 0;
  H 0 a 0`), taken verbatim from `hf_overlap_subspace.py`'s `__main__`. A different ordering gives a
  different HF determinant and may not exhibit it.

## 9. Deliverables

- `hf_overlap_certificate.py` — `REACHABLE_TOL_CERTIFIED`, no behaviour change.
- `tests/test_reachability_tolerance_spec.py` — G1–G5.

## 10. Platform dependence — G1–G3 are Linux-reference-only (chem-7rb, 2026-10-02)

On macOS 27 (Apple M3, scipy 1.15.3 on Accelerate, qiskit-nature 0.8.0, pyscf 2.13.1) G1, G2 and G3
fail and the other nine gates pass (`3 failed, 9 passed`). The cause is neither LAPACK noise in a
degenerate eigenspace nor a symmetry-broken SCF; it is where the SCF stops:

| conv_tol | p₀ Linux (§2b) | p₀ macOS 27 (measured) |
|---|---|---|
| 1e-6 | 1.49e-9 | 2.66e-6 |
| 1e-7 | — | 1.34e-6 |
| 1e-8 | — | 8.69e-8 |
| **1e-9 (default)** | **5.07e-10** | **1.2e-29** |
| 1e-11 | 1.14e-10 | 1.2e-29 |
| 1e-13 | 1.53e-28 | 1.2e-29 |

At a = 1.1 and the default, the lowest level's HF population on macOS is machine zero (the floor
moves between 1.2e-29 and 3.9e-29 from run to run), not degenerate: both thresholds select the Ag
level (overlap 0.6672), so G1 and G2 see a = b = 0.667 and G3's (1e-10, 1e-8) window is empty.
`scf_symmetry_status(1.1) = (False, −5.3e-15)`: the symmetric branch, as in the Linux session
(`test_symmetry_reachability_spec` passes 14/14 here). No `conv_tol` puts the lowest-level population
in the window on macOS — 2.7e-6, 1.3e-6, 8.7e-8, then machine zero — so no tolerance reproduces G1–G3
there. Regenerate: `t.residue_sweep()` in `tests/test_scf_conv_tol_spec.py` (docstring has the
command) and `uv run --no-sync python reachability.py` (p₀ and SCF branch over the a-sweep).

**Consequence.** G1–G3 assert a property of a platform-specific SCF stopping point, which §2b(i)
and G6 establish is unphysical. They are `skipif(sys.platform != "linux")` with the measured macOS
numbers in the reason; no threshold changed and Linux behaviour is untouched. The platform-independent
claim is carried by G6 (residue moves 19 orders with `conv_tol`), G7 (symmetry mechanism) and G8
(no fixed constant is safe), all of which pass on macOS. G4 and G5 also pass everywhere.

**G8 weakness (not changed).** On macOS a = 1.185 is a symmetry-broken solve (ΔE = −7.6e-2 Ha,
p₀ = 0.4625) and also counts as an offender (p₀ > 1e-8), so G8 would pass even if a = 1.190 did not
breach. The symmetric offender carries §2b(iii): a = 1.190 gives p₀ = 1.401e-8 (symmetric, ΔE =
−8.9e-16 Ha), matching the Linux 1.4e-8; a = 1.195 (symmetric) gives 3.8e-12. Requiring a
symmetric-SCF offender would make G8 stricter; not done because it cannot be verified on Linux here.
(Source: throwaway probe, `scf_symmetry_status` plus the default-build p₀ at these three a; not committed.)

**Not verified:** Linux. The Linux column rests on the numbers this spec recorded in §2b.
