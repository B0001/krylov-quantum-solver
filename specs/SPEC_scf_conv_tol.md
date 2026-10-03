# SPEC: does the SCF stopping tolerance move any number the gates consume? (chem-ayr)

**Status:** PRE-REGISTERED 2026-10-02. Every threshold in §5 was fixed before
`tests/test_scf_conv_tol_spec.py` was written or run; results go in a later section and may not
edit §5. A threshold that turns out wrong is recorded as a finding, not retuned.

---

## 1. Goal

chem-ayr was filed (2026-09-28) as *"`build_molecular_hamiltonian` hardcodes `PySCFDriver` with no
`conv_tol`; thread it through, re-run the affected gates at 1e-13, die if no recorded number
moves."* The threading has since landed (`conv_tol: float = 1e-9`, the driver's own default, so no
call site changed), which turns the question into a measurable one:

**Claim.** Through the public builder, tightening `conv_tol` from the default 1e-9 to
`reachability.TIGHT_SCF_CONV_TOL` (1e-13) (a) leaves every clean gate geometry unchanged, (b) moves
the symmetric-SCF square-H₄ artifact witnesses (a = 1.35, 1.19) by far more than the gates'
tolerances, and (c) collapses the artifact residue to machine zero there.

**Verdict if (a)+(b)+(c) hold:** a *documentation* outcome, not a default change. `conv_tol` stays an
opt-in knob; the artifact is real but confined to exact-symmetry geometries; no clean-set number in
the certified arc needs re-running. Changing the default would break by design every gate that pins
the artifact and would change chemcheck's frozen hashes.

**Kill rules (pre-registered).** (a) fails → tightening moves recorded numbers: it is a blast-radius
change and the "do not change the default" recommendation flips to "re-run every consumer first".
(b) fails → the bead's own kill rule: the artifact is practically unreachable from the public API
and this is a note, not a finding. (c) fails → `TIGHT_SCF_CONV_TOL` is not tight enough.

## 2. Background and honest framing

- `SPEC_reachability_tolerance` §2b: at square H₄ a = 1.1 the amplitude `|⟨HF|ψ₀⟩|²` of a
  symmetry-forbidden level is SCF convergence residue (moves 19 orders with `conv_tol`, eigenvalue
  unchanged). `SPEC_chained_overlap` R2b/R3 excluded a = 1.10 and 1.35 for it.
- PySCF stops when `|ΔE| < conv_tol` **and** `‖g_orb‖ < sqrt(conv_tol)` (`scf/hf.py` kernel). At the
  default the orbital-gradient norm at stop is therefore up to ~3.2e-5, at 1e-13 up to ~3.2e-7. An
  honest SCF-error effect on an overlap is O(that); the artifact is O(1) (a forbidden level admitted
  as the "ground state"). The G1 thresholds below are set from that scale, not from any measured
  value.
- **Disclosure.** An earlier unrecorded one-off probe (`git show 5c5c6db:scratch_conv_tol_probe_ayr.py`;
  its numbers were passed to this work in a hand-off, not stored in the repo) already suggested the
  outcome. The thresholds here are derived from the PySCF criteria above and are looser than what
  that probe saw; the gate re-measures everything from scratch, both tolerances in one process.
- **Stale text found while preparing this spec** (corrected in the same change, see §9):
  `tests/test_chained_overlap_spec.py` and `SPEC_chained_overlap.md` §2 say the references are
  "built at conv_tol=1e-13". The code builds at the default; only the claim was wrong.

**Can claim** if the gates pass: the verdict above, for the geometries below.
**Cannot claim:** that no number elsewhere moves (12-qubit H₆, CAS tiers, other call sites are not
re-measured); anything about Linux (the gate is run on macOS only, §8 R1).

## 3. Approach

For each geometry build `mh_default = build_molecular_hamiltonian(atom)` and
`mh_tight = build_molecular_hamiltonian(atom, conv_tol=TIGHT_SCF_CONV_TOL)` **in the same process**,
dense-diagonalize both (`reachability._dense_hf_projection`, ZHEEVR-safe) and compare ground-state
energy `mh.ground_state_energy()`, `mh.hf_energy`, `hf_overlap_certificate.exact_reachable_overlap`,
the lowest reachable level's energy, the lowest level's HF population p₀, and the Ritz leakage
‖P_unreach v‖ of a depth-6 `QuantumKrylovSolver` ground Ritz vector (the chained-overlap R2b
quantity). **Reference:** the tight build, whose p₀ is machine zero where the level is
symmetry-forbidden (`SPEC_reachability_tolerance` G6/G7 gives the symmetry-resolved ground truth).

**Geometries (subset, by cost).** Clean set (9): H₂ 0.74 / 2.0 Å; linear H₄ at 0.9 / 1.0 / 2.0 Å;
square H₄ a = 1.0 / 1.05 / 1.2 / 1.4. Witnesses (2): square H₄ a = 1.35 and 1.19 (symmetric SCF,
artifact present at the default on every platform recorded so far). **Not gated:** a = 1.10 (its
default-tolerance residue is platform-dependent, §8 R1); linear H₆ (12 qubits: minutes, and dense
`eigh` is unreliable on this Mac); LiH/N₂ CAS tiers.

## 4. Public interface

None. No source file changes: `build_molecular_hamiltonian(..., conv_tol=1e-9)` and
`reachability.TIGHT_SCF_CONV_TOL = 1e-13` already exist.

## 5. Acceptance criteria (validation gates) — thresholds FIXED

`tests/test_scf_conv_tol_spec.py`, 8 qubits and below, a few minutes.

- **G1 — clean set does not move (claim a).** For each of the 9 clean geometries:
  |ΔE₀| < 1e-9 Ha (full-space spectrum is orbital-rotation invariant; 1e-9 is the default energy
  tolerance); |ΔE_HF| < 1e-8 Ha (10× that tolerance); |Δ exact_reachable_overlap| < 1e-4 (≈ 3×
  sqrt(1e-9), the default gradient criterion); the lowest reachable level (p > 1e-10) has the same
  index in both builds.
- **G2 — witnesses move, eigenvalue does not (claim b; DEFINITION OF DONE).** At a = 1.35 and 1.19:
  |ΔE₀| < 1e-9 Ha (same eigenvalue) **and** at the default p₀ > `REACHABLE_TOL_CERTIFIED` (1e-10)
  with exact_reachable_overlap < 1e-3 (the forbidden level is admitted) **and**
  |Δ exact_reachable_overlap| > 0.1 **and** |Δ lowest-reachable energy| > 10 mHa.
- **G3 — residue collapses at TIGHT (claim c).** At a = 1.35 and 1.19: p₀(tight) < 1e-20 and tight
  Ritz leakage < 1e-12 (chained-overlap `LEAK_TOL`), while default leakage > 1e-12 (so that
  exclusion was warranted at the default).
- **G4 — the default is pinned at 1e-9.** `inspect.signature(build_molecular_hamiltonian)`'s
  `conv_tol` default `== 1e-9`. Changing it is a deliberate blast-radius change (re-freeze chemcheck
  T0–T3 hashes, move the gates that pin the artifact) and must revisit this verdict.

## 6. Implementation plan (test-first)

1. Commit this spec. 2. Write the gate encoding §5 verbatim. 3. Run it; record §10 without editing
§5. 4. Correct the stale conv_tol=1e-13 text in the chained-overlap gate and spec.

## 7. Out of scope (filed as follow-ups, not done)

- Changing the default `conv_tol`.
- Threading `conv_tol` through the callers that do not forward it: `pipeline.py:115`,
  `certchem/core.py:94`, `chemcheck/tiers.py:109`, `certkit_bridge.py:190`, `screening_loop.py:136`,
  and `build_dipole_operators` (`molecular_hamiltonian.py:212`), which runs its own SCF with no
  `conv_tol` at all.
- Re-admitting a = 1.35 to the chained-overlap gate at `TIGHT_SCF_CONV_TOL` (G3 makes it a candidate).

## 8. Caveats and risks

- **R1 — platform.** The SCF stopping point is platform-dependent: at a = 1.10 and the default the
  forbidden-level residue is 5.07e-10 on the Linux freeze and ~1e-29 on macOS (see
  `SPEC_reachability_tolerance` §10). a = 1.35 and 1.19 are used because their default-tolerance
  numbers agree with the Linux records (`SPEC_chained_overlap` R2b: leakage 1.6e-6 at 1.35;
  `SPEC_reachability_tolerance` §2b(iii): p₀ = 1.4e-8 at 1.19). This gate has **not** been run on
  Linux.
- **R2 — scale.** G1's thresholds mean "nothing at the scale of the SCF criterion moved", not
  bitwise identity. 9 + 2 geometries at ≤ 8 qubits do not cover the whole certified arc.
- **R3 — a = 1.19 sits next to a symmetry-breaking geometry** (RHF is reported to break symmetry at
  a = 1.185 on this machine; to be re-measured, not assumed). It is a witness only while the SCF
  stays on the symmetric branch, which G2's default-tolerance conditions would detect.

## 9. Deliverables

- `specs/SPEC_scf_conv_tol.md` (this file) — pre-registered, then results (§10).
- `tests/test_scf_conv_tol_spec.py` — G1–G4.
- Corrected text: `tests/test_chained_overlap_spec.py` (docstring + exclusion comment),
  `specs/SPEC_chained_overlap.md` §2, and the chem-ayr entry of `specs/BACKLOG.md`.

## 10. Results (2026-10-02; section 5 above is unedited since its pre-registration commit 6e7c8b6)

Apple M3, macOS 27.0.1, Python 3.13.14, numpy 2.4.6, scipy 1.15.3 (Accelerate), qiskit-nature
0.8.0, pyscf 2.13.1; 2 BLAS threads, shared machine. `tests/test_scf_conv_tol_spec.py`:
**14 passed in 5.2 s** (G1 ×9, G2 ×2, G3 ×2, G4). Regenerate the table below (not a gate):
`uv run --no-sync python -c "import sys; sys.path.insert(0, 'tests'); import test_scf_conv_tol_spec as t; t.table()"`.

| geometry | \|ΔE₀\| | \|ΔE_HF\| | overlap default → tight | p₀ default → tight | \|ΔE₀,reach\| | lowest reachable idx | reachable levels | Ritz leakage default → tight |
|---|---|---|---|---|---|---|---|---|
| H₂ 0.74 | 0 | 0 | 0.9936 → 0.9936 | 9.9e-1 | 0 | 0 / 0 | 2 / 2 | 3.5e-17 → 3.5e-17 |
| H₂ 2.0 | 0 | 0 | 0.8437 → 0.8437 | 7.1e-1 | 0 | 0 / 0 | 2 / 2 | 0 → 0 |
| H₄ chain 0.9 | 4.4e-15 | 5.5e-14 | 0.9769 (Δ 2.7e-9) | 9.5e-1 | 1.8e-15 | 0 / 0 | 12 / 12 | 2.8e-15 → 2.7e-15 |
| H₄ chain 1.0 | 3.6e-15 | 3.2e-13 | 0.9677 (Δ 6.8e-9) | 9.4e-1 | 6.2e-15 | 0 / 0 | 12 / 12 | 2.2e-15 → 5.4e-15 |
| H₄ chain 2.0 | 4.4e-15 | 0 | 0.6941 (Δ 1.3e-10) | 4.8e-1 | 4.4e-15 | 0 / 0 | 12 / 12 | 2.2e-13 → 3.0e-14 |
| H₄ square 1.0 | 5.3e-15 | 0 | 0.6898 (Δ 9.0e-15) | 4.8e-1 | 0 | 0 / 0 | 8 / 8 | 1.1e-14 → 6.5e-15 |
| H₄ square 1.05 | 3.4e-14 | 1.2e-14 | 0.6739 (Δ 7.7e-10) | 3.7e-32 → 1.7e-29 | 4.4e-15 | 4 / 4 | 8 / 8 | 9.3e-15 → 1.5e-14 |
| H₄ square 1.2 | 4.4e-15 | **2.4e-10** | 0.6791 (Δ 1.0e-9) | 4.6e-1 | 2.7e-15 | 0 / 0 | **9 / 8** | **5.1e-7** → 1.1e-14 |
| H₄ square 1.4 | 1.8e-15 | 8.9e-16 | 0.6637 (Δ 4.9e-15) | 4.4e-1 | 2.7e-15 | 0 / 0 | 8 / 8 | 1.4e-14 → 9.8e-15 |
| **H₄ square 1.35** | 8.9e-16 | 1.1e-9 | **7.80e-5 → 0.6128** | **6.1e-9 → 4e-30** | **0.1376 Ha** | 0 / **4** | 10 / 8 | **1.6e-6 → 1.8e-14** |
| **H₄ square 1.19** | 5.3e-15 | 2.3e-9 | **1.18e-4 → 0.6518** | **1.4e-8 → 3e-28** | **0.1496 Ha** | 0 / **4** | 10 / 8 | **2.6e-6 → 2.9e-14** |

**Verdict: (a), (b), (c) and the G4 pin all hold; no kill rule fired.** Margins against the
pre-registered thresholds: clean-set max |ΔE₀| 3.4e-14 (limit 1e-9), max |ΔE_HF| 2.4e-10 (1e-8),
max |Δ overlap| 6.8e-9 (1e-4), lowest reachable index unchanged in 9/9. Witnesses: overlap moves
0.61 / 0.65 (limit > 0.1), reachable ground energy moves 137.6 / 149.6 mHa (limit > 10 mHa) while
E₀ is unchanged to 9e-16 / 5e-15 Ha (the eigenvalue is untouched, only the HF-reachable
classification moves), p₀(tight) is machine zero (limit 1e-20; the floor varies run to run between
~1e-31 and ~3e-28), and the Ritz leakage drops from ≥ 1.6e-6 to ≤ 2.9e-14 (limit 1e-12).

**Unregistered observation, reported not gated (it does not touch the claims as pre-registered).**
At a = 1.2 (broken-symmetry RHF, lowest level physical) the default tolerance still admits one
symmetry-forbidden *excited* level into the reachable sector: 9 vs 8 reachable levels and Ritz
leakage 5.1e-7 vs 1.1e-14. The lowest reachable level, the overlap and every G1 observable are
unchanged, so G1 passes; but "only the witnesses move" is true of the pre-registered observables,
not of the sector size or R2b's leakage. This is a second instance of `SPEC_chained_overlap` R2b's
"the margin is a property of this set". a = 1.2 is not in the chained-overlap gated set.

**Not gated, a = 1.10 (platform-dependent).** Lowest-level p₀ vs `conv_tol` on this Mac
(`t.residue_sweep()`): 1e-6 → 2.7e-6, 1e-7 → 1.3e-6, 1e-8 → 8.7e-8, then 1.2e-29 for every
`conv_tol` ≤ 1e-9. The default build therefore shows no artifact here, whereas the Linux freeze
records 5.07e-10 (`SPEC_reachability_tolerance` §2b). See that spec's §10.

**Recorded verdict (chem-ayr).** Document, do not change the default. The builder already exposes
`conv_tol`; the artifact is real and large (reachable overlap 7.8e-5 → 0.613 at a = 1.35) but confined
to exact-symmetry geometries; no clean-set number needs re-running; the default stays 1e-9 (pinned by
G4) because changing it would break by design the gates that pin the artifact and change chemcheck's
frozen hashes. Follow-ups are listed in §7. **Not verified:** Linux; 12-qubit H₆ and CAS tiers.
