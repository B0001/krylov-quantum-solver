# SPEC: QPE readout — the precision law is a staircase, and the state-prep law is overlap-independent

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

---

## 1. Goal

`qpe_walk_readout.py`'s docstring makes two quantitative claims about textbook QPE on the
qubitization walk operator, neither ever CI-gated: (1) "energy resolution ~ lambda / 2^t" and (2)
"a trial state's ground-state overlap sets the QPE success probability." This spec gates BOTH,
correcting the first (the literal reading — smooth `~1/2^t` scaling of the realized point-estimate
error — is false; the truth is a monotonically non-increasing staircase, bounded above but not
smooth) and sharpening the second into a precise, falsifiable law: the ratio of window success
probability to overlap depends almost entirely on `t`, not on the overlap value itself, across a
50x range of overlaps. False if the error is ever non-monotonic in `t`, or exceeds the bound; false
if the success/overlap ratio varies substantially across overlap values at fixed `t`.

## 2. Background and honest framing

- `qpe_walk_readout.py` already reuses validated primitives (`qubitization_blueprint`'s exact JW
  Hamiltonian and walk-operator eigenphase relation) — no new physics, only falsifiers around a
  simulation the module already runs and prints, never gates.
- **What you can claim if the gates pass:** the argmax point-estimate error from exact-Fejer-kernel
  QPE is a MONOTONICALLY NON-INCREASING function of phase bits `t` (adding resolution never hurts),
  bounded above by a small constant multiple of `lambda/2^t` at every tested `t` — a real, provable
  guarantee, not the docstring's smoother-sounding `~` claim; and the state-prep bottleneck obeys a
  clean multiplicative law `p_success(window) ~= f(t) * |<g|psi>|^2` where `f(t)` is essentially
  INDEPENDENT of the overlap value (verified across overlaps spanning 0.02 to 0.99, a 50x range).
- **What you cannot claim:** that `f(t)` itself is monotonic in `t` — probing this spec found it is
  NOT (t=8: 0.960-0.989; t=10: 0.928-0.930, WORSE than t=8; t=12: 0.997) — a genuine artifact of how
  the fixed 3-bin window aligns with the dyadic phase grid at different `t`, recorded as a boundary
  rather than smoothed into a false monotonic story; nor that the precision bound's constant
  (measured ~2.2x here) is universal — it is a property of the argmax-over-Fejer-kernel estimator on
  this system, not derived analytically.
- **Reference:** the exact CASCI ground energy and exact ground eigenvector (dense diagonalization
  of the qubit Hamiltonian) — the same reference `qpe_walk_readout.py`'s own `__main__` already uses.

## 3. Approach

Reuse `run_qpe`, `qpe_distribution`, `hartree_fock_vector` unmodified. Precision: run QPE with the
EXACT ground state as trial across `t = 4..15`, record `|E_est - CASCI|` at each `t`; check monotone
non-increasing and bounded by `C * lambda/2^t`. State-prep: build synthetic trial states
`sqrt(p)*ground + sqrt(1-p)*(random state orthogonal to ground)` at `p in {0.99, 0.9, 0.7, 0.5, 0.3,
0.1, 0.05, 0.02}`, run QPE at fixed `t in {8, 10, 12}`, record `p_success(window) / measured_overlap`
and check the band width across the `p` sweep at each fixed `t`.

## 4. Public interface

No new library code — this spec adds only test-file assertions around `qpe_walk_readout.py`'s
existing public functions (`run_qpe`, `hartree_fock_vector`), reused unchanged.

## 5. Acceptance criteria (validation gates)

- **G1 — the precision error is monotonically non-increasing in `t`.** On H2 CAS(2,2), exact ground
  state trial, `t = 4..15`: `err(t+1) <= err(t) + 1e-15` at every step.
  *Measured: strictly staircase-shaped, e.g. 5.6e-2 (t=4-6) -> 4.0e-2 (t=7) -> 8.9e-3 (t=8-9) ->
  3.2e-3 (t=10) -> ... -> 1.4e-4 (t=12-15); never increases.*
- **G2 — a real, provable upper bound (not the docstring's bare `~`).** `err(t) <= 3 * lambda/2^t`
  at every tested `t`. *Measured max ratio `err(t) / (lambda/2^t)` = 2.175, comfortably inside the
  3x margin.*
- **G3 — THE FINDING (definition of done): the success/overlap ratio is essentially
  overlap-independent at fixed `t`.** Across the overlap sweep (0.02 to 0.99, a 50x range) at each
  of `t = 8, 10, 12`, the band width of `p_success(window)/measured_overlap` is `< 0.05`.
  *Measured band widths: t=8: 0.960-0.989 (0.029); t=10: 0.928-0.930 (0.002); t=12: 0.997-0.997
  (0.0001).*
- **G4 — boundary, recorded not smoothed over: the ratio is NOT monotonic in `t`.** The
  representative ratio at each tested `t` stays bounded in `[0.85, 1.05]` (a sanity regime bound),
  but t=10's ratio is measurably LOWER than t=8's — proving G3's "independent of overlap" is not
  also "independent of/monotonic in t," a distinct and not-yet-explained wrinkle.
  *Measured representative (min) ratios: t=8: 0.960; t=10: 0.928; t=12: 0.997.*

> Definition of done: **G3**. G1/G2 replace the docstring's imprecise prose with an actually-checked
> bound; G4 keeps G3 honest about what it does NOT also claim.

## 6. Implementation plan (test-first)

1. Write `tests/test_qpe_readout_laws_spec.py` encoding G1-G4 (RED in the sense that these checks
   are new, even though `qpe_walk_readout.py`'s functions are not).
2. No changes to `qpe_walk_readout.py` — a genuine external verification of behavior the module
   already exhibits but never gated.
3. Targeted pytest to green; ruff clean.

## 7. Out of scope

- An analytic derivation of the precision bound's constant or the `f(t)` non-monotonicity mechanism
  (both recorded as measured phenomena, not derived).
- Circuit-level QPE (ancilla-controlled walk-operator powers) — `run_qpe` is the exact
  Fejer-kernel-distribution simulation oracle, not a circuit.
- Systems beyond H2 CAS(2,2) — a natural follow-up once the pattern is established here.

## 8. Caveats and risks

- **R1 — one system only.** The precision-bound constant (3x margin over a measured 2.175x) and the
  overlap-independence band widths are measured on H2 CAS(2,2); a different `lambda`/spectral gap
  could shift both. The falsifiable claims (G1's monotonicity, G3's overlap-independence) are
  structural properties of the argmax-Fejer-kernel construction and are expected to generalize, but
  that generalization is not tested here.
- **R2 — G4's non-monotonicity is an open, unexplained wrinkle**, not a bug fix. It is recorded as a
  boundary precisely because smoothing it into a false "f(t) increases with t" claim would be the
  kind of unfalsifiable overreach this repo's culture exists to prevent.

## 9. Deliverables

- `tests/test_qpe_readout_laws_spec.py` — gates G1-G4 (no library code changes).
- Results summary (with R1/R2 caveats) in the PR description / BACKLOG entry.
