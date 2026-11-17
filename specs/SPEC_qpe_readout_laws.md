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
  bounded above by `pi * lambda/2^t` at every tested `t` — a real, provable guarantee (the argmax
  decode picks the dyadic bin nearest the true phase, so its phase error is at most half a bin,
  `pi/2^t` out of a full `2*pi` circle of `2^t` bins; `E = lambda*cos(theta)` has
  `|dE/dtheta| <= lambda` everywhere by the mean value theorem, so `err(t) <= lambda*pi/2^t`
  unconditionally — see G2), not the docstring's smoother-sounding `~` claim; and the state-prep
  bottleneck obeys a clean multiplicative law `p_success(window) ~= f(t) * |<g|psi>|^2` where `f(t)`
  is essentially INDEPENDENT of the overlap value (verified across overlaps spanning 0.02 to 0.99, a
  50x range).
- **What you cannot claim:** that `f(t)` itself is monotonic in `t` — probing this spec found it is
  NOT (t=8: 0.960-0.989; t=10: 0.928-0.930, WORSE than t=8; t=12: 0.997) — a genuine artifact of how
  the fixed 3-bin window aligns with the dyadic phase grid at different `t`, recorded as a boundary
  rather than smoothed into a false monotonic story; nor that the measured *realized* ratio
  (`err(t)/(lambda/2^t)` tops out at 2.175 here) is itself derived or universal — closing the gap
  between that measurement and the provable `pi` ceiling to a tighter, per-system `pi*sin(theta_0)`
  envelope is `SPEC_qpe_precision_bound` (chem-0ke)'s job, not this spec's. (Revision note,
  chem-aob: an earlier version of this file called the then-constant `3` in G2 "provable" — it
  wasn't; `3` was a bare fit to this one measured sweep, with no derivation, and `3 < pi` means it
  was not even a valid ceiling in general. G2 now uses the actually-derived `pi`.)
- **Reference:** the exact CASCI ground energy and exact ground eigenvector (dense diagonalization
  of the qubit Hamiltonian) — the same reference `qpe_walk_readout.py`'s own `__main__` already uses.

## 3. Approach

Reuse `run_qpe`, `qpe_distribution`, `hartree_fock_vector` unmodified. Precision: run QPE with the
EXACT ground state as trial across `t = 4..15`, record `|E_est - CASCI|` at each `t`; check monotone
non-increasing and bounded by `pi * lambda/2^t`. State-prep: build synthetic trial states
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
- **G2 — a real, provable upper bound (not the docstring's bare `~`, and not merely an empirical
  fit).** `err(t) <= pi * lambda/2^t` at every tested `t`. This is a hard ceiling, not a measured
  constant: the argmax decode's phase error is at most half a dyadic bin, `pi/2^t` (a full circle
  of `2^t` bins spans `2*pi`), and `|cos(a) - cos(b)| <= |a - b|` for all `a, b` (mean value
  theorem, `|cos'| <= 1`), so `err(t) = lambda*|cos(theta_est) - cos(theta_0)| <= lambda*pi/2^t`
  always — true for any system, not fit to this one. *Measured max ratio
  `err(t) / (lambda/2^t)` = 2.175, comfortably inside the `pi` ≈ 3.1416 margin.*
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

- Deriving the *tight*, per-system realized ratio (2.175 here, not `pi`) or the `f(t)`
  non-monotonicity mechanism (both recorded as measured phenomena, not derived; the former is
  `SPEC_qpe_precision_bound`'s job, via a per-system `pi*sin(theta_0)` envelope — tighter than the
  universal `pi` ceiling G2 now proves, but system-dependent).
- Circuit-level QPE (ancilla-controlled walk-operator powers) — `run_qpe` is the exact
  Fejer-kernel-distribution simulation oracle, not a circuit.
- Systems beyond H2 CAS(2,2) — a natural follow-up once the pattern is established here.

## 8. Caveats and risks

- **R1 — one system only.** G2's `pi` ceiling is derived (G2's own entry) and therefore universal —
  it needs no generalization. What's measured on H2 CAS(2,2) only, and not shown to generalize here,
  is how close the realized ratio (2.175) sits to that ceiling, and the overlap-independence band
  widths; a different `lambda`/spectral gap could shift both. `SPEC_qpe_precision_bound` tests the
  realized-ratio generalization on LiH/N2 CAS(3,4). The falsifiable claims (G1's monotonicity, G3's
  overlap-independence) are structural properties of the argmax-Fejer-kernel construction and are
  expected to generalize, but that generalization is not tested here.
- **R2 — G4's non-monotonicity is an open, unexplained wrinkle**, not a bug fix. It is recorded as a
  boundary precisely because smoothing it into a false "f(t) increases with t" claim would be the
  kind of unfalsifiable overreach this repo's culture exists to prevent.

## 9. Deliverables

- `tests/test_qpe_readout_laws_spec.py` — gates G1-G4 (no library code changes).
- Results summary (with R1/R2 caveats) in the PR description / BACKLOG entry.
