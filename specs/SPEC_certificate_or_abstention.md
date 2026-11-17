# SPEC: Certificate-or-abstention is a non-bypassable invariant on the energy return path

**Status:** DRAFT — gates written first (`tests/test_certificate_or_abstention_spec.py`).
Bead `portfolio-0pc`. Builds on `SPEC_certkit_return_paths.md` (`portfolio-cir`), which made
every return path *able* to carry a certificate and left it opt-in.

---

## 1. Goal

Claim: every public path of `hybrid_quantum_solver` that returns a ground-state energy hands
the caller one of two things and nothing else — the energy together with a certkit certificate
the pinned checker VERIFIED out of process, or an explicit abstention. No path returns a bare
number, and no argument, default or environment variable makes a path skip the checker while
still returning something that reads as a certified energy. Falsifiable: one public path that
returns a float, one abstained result whose number is readable under a name that does not say
`unverified`, one VERIFIED verdict the gate cannot reproduce from the certificate files, or
one unlisted public energy path, fails the gate.

## 2. Background and honest framing

- `portfolio-cir` added `certify=<dir>` to every return path. Without it energies came back
  uncertified, exactly as before: limit 1 of `SPEC_certkit_return_paths.md`, "opt-in, so
  bypassable". This spec removes the option.
- House style for invariants is `SPEC_certified_overlap_bounds.md` §4: non-bypassable, loud
  failure, an explicit VACUOUS result instead of a fabricated number. The same shape is used
  here: a missing verdict raises at construction, and a refused energy raises on read.
- **What we can claim if the gates pass:** a caller holding a result from any listed path
  can read its energy through exactly two kinds of name. `energy` (and `energies`,
  `computed_energy`, `convergence`) return the number only when the checker VERIFIED a
  certificate proving it an upper bound on the ground state, and raise `EnergyAbstained`
  otherwise. `unverified_energy` (and its siblings) always return the solver's number and
  say in their name that nothing vouches for it. Certification runs by default on every
  listed path. The energies themselves are the numbers the solver returned before this
  spec (§5, G2, and the suite-wide comparison in §8).
- **What we cannot claim:**
  1. **Python cannot stop a caller from reading `unverified_energy` and calling it an
     energy,** from forging a `Verdict`, or from calling a private method. The invariant is
     about what this package returns and how it names it, not about what a caller does
     next. G1 pins the public surface so the boundary cannot move unnoticed.
  2. **The listed paths, not every number in the package.** Out of scope, each for a stated
     reason (§7): classical references (`ground_state_energy`, `hf_energy`,
     `dmrg_reference.*`, the analytic model energies), `estimate_energy_aer` (the
     expectation value of a circuit the caller supplies), and the linear-algebra kernels
     `ritz_pairs`, `ritz_spectrum` and `solve_generalized_eig` (eigenvalues of matrices the
     caller supplies; they see no operator and no state, so nothing could certify them).
  3. **There is one uncertified door, and it is deliberate.**
     `QuantumKrylovSolver.uncertified_eigenstates` and the pipeline's
     `uncertified_estimates_only=True` return abstentions without running the checker. §3
     (d) says why they exist and why they are not a bypass.
  4. **Ground state only.** `ExcitedKrylovStep.energies[1:]` are Ritz values no certificate
     covers; certkit v1 has no claim kind for excited levels. They are released together
     with `energies[0]` and withheld together with it.
  5. **In a sweep the steps share one proof** (§3 c). A step's certificate may have been
     issued for a lower energy of the same sweep. It still proves that step's energy an
     upper bound; its ceiling is then below the step's energy, not equal to it.
  6. Everything `SPEC_certkit_return_paths.md` §2 lists still holds for a VERIFIED
     certificate: it is about the global minimum over every particle sector, the floor is
     only as good as the route, the hardware witness is rebuilt classically, and the checker
     runs out of process but in the producer's own venv.
  7. **Callers outside the package keep their own contracts.** `certchem`'s `Mode.FAST`
     returns a bare float by design (ADR-0004); it now takes that float from the uncertified
     door, so the choice is visible at its call site. Study scripts that write CSV keep
     their column names and gain a column with the checker's verdict.

## 3. Approach

One module, `hybrid_quantum_solver/certification.py`, owns the producer and the gate on the
number. Every listed path calls `certify_energies`, and every result type reads its energy
through `certified`.

**(a) How an abstention reaches the caller.** The result types (`KrylovStep`,
`ExcitedKrylovStep`, `SKQDStep`, `PipelineResult`) store the solver's number as
`unverified_energy` and the checker's `Verdict` as `certificate`. `certificate` is never
`None`: constructing a result without a `Verdict` raises `TypeError`. `energy` is a property.
It returns the number when the verdict is VERIFIED and raises `EnergyAbstained`, carrying the
checker's own line, when it is not. So code written against a certified energy cannot run on
an abstention by accident, and code that wants the raw estimate has to say `unverified`.
There is no third state: a crash of the checker is an abstention whose line starts `CRASH`.

**(b) Where the producer lives.** `certify_ground_state` and what it needs moved from the
top-level `certkit_bridge.py` into the package, and certkit moved from the `certkit` extra to
a core dependency. A lazily imported optional checker is a switch: without it installed,
every path would have had to choose between raising and returning an uncertified number. The
package now imports certkit at import time and fails loudly if it is missing. The import
rule of certkit's `INTEGRATION.md` is unchanged: interval, operators, schema and `pad_claim`
only, and the verdict comes from `certkit check` in a subprocess. `certkit_bridge.py` keeps
the regression cases (`run_case`) and re-exports the moved names.

**(c) Cost: one checker run per sweep, not one per energy.** A VERIFIED certificate for
energy `E` proves `lambda_min <= E + pad`, so it also proves every energy `E' >= E` of the
same operator an upper bound. `certify_energies` therefore submits the lowest energy of a
call first and, once one verifies, gives every energy at or above it the same verdict.
Energies below the first verified one each keep their own abstention. `convergence` on every
solver and the pipeline's convergence trace go through this, so a noiseless sweep costs one
checker run where `portfolio-cir` paid one per step. Measured cost of one run: 0.06 s at 4
qubits, 4.4 s at 8 (Temple, an O(n³) interval factorisation), 3.6 s at 12 (Gershgorin). The
baseline suite makes 11 eight-qubit sweeps returning 352 energies: 25 minutes of checker time
per energy, under a minute per sweep. A single `solve` still pays one run.

**(d) The uncertified door, and why it is not a bypass.** Two kinds of caller cannot pay
for a checker run per call. Analysis code reads the Ritz decomposition, not an energy claim:
the baseline suite makes 1,857 `eigenstates` calls, 1,635 of them at 8 to 12 qubits, which
would be about two hours of checker time for numbers nobody reads as certified. And the
checker is pure Python: one run is minutes at 16 qubits and hours at 20, so a study on a
large active space cannot certify every point. For both there is
`QuantumKrylovSolver.uncertified_eigenstates(krylov_dim, n_states)`, and the pipeline's
`uncertified_estimates_only=True`, which is built on it. It is not a bypass because of what
it returns: the same result types, already abstained, with a verdict line that reads
`ABSTAIN   not submitted to the checker: ...`. Their `energy` raises. The number is there
only as `unverified_energy`. So the door cannot produce anything that passes for a certified
energy, which is the only thing a bypass would be for. It is one method and one keyword,
both named `uncertified`, and G3 covers both. The certified paths themselves take no
argument that changes whether the checker runs (G1 pins their parameters).

**Reference:** exact `lambda_min` by dense diagonalisation, `mh.ground_state_energy() -
mh.energy_offset`, as in the sibling gates. Every VERIFIED verdict is also re-checked by the
gate itself from the certificate files.

## 4. Public interface

```
hybrid_quantum_solver.certification
    Verdict(name, rule, ok, line, lo, hi)             the checker's word on one certificate
    EnergyAbstained(RuntimeError)                     raised on reading a refused energy
    certified(value, verdict)                         -> value, or raises EnergyAbstained
    certify_energies(mh, candidates, out=None, names=None) -> [Verdict], lowest energy first
    certify_ground_state(mh, psi, energy, out)        -> Verdict   (moved from certkit_bridge)
    not_submitted(reason)                             -> the door's abstention Verdict

QuantumKrylovSolver.solve(krylov_dim, certificate_dir=None)               -> KrylovStep
QuantumKrylovSolver.convergence(max_dim, certificate_dir=None)            -> [KrylovStep]
QuantumKrylovSolver.solve_excited(krylov_dim, n_states=None, certificate_dir=None)
                                                                          -> ExcitedKrylovStep
QuantumKrylovSolver.eigenstates(krylov_dim, n_states=None, certificate_dir=None)
                                                              -> (ExcitedKrylovStep, states)
QuantumKrylovSolver.uncertified_eigenstates(krylov_dim, n_states=None)
                                                   -> (ExcitedKrylovStep [abstained], states)
TrotterKrylovSolver.solve / .convergence(..., certificate_dir=None)       -> KrylovStep(s)
HardwareKrylovSolver.solve / .convergence(..., certificate_dir=None)      -> KrylovStep(s)
SampleKrylovSolver.solve / .convergence(..., certificate_dir=None)        -> SKQDStep(s)
run_geometry / run_from_integrals(..., certificate_dir=None, uncertified_estimates_only=False)
                                                                          -> PipelineResult

KrylovStep / SKQDStep        .unverified_energy  .certificate  .energy (gated)
ExcitedKrylovStep            .unverified_energies  .certificate  .energies (gated)
PipelineResult               .steps  .unverified_energy  .certificate  .certificates
                             .computed_energy  .error_vs_reference  .convergence (all gated)
```

`certificate_dir` replaces `certify`. It only says where the certificate files are kept;
without it they are written to a temporary directory and removed after the check. It does
not decide whether the checker runs. `eigenstates` returns the spectrum step in place of the
bare energy list.

## 5. Acceptance criteria (validation gates)

- **G1 — the public surface is pinned.** Every public callable of the four solver classes
  and of `pipeline` is listed in the gate as an energy path or as not one, with its parameter
  names. A new public method, or a new parameter on an energy path, fails until it is
  classified. No module on the return path reads the environment.
- **G2 — DEFINITION OF DONE. Every path returns a certified energy, and it is the same
  number.** On H2, every energy path returns results whose `certificate` is a VERIFIED
  `Verdict`, whose `energy` equals `unverified_energy`, and whose certified ceiling is not
  above that energy by more than float padding (1e-7 Ha) with the exact `lambda_min` inside
  the enclosure. The lowest energy of each call has its own certificate, with the ceiling
  equal to it. The gate re-runs the checker on the certificate file behind every VERIFIED
  verdict and gets exit 0. Energies agree with the uncertified door to 1e-12 Ha.
- **G3 — an abstention cannot be mistaken for an energy.** With the checker made to refuse
  every certificate, every path returns a result whose `certificate` is not ok, whose every
  gated accessor raises `EnergyAbstained`, and in which the number is reachable only under a
  name containing `unverified`. The same holds for the door: `uncertified_eigenstates` and
  `uncertified_estimates_only=True` return abstentions that say `not submitted`, and start
  no checker process.
- **G4 — the checker's own refusal reaches the caller.** A noisy H2 solve whose energy is
  below `lambda_min` comes back abstained by the real checker, by verdict and not by crash,
  and reading its `energy` raises.
- **G5 — non-vacuity, demonstrated.** The gate's own detector is run on four mutants and
  must reject each: a result with a bare float `energy` and no certificate, an abstained
  result that also exposes the number under an innocent name, a VERIFIED verdict with no
  certificate file behind it, and a path patched to return the bare Ritz value.

## 6. Implementation plan (test-first)

1. `tests/test_certificate_or_abstention_spec.py` encoding G1–G5 (failing).
2. `hybrid_quantum_solver/certification.py`: the producer moved in, `certified`,
   `certify_energies`.
3. Result types and the twelve return paths; `pyproject.toml` (certkit a dependency).
4. Callers: tests, study scripts, `certchem`, `certkit_bridge.py`.
5. Sibling specs and gates updated; the new gate added to CI.

## 7. Out of scope

- Certifying excited levels, classical references, `estimate_energy_aer`, or the
  linear-algebra kernels (§2, item 2).
- An uncertified door on `TrotterKrylovSolver`, `HardwareKrylovSolver` and
  `SampleKrylovSolver`. No caller needs one today. A run on real hardware will: its witness
  cannot be rebuilt classically (`SPEC_certkit_return_paths.md` §2, item 5).
- Changing what `certchem`'s `Mode.FAST` returns (ADR-0004).
- A faster checker, a checker in its own environment, or a shared inertia certificate.

## 8. Caveats and risks

- **Cost.** Every single `solve` on 8 or more qubits now takes seconds where it took
  milliseconds. Sweeps are cheap (§3 c); loops of single solves are not.
- A certified energy is certified as an upper bound with a route-dependent floor. Read
  `certificate.lo` before quoting a bracket.

## 9. Deliverables

- `hybrid_quantum_solver/certification.py`, the four solver modules, `pipeline.py`.
- `tests/test_certificate_or_abstention_spec.py` — gates G1–G5.
- `certkit_bridge.py`, `pyproject.toml`, `uv.lock`, `.github/workflows/certkit-regression.yml`.
- Updated `SPEC_certkit_return_paths.md`, `SPEC_certkit_regression_gate.md`, `BACKLOG.md`.
