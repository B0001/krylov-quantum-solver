# SPEC: Every emitted certificate must pass the independent checker, in CI

**Status:** IMPLEMENTED — gates G1–G5 green (`tests/test_certkit_regression_gate_spec.py`,
`.github/workflows/certkit-regression.yml`).

---

## 1. Goal

`certkit_bridge.py` turns a QKSD Krylov solve into certificates and hands them to an
independent checker. Nothing enforced that the answer stayed *yes*. Claim: every certificate
this repo's producer emits still verifies under the exact certkit release pinned in
`pyproject.toml`, and the energy each verified certificate certifies is still within chemical
accuracy of the exact reachable minimum — checked automatically, on every push, by a checker
that never shares a process with the solver. Falsifiable at machine precision: one certificate
that stops verifying, one route that silently stops being emitted, or one certified upper
bound that drifts past 1.6 mHa fails the gate.

## 2. Background and honest framing

- The sibling gate `certkit-conformance.yml` proves the *pinned checker* still passes
  certkit's own conformance suite. It deliberately does not run the producer. This spec is
  the other half: the checker is known-good, so what does it say about *our* certificates.
- **What we can claim if the gates pass:** the Krylov ground-state energies this repo
  produces carry enclosures that an independent, out-of-process checker re-derives from the
  certificate and the operator alone, at a pinned release, with the exact energy inside every
  one of them and the certified upper bound within chemical accuracy.
- **What we cannot claim.** Three limits, all load-bearing:
  1. **The checker alone proves almost nothing about the solver.** The `gershgorin_rayleigh`
     route certifies `[gershgorin_lower(H), ⟨x|H|x⟩]`, a true enclosure for *any* unit vector
     `x`. A witness of pure random noise, 0.84 Ha above the ground state, is VERIFIED. G5
     demonstrates this rather than assuming it. Soundness comes from the checker; correctness
     comes from G3/G4 comparing against a reference. Neither is sufficient alone.
  2. **Coverage is the producer's four cases, not the test suite.** `certkit_bridge.py` emits
     for H2, H4, H4-stretched and N2. No regression test in `tests/` emits a certificate,
     because no solver return path produces one — that is `portfolio-cir` (emit for every
     solved ground-state energy) and `portfolio-0pc` (certificate-or-abstention as a
     non-bypassable invariant), both still open. This gate is written so those land *into* it:
     new cases join `EXPECTED` and are gated on arrival.
  3. **Electronic frame.** Certificates are about the qubit matrix; `mh.energy_offset` is a
     trusted addition made outside it. The pinned references here are therefore not the total
     energies the rest of the suite asserts.
- The checker is consumed as a protocol, never imported: `certkit/INTEGRATION.md`.

## 3. Approach

`certkit_bridge.run_case(case, out)` emits every certificate a case supports and checks each
one by running `python -m certkit.cli check <certificate> <operator>` as a separate process,
returning the exact electronic-frame `lambda_min` and one `Verdict` per certificate. The gate
pins, per case, the reference energy and the exact set of certificates with each one's
expected status, then asserts the four properties below. CI installs the `certkit` extra —
the exact git tag from `pyproject.toml` — and runs this gate file.

**Reference:** exact `lambda_min` of the qubit Hamiltonian by dense diagonalization
(`mh.ground_state_energy() - mh.energy_offset`), itself cross-checked against independent
PySCF FCI in `tests/test_reference_energies.py`. Pinned as a literal *and* re-derived at test
time, so neither the pin nor the derivation can drift alone.

## 4. Public interface

```
certkit_bridge.Verdict                                   -> (name, rule, ok, line, lo, hi)
certkit_bridge.check_certificate(cert_path, op_path)     -> (ok, first stdout line)
certkit_bridge.emit(enc, rule, x, lo, hi, out, name, beta=None) -> Verdict
certkit_bridge.run_case(case, out)                       -> (lambda_min, [Verdict])
.github/workflows/certkit-regression.yml                 -> the CI gate
```

Certificate filenames are now per route (`certificate_sector`, `certificate_temple`,
`certificate_gershgorin`). The old mode-dependent `certificate.json`, which meant
`temple_inertia` or `gershgorin_rayleigh` depending on whether the dense route succeeded,
cannot be pinned by a gate and is gone.

## 5. Acceptance criteria (validation gates)

- **G1 — verdicts match the pin.** For every emitted certificate, the independent checker's
  status equals the expected status, *and the checker's own output line agrees with its exit
  code*. That second half matters: exit 1 covers both ABSTAIN and a crash, so without it a
  checker that started throwing on the sector certificates would be accepted as the ABSTAIN
  they are pinned to. Expected ABSTAINs are pinned as ABSTAIN: the two sector certificates
  take `eps` from the solver's own second Ritz pair, which the checker cannot discharge by
  inertia count. Only the status is pinned, not the reason — the reason varies with Pauli
  term ordering, which is not stable across processes.
- **G2 — the emitted set is exactly the pinned set.** A route that silently stops emitting is
  a loss of coverage, not a pass.
- **G3 — soundness. DEFINITION OF DONE.** The exact `lambda_min` lies inside every VERIFIED
  enclosure, and the re-derived reference agrees with the pinned literal to `< 1e-6` Ha. One
  escape kills the claim.
- **G4 — the enclosure is still tight at both ends.** For every VERIFIED certificate,
  `hi - lambda_min <= 1.6e-3` Ha (chemical accuracy), *and* `lambda_min - lo` stays under a
  per-certificate pin. Both halves are needed, and for different reasons. `hi` is the
  solver's own energy, so gating it is what makes this a *regression* gate rather than a
  soundness check. `lo` needs its own pin because the checker accepts **any** lower bound
  beneath the one it re-derives: a certificate claiming `lambda_min > -1e6` Ha is VERIFIED and
  worthless, which would otherwise leave N2 — whose only route is Gershgorin, and whose entire
  content is a 2.6 Ha-wide floor — gated on nothing.

  Calibration, so the headroom is not mistaken for more than it is. Today's worst upper-bound
  margin is 1.6e-4 Ha (H4-stretched), a decade below the threshold. But `hi` is a Rayleigh
  quotient and therefore *quadratically* insensitive to witness error, so in eigenvector terms
  1.6 mHa is a 2–4% wrong witness depending on the case (N2 is tightest at 2.3%). This gate
  catches a broken solver, not a slightly degraded one.
- **G5 — non-vacuity, demonstrated.** A `gershgorin_rayleigh` certificate built on a random
  noise witness is VERIFIED by the checker *and* is caught by G4's criterion. If the checker
  ever gets strong enough to reject it, this gate fails and says so — that is a finding about
  certkit, not a regression here.

Not a gate: certificate *width* as such. H2's certified width is 3.7e-9 Ha, which is entirely
`pad_claim`'s floating-point padding at the hardcoded `rel=1e-9`; a width threshold there would
measure the pad, not the solver. G4's two ends are pinned separately instead, which is the same
information where it means something and avoids the pad where it does not.

## 6. Implementation plan (test-first)

1. `tests/test_certkit_regression_gate_spec.py` encoding G1–G5.
2. Minimum producer changes to make them expressible: per-route certificate filenames,
   `run_case` returning verdicts instead of printing them and returning a disjunction, and no
   certificate at all where Temple has no finite bound (`-inf` is not representable, and the
   honest translation of "vacuous" is abstention).
3. `.github/workflows/certkit-regression.yml` runs the gate on ubuntu-latest.

## 7. Out of scope

- Emitting certificates from the solver's return path, or from existing regression tests
  (`portfolio-cir`, `portfolio-0pc`).
- Running the checker in a *separate environment* from the producer. It already runs
  out of process from a pinned release; a dedicated checker venv — certkit has zero runtime
  dependencies, so it is cheap — would be strictly stronger and is a good follow-up.
- Claim kinds certkit does not have: gaps, dipoles, relative energies, overlaps. The
  `certified_*` gates assert quantities no single certificate can express.
- Tying the encoded operator to the matrix the reference comes from. `encode_hamiltonian` is
  a second, independent translation of `mh.qubit_hamiltonian`, and the checker re-derives
  everything from that encoding — so a broken encoding is checked against itself. On the dense
  cases the Temple route pins it tightly by accident (a sign flip, a 0.999 coefficient scale
  and dropping the LSB/MSB reversal all turn G1 red), and G4's pinned floor catches a shifted
  encoding on the Gershgorin-only cases. Making it airtight means having `run_case` return
  `<x|decode_operator(enc)|x>` and asserting it equals the certified `hi`; worth doing when a
  second Gershgorin-only case appears.
- Byte-reproducible certificates. The producer's Pauli term ordering varies with
  `PYTHONHASHSEED`, so certificate and operator are only valid as a same-run pair. The gate
  regenerates both together, which sidesteps it; committing golden certificates would not.
