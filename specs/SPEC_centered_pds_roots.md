# SPEC: The centered moment frame removes PDS's excited-root bound violation on stretched H₄

**Status:** CLOSED. Gates G1–G3 (+G1b) PASS (2026-10-01, chem-dcz). **The claim survives**, now on a
pinned configuration. The raw −mHa violation reproduces deterministically. Its "bimodality" was
never BLAS threading. It is **`PYTHONHASHSEED`**: hash-randomized iteration inside the
qiskit-nature operator build changes the Pauli-term order and the last bits of the coefficients of
the qubit Hamiltonian. That changes the last bits of the moments, and at cond(M) ≈ 1.2e17 those bits
pick the path. On both paths the centered roots stay above every reachable target for K=3..8, and
they agree to < 1e-13 Ha between the paths. **Amendment:** this machine's violation is
**−3.548 mHa**, not `SPEC_pds_excited_roots`' −7.49 mHa. The sign and the cond(M) match; the size
depends on LAPACK.

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

`SPEC_centered_pds` §7 left a lead: "centering fixes `SPEC_pds_excited_roots` R3's root-bound
violation". It was unproven because the raw violation did not reproduce there, and an earlier attempt
never pinned it. Claim: (a) the raw K=7 violation reproduces in one fully pinned configuration, and
(b) on that configuration, every centered root_j ≥ reachable E_j − 1e-9 Ha for K=3..8. **Dies if**
no pinned configuration reproduces the raw violation (then "fixes it" cannot be tested, so the spec is
revised), **or** any centered root falls below its target on a configuration where raw violates.

## 2. Background and honest framing

- **Reference.** The reachable exact spectrum (dense eigh, |⟨v|HF⟩|² > 1e-8), built the same way as
  in `SPEC_pds_excited_roots`. Stretched H₄ (0/2/4/6 Å, sto3g, 8 qubits).
- **What we can claim.** On this geometry and configuration, the raw-frame bound violation is a
  deterministic function of the moment bits. Centering removes it and also removes the sensitivity:
  the centered margins are the same on both bit-paths.
- **What we cannot claim.** (a) That centering guarantees the bound in general. Centered cond(M) is
  7.2e7 at K=8 here, but `SPEC_centered_pds` R2 already shows LiH reaching 8.6e12; the next wall
  exists. (b) That the violation's *size* is portable, since it depends on LAPACK. (c) Anything
  past K=8, or for other systems.

## 3. Approach

No library change. Each configuration runs in a child process (the test file as `__main__`).
`PYTHONHASHSEED` has to be fixed before the interpreter starts, so it cannot be set from inside the
test process. The child computes raw and `centered_frame` moments with `hamiltonian_moments`,
`pds_roots` for K=3..8, the per-K min_j(root_j − E_j), and cond(M).

**Pinned configuration (recorded):** `PYTHONHASHSEED=2` (violating) / `0` (clean);
`OMP/OPENBLAS/MKL/VECLIB_MAXIMUM_THREADS=1`; macOS 27.0.1 arm64, Python 3.13.14, numpy 2.4.6
(Accelerate), scipy 1.15.3, qiskit 2.5.0, qiskit-nature 0.8.0, pyscf 2.13.1; geometry and basis as
above; K=3..8; no RNG involved.

## 4. Public interface

None new.

## 5. Acceptance criteria

Gates in `tests/test_centered_pds_roots_spec.py`.

- **G1: raw violation reproduced, pinned.** Seed 2, 1 thread: raw PDS(7) min margin < −1e-3 Ha, and
  raw cond(M₇) > 1e14. A second fresh process gives a bit-identical result.
- **G1b: hidden variable.** Seed 0 produces different raw moment bits, and its raw PDS(7) respects the
  bound. Seed 2 at 4 threads equals seed 2 at 1 thread exactly.
- **G2: the claim.** On seeds 0 and 2, every centered root_j ≥ E_j − 1e-9 at K=3..8, and the centered
  margins agree across the seeds to < 1e-8 Ha.
- **G3: why.** At K=7, raw cond(M) > 1e14 and centered cond(M) < 1e10 on both seeds.

## 6. Measured data (min_j root_j − E_j, mHa)

| K | raw, seed 0 | raw, seed 2 | centered (both) | cond raw | cond centered |
|---|---|---|---|---|---|
| 3 | +2.4146 | +2.4146 | +2.4146 | 2.3e5 | 93 |
| 4 | +1.0780 | +1.0780 | +1.0780 | 9.9e8 | 6.1e3 |
| 5 | +0.8624 | +0.8624 | +0.8624 | 4.5e11 | 4.0e4 |
| 6 | +0.6120 | +0.6119 | +0.6123 | 1.5e14 | 2.1e5 |
| 7 | +0.0259 | **−3.5478** | +0.1198 | 1.2e17 | 2.4e6 |
| 8 | +0.0926 | +0.0987 | +0.0042 | 4.4e19 (s0) / 3.3e19 (s2) | 7.2e7 |

Seed scan 0..15, 1 thread: seeds {2, 4, 5, 6, 8, 10, 14} give moment hash `b680653a…` and −3.548 mHa
at K=7. The other 9 seeds give `6986dfc5…` and +0.026 mHa. Thread count (1 vs 4) never changed
a result for a fixed seed. Across the two seeds the sorted Pauli label set is identical (185 terms),
but both the term order and the coefficient bits differ. Raw K=6 already departs from centered by
0.4 µHa, and raw K=7 on the clean path is 0.09 mHa looser than centered.

## 7. Out of scope

- Making `build_molecular_hamiltonian` hash-seed deterministic (canonical term order and summation).
  That would change bits repo-wide and could move other gates that sit near singular matrices. It
  should be its own bead.
- Changing `moment_expansion` to center by default. Nothing here justifies an API change.

## 8. Caveats and risks

- **R1: the seed→path map is platform-specific.** The seed-2 path was measured on the configuration
  in §3. A different qiskit-nature or Python version could reshuffle which seeds violate, and if G1
  fails elsewhere, scan seeds before concluding anything. The general mechanism, bits selecting the
  path at cond > 1e17, is what transfers.
- **R2: one geometry.** Stretched H₄ only.
