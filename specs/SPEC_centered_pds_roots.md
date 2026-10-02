# SPEC: The centered moment frame removes PDS's excited-root bound violation on stretched H₄

**Status:** CLOSED, revised for chem-8tl. The revised gates G1–G3 (+G1b) PASS, 4/4, on macOS 27
arm64 / Accelerate (2026-10-02). **The claim survives** on all 33 paths measured there. These gates
have not been run on Linux in this form. The Linux/OpenBLAS numbers in §6 come from the chem-8tl
authoring run (commit 1ff123a), whose G1 pinned shuffle 15 and had no term-order paths.

**Original (2026-10-01, chem-dcz):** gates G1–G3 (+G1b) PASS. The raw −mHa violation reproduces
deterministically. Its "bimodality" was never BLAS threading. It is **`PYTHONHASHSEED`**:
hash-randomized iteration inside the qiskit-nature operator build changes the Pauli-term order and
the last bits of the coefficients of the qubit Hamiltonian. That changes the last bits of the moments,
and at cond(M) ≈ 1.2e17 those bits pick the path. On both paths the centered roots stay above every
reachable target for K=3..8, and they agree to < 1e-13 Ha between the paths. **Amendment:** this
machine's violation is **−3.548 mHa**, not `SPEC_pds_excited_roots`' −7.49 mHa. The sign and the
cond(M) match; the size depends on LAPACK.

**Revision (2026-10-02, chem-8tl):** `build_molecular_hamiltonian` is now hash-seed independent
(`map_canonical`), so `PYTHONHASHSEED` no longer selects a path, and the old G1 and G1b were killed by
design. The mechanism underneath turned out to be narrower than "hashing". The fermionic coefficients
were always bit-identical across seeds. Only the operator's key order varied. The Jordan–Wigner mapper
sums colliding Pauli terms in that order, which sets the low coefficient bits, and the qubit operator
inherits it as its term order, which sets the order the moments accumulate in. The child now builds
the paths explicitly: the canonical build, the same fermionic terms mapped in 16 **seeded shuffled
orders**, and the canonical operator with only its **term order** permuted (16 seeds). The claim
survives on both platforms: on Linux over the authoring run's 17 paths, and on macOS over all 33.
Some paths violate raw: shuffles 4 and 15 on Linux/OpenBLAS (−1.33 and −5.91 mHa at K=7); on
macOS/Accelerate 10 of 16 shuffles and 10 of 16 term orders, by > 1 mHa (worst: shuffle 7,
−12.00 mHa). On each platform the centered margins are the same on every path, to < 1e-8 Ha, and none
falls below its target.

Two findings came out of re-running it on macOS. (1) **The pinned shuffle did not transfer.**
Shuffle 15 was chosen because it was the worst path on Linux; on macOS it gives **+0.0549 mHa**, so
the first revision's G1 failed there. This is R1 happening. G1 now takes the worst of the
pre-registered paths, which is what §1(a) claims (some pinned configuration reproduces the
violation). The thresholds are unchanged. (2) **Term order alone flips raw.** With the coefficient
bits fixed, permuting only the qubit term order moves raw K=7 between **−11.83 and +0.0097 mHa** on
macOS. So the hidden variable is the operator's bits *and* its term order, not summation order
alone, and `map_canonical` has to fix both. The canonical build itself respects the bound on both
platforms (+0.1595 mHa Linux, +0.0068 mHa macOS), but that is one draw of the bits, not a property.

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
  the centered margins are the same on every bit-path.
- **What we cannot claim.** (a) That centering guarantees the bound in general. Centered cond(M) is
  7.2e7 at K=8 here, but `SPEC_centered_pds` R2 already shows LiH reaching 8.6e12; the next wall
  exists. (b) That the violation's *size* is portable, since it depends on LAPACK. (c) Anything
  past K=8, or for other systems.

## 3. Approach

No library change. Each configuration runs in a child process (the test file as `__main__`),
with `PYTHONHASHSEED` and the thread count fixed before the interpreter starts. Since chem-8tl the
child builds 33 bit-paths: the canonical Hamiltonian (`"canonical"`); the same sorted fermionic terms
mapped in `np.random.default_rng(s).permutation` order (`"0"`..`"15"`, new coefficient bits and term
order); and the canonical qubit operator with only its terms permuted by
`np.random.default_rng(s).permutation` (`"order0"`..`"order15"`, same coefficient bits). The child
computes raw and `centered_frame` moments with `hamiltonian_moments`, `pds_roots` for K=3..8, the
per-K min_j(root_j − E_j), and cond(M).

**Pinned configuration (revised):** all 33 paths at `PYTHONHASHSEED=2` (any seed gives the same bits
now) and 1 thread. G1 takes the worst path at runtime, because which one that is depends on the
platform (§6). Measured on macOS 27.0.1 arm64, Python 3.13.14, numpy 2.4.6 (Accelerate), scipy
1.15.3, qiskit 2.5.0, qiskit-nature 0.8.0, pyscf 2.13.1; the Linux columns are from Linux x86_64,
Python 3.13, numpy/OpenBLAS. **Original (pre-chem-8tl):** `PYTHONHASHSEED=2` (violating) / `0` (clean);
`OMP/OPENBLAS/MKL/VECLIB_MAXIMUM_THREADS=1`; macOS 27.0.1 arm64, Python 3.13.14, numpy 2.4.6
(Accelerate), scipy 1.15.3, qiskit 2.5.0, qiskit-nature 0.8.0, pyscf 2.13.1; geometry and basis as
above; K=3..8; no RNG involved.

## 4. Public interface

None new.

## 5. Acceptance criteria

Gates in `tests/test_centered_pds_roots_spec.py`.

- **G1: raw violation reproduced, pinned.** Hash seed 2, 1 thread: on the worst of the 33 paths, raw
  PDS(7) min margin < −1e-3 Ha and raw cond(M₇) > 1e14. A second fresh process gives a bit-identical
  result on every path. *Superseded:* "shuffle 15 violates". It held on Linux, where shuffle 15 was
  measured and then pinned, but not on macOS (+0.0549 mHa).
- **G1b: hidden variable (revised).** Hash seeds 0 and 2 give identical results on every path
  (chem-8tl). Across the paths the raw moment bits differ, and raw PDS(7) both violates
  (< −1e-3) and respects the bound (≥ −1e-9). *Killed original:* "seed 0 vs seed 2 differ", and "4
  threads = 1 thread", which held on Accelerate but not on Linux OpenBLAS, where the multi-threaded
  SCF is not bit-reproducible run to run.
- **G2: the claim.** On all 33 paths, every centered root_j ≥ E_j − 1e-9 at K=3..8, and the centered
  margins agree across paths to < 1e-8 Ha.
- **G3: why.** At K=7, raw cond(M) > 1e14 and centered cond(M) < 1e10 on every path.

## 6. Measured data (min_j root_j − E_j, mHa)

**Revised (chem-8tl), Linux x86_64 / OpenBLAS, 1 thread.** From the authoring run (commit 1ff123a,
17 paths, no term-order paths); not re-run for this revision:

| K | raw, canonical | raw, shuffle 15 | centered (all 17 paths) | cond raw (canon / s15) | cond centered |
|---|---|---|---|---|---|
| 3 | +2.4146 | +2.4146 | +2.4146 | 2.3e5 | 93 |
| 4 | +1.0780 | +1.0780 | +1.0780 | 9.9e8 | 6.1e3 |
| 5 | +0.8624 | +0.8624 | +0.8624 | 4.5e11 | 4.0e4 |
| 6 | +0.6127 | +0.6123 | +0.6123 | 1.5e14 | 2.1e5 |
| 7 | +0.1595 | **−5.9075** | +0.1198 | 1.2e17 | 2.4e6 |
| 8 | +0.1031 | +0.1115 | +0.0042 | 5.4e19 / 2.1e19 | 7.2e7 |

Raw K=7 across the 16 shuffles: −5.9075 (s15), −1.3318 (s4), then +0.0271 to +0.1708. *Removed:* a
sentence here gave a term-order-only range on Linux (+0.029 to +0.178 mHa), but no code in the repo
produced it. The child now measures that experiment (`"order*"` paths); its Linux numbers are a
follow-up.

**Revised (chem-8tl), macOS 27 arm64 / Accelerate, 1 thread, all 33 paths** (§3 configuration;
regenerate with `PYTHONHASHSEED=2 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 uv run python
tests/test_centered_pds_roots_spec.py`, which prints the child's JSON):

| K | raw, canonical | raw, shuffle 7 (worst) | centered (all 33 paths) | cond raw (canon / s7) | cond centered |
|---|---|---|---|---|---|
| 3 | +2.4146 | +2.4146 | +2.4146 | 2.3e5 | 93 |
| 4 | +1.0780 | +1.0780 | +1.0780 | 9.9e8 | 6.1e3 |
| 5 | +0.8624 | +0.8624 | +0.8624 | 4.5e11 | 4.0e4 |
| 6 | +0.6122 | +0.6118 | +0.6123 | 1.5e14 | 2.1e5 |
| 7 | +0.0068 | **−11.9992** | +0.1198 | 1.2e17 / 1.3e17 | 2.4e6 |
| 8 | +0.1045 | +0.0997 | +0.0042 | 3.4e19 / 1.6e19 | 7.2e7 |

Raw K=7 across the 16 shuffles: 10 violate by more than 1 mHa (shuffles 0, 2–9, 13; −3.76 to −12.00
mHa), 3 violate by less (10, 11, 12: −0.53, −0.98, −0.50 mHa), and 3 respect the bound (1, 14, 15:
+0.080, +0.103, +0.055 mHa). Across the 16 term-order-only permutations: 10 violate by more than 1 mHa
(worst `order14`, −11.83 mHa), 5 by less, and 1 respects it (`order1`, +0.0097 mHa). The centered
margins agree across all 33 paths to 5.7e-14 Ha. Under the old library build (`mapper.map` of
`second_q_op()` without `map_canonical`), the original gates at commit 067469e still pass on this
machine (4/4, 2026-10-02): seed 2 violates and seed 0 respects.

**Original (pre-chem-8tl, macOS/Accelerate):**

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

- ~~Making `build_molecular_hamiltonian` hash-seed deterministic~~: done in chem-8tl
  (`molecular_hamiltonian.map_canonical`, gate `tests/test_hamiltonian_determinism.py`).
- Changing `moment_expansion` to center by default. Nothing here justifies an API change.

## 8. Caveats and risks

- **R1: the path map is platform-specific, and this has now been seen.** Shuffle 15 is the worst path on
  Linux/OpenBLAS and respects the bound on macOS/Accelerate (+0.0549 mHa), where shuffle 7 is worst.
  Two violating Linux shuffles out of 16 became 10 on macOS. This is why G1 takes the worst
  pre-registered path instead of a named one. If G1 fails somewhere, that means no path among the 33
  violates by 1 mHa there, and §1's claim cannot be tested on that platform; widen the
  pre-registered set in the spec before re-running, never after seeing the result. The general
  mechanism, bits selecting the path at cond > 1e17, is what transfers.
- **R2: one geometry.** Stretched H₄ only.
