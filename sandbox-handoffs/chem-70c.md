# chem-70c handoff: is PDS's "det M → 0 at large K" floor a frame artifact?

## Verdict: the claim SURVIVES, with one amendment

Neither kill condition fired:

| kill condition | result |
|---|---|
| centered ≠ raw by > 1e-6 Ha where raw cond(M) < 1e10 | **max diff 8.2e-12 Ha** (H₄ K=4); H₂ 0, LiH 2.3e-12, N₂ 3.6e-12 |
| < 8 orders of cond(M) reduction at K=8 | **13.0 (H₄), 11.3 (LiH), 17.4 (N₂ CAS(6,6))** orders |

**Amendment (the finding):** `SPEC_moment_pds` R1 mixes two causes.
1. The large-K ill-conditioning is mostly the **frame**: raw moments grow like ‖H‖ⁿ. Centering
   removes it and gives K=7–8.
2. "Reference → exact state" is **fundamental**. Once K exceeds the number of HF-reachable
   eigenstates, M is exactly singular in every frame. H₂ (2 reachable states): centering cuts only
   **2.3 orders** at K=8, and centered cond(M) is already 3.8e16 at K=3. The claim's 8-order bar
   *fails* on H₂. This is harmless (PDS is exact at K=2 there) but real, so it is gated as the
   boundary (G2b) rather than hidden.

## N₂ CAS(6,6) non-monotonicity: reproduced, and worse than the scout reported

Raw N₂ PDS at K≥6 has cond(M) of 4.7e20 to 2.3e26, which is past 1/eps. The outcome depends on the
BLAS floating-point path. Errors are in mHa, K=6/7/8, from six fresh processes:

| BLAS threads | K=6 | K=7 | K=8 | failure |
|---|---|---|---|---|
| 1 | 0.0340 | 0.0464 | 0.0287 | non-monotone 6→7 (matches the scout's 0.044→0.049) |
| 2 | 0.0416 | 0.0477 | 0.0262 | non-monotone 6→7 |
| 4 | 0.0564 | 0.0378 | −337 045 | non-variational, 337 Ha below FCI |
| default | 0.0572 | −18 754 | 0.0315 | non-variational |
| default | 0.0686 | 0.0371 | −23 893 | non-variational |
| default | 0.0581 | 0.0417 | −977.7 | non-variational |

All six raw runs violated `SPEC_moment_pds` G2 or G3 somewhere in K=6..8. **Centered results were
identical in every run**, at 0.0540 / 0.00587 / 0.00106 mHa: monotone, variational, and matching
the scout's centered numbers. The fix is to build the moments in `centered_frame`. The raw failure
mode is path-dependent, so G4 gates only the robust part: cond > 1e17, and the deviation from
centered exceeds 1e-5 Ha. The table is recorded in the spec.

## μ sensitivity (μ from dense diagonalization only works at validation scale)

- N₂, μ = μ_center + δ: for δ ∈ [−5, +2] Ha, PDS(6..8) moves by ≤ 0.07 µHa. The reachable
  half-width is 2.39 Ha. At δ = +5, PDS(8) dips 1.8 µHa below FCI. At δ = +10 it breaks (PDS(7) is
  311 mHa below FCI). The response is asymmetric: shifting up past the band hurts.
- **μ = ⟨H⟩ (the first moment, which is free and needs no diagonalization)** reproduces centered PDS to
  < 1e-6 Ha at every K ≤ 8 on H₄/LiH/N₂, and still cuts cond(M) by 9.5 / 12.2 / 14.1 orders at K=8.
  So the validation-scale μ is not what makes this work.

## What changed

- `specs/SPEC_centered_pds.md` (new): the spec, gates, the full data table and caveats.
- `tests/test_centered_pds_spec.py` (new): G1, G2, G2b, G3 (**the K=7–8 gates**: centered PDS is
  variational and monotone for K=1..8; PDS(8) error < 1e-6 Ha on H₄, < 5e-6 Ha on N₂, < 1.6 mHa on
  LiH), G4 (raw breakdown on N₂), and G5 (μ sensitivity and μ = ⟨H⟩). BLAS is pinned to 1 thread,
  as in `test_pds_excited_roots_spec.py`.
- `specs/SPEC_moment_pds.md`: R1 and §7 annotated with the split cause and a pointer to the new spec.
  The K ≤ 4 raw-frame gates are unchanged.
- `specs/BACKLOG.md`: the chem-70c entry is marked CLOSED with numbers.
- **No library code changed.** It is pure composition of `hamiltonian_moments`, `pds_energy` and
  `device_odmd.centered_frame`. The test's `_shifted` helper is the same `dataclasses.replace`
  shift that `centered_frame` does, at an arbitrary μ (used for G5 only).

## Validation

```
uv run pytest -q tests/test_centered_pds_spec.py tests/test_moment_pds_spec.py \
    tests/test_pds_excited_roots_spec.py
15 passed in 397s
uv run ruff check tests/test_centered_pds_spec.py   # All checks passed
```

I didn't run the full `make gates` / `make test` (no library code changed). The new gate takes about 3–4 min,
mostly from the 4096×4096 dense eigh for N₂, which is called twice.

## Honest caveats

- This is numerical hygiene, not physics. The convergence rate is unchanged.
- LiH's centered cond(M) at K=8 is already 8.6e12, so expect the next float64 wall around K≈10–11.
- Exact statevector moments only. Shot-noise propagation in ⟨(H−μ)ⁿ⟩ versus ⟨Hⁿ⟩ is untested.

## Follow-up leads (not filed: no bd writes per instructions)

1. **Excited roots (chem-pc1):** on stretched H₄, centered `pds_roots` stayed above every reachable
   target at K=3..8, at both 1 and 4 threads (K=8 margin +0.0042 mHa versus +0.09 for raw). But the
   raw run here did not reproduce that spec's −7.49 mHa violation at K=7, so "centering fixes it" is
   unproven. This needs the bimodal raw path in hand first.
2. Consider having `moment_expansion` shift by ⟨H⟩ internally. G5 shows it's free and nearly as good
   as the exact center. That would change a closed spec's API, so it needs its own spec.
3. Shot-noise study: does centering change how ⟨Hⁿ⟩ measurement noise amplifies into PDS(K)?
