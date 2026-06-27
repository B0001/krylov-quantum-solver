# SPEC: single-ramp DMRG bond-dimension extrapolation

**Status:** APPROVED-BY-BACKLOG (gate pre-agreed in `BACKLOG.md`); implementing test-first.

---

## 1. Goal

The current `dmrg_energy_extrapolated` runs a **separate converged DMRG per bond dimension**
(clean truncation points, but ~`len(bond_dims)`× the sweeps). Claim: reading block2's per-stage
`get_dmrg_results()` from **one ramping DMRG run** yields the same `E(D→∞)` at a fraction of the
cost — unblocking the large-`n` Hₙ study that was too slow with the per-D protocol.

## 2. Background and honest framing

block2 reports `(bond_dims, discarded_weights, energies)` per distinct bond dimension after a single
sweep schedule. A schedule that holds each target `D` for a few sweeps gives converged-enough
per-`D` points for the discarded-weight extrapolation in one call.

**Claim:** same physics, less compute. **Not claimed:** the ramp points are *as* converged as the
per-D protocol at intermediate `D`; the gate tests whether the *extrapolated* energy still agrees.

## 3. Approach

Add `protocol="perD" | "ramp"` to `dmrg_energy_extrapolated` (default `"perD"`, unchanged). The
`"ramp"` path: one `drv.dmrg` with schedule `[D]*sweeps_per_stage for D in bond_dims`, then
`drv.get_dmrg_results()` → `(D, discarded_weight, E)` per stage → the same `E` vs `δ` linear fit.

## 4. Public interface

```
dmrg_reference.dmrg_energy_extrapolated(..., protocol="perD"|"ramp",
                                        sweeps_per_stage=4) -> ExtrapResult
```
`ExtrapResult` unchanged. Default behaviour and all existing call sites are byte-for-byte unchanged.

## 5. Acceptance criteria (validation gates)

In `tests/test_singleramp_spec.py` (block2; runs isolated via `make gates`):

- **G1 — agreement (the "done" gate).** On H₁₀ and H₁₂, `|E_ramp − E_perD| < 0.1 mHa`.
- **G2 — soundness.** On H₁₂, `|E_ramp − E_FCI| < 5e-4 Ha` (the ramp is itself accurate, not just
  close to perD).
- **G3 — speedup.** `ramp` does strictly fewer total DMRG sweeps than `perD`
  (`sweeps_per_stage * len(dims)` vs `n_sweeps_per * len(dims)`), and measured wall-time is
  `< 0.7×` perD (generous bound; the deterministic sweep-count ratio is the real proof, wall-time
  is reported).

## 6. Implementation plan (test-first)

1. `tests/test_singleramp_spec.py` with G1–G3 (failing).
2. Add the `"ramp"` branch + `sweeps_per_stage` to `dmrg_energy_extrapolated`.
3. `make gates` until green.

## 7. Out of scope

- Changing the default protocol (stays `"perD"`).
- Re-running the large-`n` Hₙ study (separate backlog item, now unblocked).

## 8. Caveats and risks

- **R1:** ramp intermediate-`D` points are less converged; if G1 fails, raise `sweeps_per_stage`
  (trading back some speedup) — the spec is met as long as G1 holds at a ratio that still beats perD.
- The ramp's `stderr` may be larger than perD's (less-converged points); acceptable.

## 9. Deliverables

- `dmrg_reference.dmrg_energy_extrapolated` — `protocol`/`sweeps_per_stage`.
- `tests/test_singleramp_spec.py` — G1–G3.
- `BACKLOG.md` — mark item done.
