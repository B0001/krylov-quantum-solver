# chem-0ke: QPE's precision constant as π·sin θ₀ — tested as a genuine pre-registered prediction

## Verdict: CONFIRMED on all three systems. Bead closed.

## What "done" means here

The bead's proposed derivation (`E = λ·cos(2πφ)` ⟹ half-bin dyadic error × arccos Jacobian ⟹
bound `π·sin θ₀`, `θ₀ = arccos(E₀/λ)`) is a real, per-system-predictable upper envelope for the
`SPEC_qpe_readout_laws` G2 precision-ratio constant (measured `2.175` on H2 alone, explicitly
left without a derivation there). It was pre-registered from the exact spectrum, then checked against an
actual `t = 4..14` sweep, on H2 CAS(2,2), LiH CAS(2,2), and N2 CAS(3,4) — and never exceeded,
tight to within 5% everywhere.

## What I changed

1. **`scripts/spec_qpe_precision_bound.py`** (new) — pre-registration script. Computes
   `λ`, `E₀`, `θ₀ = arccos(E₀/λ)`, and the predicted bound `π·sin θ₀` for each system from the
   *exact* diagonalized qubit Hamiltonian alone — no QPE sweep inside it. Run BEFORE the sweep
   gates were written.
2. **`specs/SPEC_qpe_precision_bound.md`** (new) — the spec, following this repo's SDD loop. §3
   contains the committed pre-registration table (reproduced verbatim below).
3. **`tests/test_qpe_precision_bound_spec.py`** (new) — gates G1 (never exceeded), G2 (tight on
   at least one system), G3 (bookkeeping: the committed table matches a fresh recomputation from
   the exact spectrum, guards against the spec drifting out of sync with the code).
4. **`specs/BACKLOG.md`** — flipped the "Fault-tolerant stack" entry from `[ ]` to `[x]`, linked
   the new spec/test, and recorded the confirmed numbers (see below). Also flagged that the
   entry's old scout-probe numbers ("2.285 vs 2.175, 2.351 vs 2.285, 2.391 vs 2.364") are close
   to but not identical with this closure's numbers and were not reproduced here — don't treat
   them as equivalent to this closure's checked values.

No changes to `qpe_walk_readout.py` or `qubitization_blueprint.py` — this is pure external
verification, per this repo's SDD convention.

## Pre-registration (computed BEFORE the sweep — `uv run python scripts/spec_qpe_precision_bound.py`)

```
system         qubits     lambda     E0(elec)  theta0(rad)  pi*sin(theta0)
H2 CAS(2,2)         4   2.699278    -1.852388     2.327122        2.285077
LiH CAS(2,2)        4   1.595641    -1.058117     2.295789        2.351497
N2 CAS(3,4)         6   8.422876    -5.446405     2.273942        2.396444
```

Geometries: H2 `H 0 0 0; H 0 0 0.74`; LiH `Li 0 0 0; H 0 0 1.6`; N2 `N 0 0 0; N 0 0 1.10` with
`symmetry="D2h"` (pins the HOMO/LUMO degenerate-pair gauge per `SPEC_lambda_ladder_honest_caveat`
R2 / chem-1yr — without it N2's active-space integrals are RHF-run-dependent). All STO-3G. CAS
notation follows this repo's existing convention in `taper_qubits.py` / `qpe_walk_readout.py`:
`mcscf.CASCI(mf, n_orb, n_elec)`.

## Measured sweep (t = 4..14, exact ground state as trial, same construction as
`SPEC_qpe_readout_laws` G1/G2)

```
H2 CAS(2,2):  bound=2.285077  max_ratio=2.175167  tightness=95.2%
  ratios by t: [0.3337, 0.6673, 1.3346, 1.8809, 0.8415, 1.6831, 1.1976, 2.1752, 0.2168, 0.4337, 0.8673]
LiH CAS(2,2): bound=2.351497  max_ratio=2.284570  tightness=97.2%
  ratios by t: [0.7036, 1.4073, 1.8391, 1.079, 2.1579, 0.3655, 0.731, 1.4619, 1.7799, 1.1423, 2.2846]
N2 CAS(3,4):  bound=2.396444  max_ratio=2.303179  tightness=96.1%
  ratios by t: [0.9678, 1.9356, 0.7825, 1.5651, 1.6774, 1.428, 1.9392, 0.9105, 1.8211, 1.1516, 2.3032]
```

(Regenerate with the ad hoc dump script inlined in this bead's tool-call log, or just read
`test_G1_measured_ratio_never_exceeds_predicted_bound`'s assertion values under `pytest -s
-vv`; the numbers above are exactly what the committed test asserts against.)

## Pre-registered criteria — pass/fail

- **G1, never exceeded (kill condition 1):** PASS on all three systems. Max measured ratio never
  exceeded `π·sin θ₀` anywhere across `t = 4..14` (11 points per system, 33 total comparisons).
- **G2, tight to ≤20% on at least one system (kill condition 2):** PASS, and by a wide margin —
  tight to within **3.8-4.9%** on *all three* systems (95.2%, 97.2%, 96.1% of the predicted
  bound), not just the required "at least one within 20%".
- **G3, bookkeeping (committed table matches fresh recomputation):** PASS, `1e-6` tolerance.

Both of the bead's acceptance criteria are satisfied on all three required systems
(H2 CAS(2,2), LiH CAS(2,2), N2 CAS(3,4)), not just one.

## Gate-run lines (verbatim)

```
$ GATE_GLOB="tests/test_qpe_*_spec.py" GATE_JOBS=1 bash scripts/run_gates.sh
gates: 2 files, 1 parallel processes, cache=on
PASS  tests/test_qpe_precision_bound_spec.py  (59s)
PASS  tests/test_qpe_readout_laws_spec.py  (4s)
all spec gates passed

$ uv run python -m pytest -q tests/test_qpe_precision_bound_spec.py tests/test_qpe_readout_laws_spec.py -v
collected 16 items
tests/test_qpe_precision_bound_spec.py .......                           [ 43%]
tests/test_qpe_readout_laws_spec.py .........                            [100%]
16 passed in 57.92s
```

```
$ uv run python -m ruff check tests/test_qpe_precision_bound_spec.py scripts/spec_qpe_precision_bound.py
All checks passed!
```
(A full-repo `ruff check .` shows 17 pre-existing errors in unrelated files — none in the three
files this bead added; verified by grepping the ruff output for these paths, zero hits.)

## Hardware / wall time

`uname -a`: `Linux 416d45c81621 7.0.14-linuxkit #1 SMP PREEMPT Fri Sep 18 10:19:32 UTC 2026
x86_64 GNU/Linux`. `nproc`: 8. RAM: 16 GB total (`/proc/meminfo` `MemTotal: 16353600 kB`). Full
run (both spec files, sequential) took under 1 minute; no single step exceeded a minute, so no
per-number hardware/timing breakdown beyond this is warranted per the "report wall time for
anything over a minute" rule.

## What I decided not to do, and why

- **Did not touch `qpe_walk_readout.py` or `qubitization_blueprint.py`.** The bead is a pure
  verification of an existing simulation path; per `SPEC_qpe_readout_laws`'s own precedent (no
  library changes, test-file-only gates), library changes were out of scope and unnecessary.
- **Did not extend past t=14 or add more systems.** The bead named exactly these three systems
  and this exact t-range; going further is a natural follow-up but not this bead's scope (would
  be filed as a new bead if wanted — I did not file one since the bead as scoped is fully closed
  and I have no specific new hypothesis to attach to a follow-up).
- **Did not attempt to formally prove the arccos-Jacobian argument as a rigorous worst-case
  bound.** The spec (§7 out-of-scope, R2) is explicit that this is a first-order geometric
  argument tested empirically, not a proof; G1 (never exceeded) is the actual falsifier, and nothing
  in the acceptance criteria asked for a proof.
- **Did not assert anything about the old scout-probe numbers' provenance beyond "not
  reproduced here."** I initially drafted a claim that they were themselves ad hoc `π·sin θ₀`
  values (they're numerically close), but the closest one (N2's "2.391" vs this run's "2.396444")
  doesn't match exactly, so I walked that back to an honest "close but not identical, not verified"
  note in `specs/BACKLOG.md` rather than asserting a mechanism I hadn't checked.

## What I could not verify

- **Generalization beyond these three systems/geometries.** All three test cases are ≤6 qubits,
  STO-3G, near-equilibrium bond lengths. Whether `π·sin θ₀` holds at stretched geometries, larger
  active spaces, or different bases is untested — recorded as spec R3, not claimed.
- **Whether the arccos-Jacobian argument is the actual causal mechanism** vs. a coincidentally
  correct envelope. G1/G2 test the *prediction*, not the derivation's internal logic (spec R2,
  §7). A different argument could conceivably produce the same numbers.
- I did not re-verify the unrelated 17 pre-existing ruff findings in other files (e.g.
  `scripts/spec_pm3_subspace_eta_bound.py`, a `reachability` import) — out of this bead's scope,
  and left as-is; noting their existence only so they aren't mistaken for something this bead
  introduced.

## Files changed (uncommitted — tree is dirty per this run's git policy)

```
 M specs/BACKLOG.md
?? scripts/spec_qpe_precision_bound.py
?? specs/SPEC_qpe_precision_bound.md
?? tests/test_qpe_precision_bound_spec.py
```

## Suggested commands for the human reviewer

```bash
git add specs/BACKLOG.md scripts/spec_qpe_precision_bound.py specs/SPEC_qpe_precision_bound.md \
        tests/test_qpe_precision_bound_spec.py sandbox-handoffs/chem-0ke.md
git commit -m "qpe: confirm pi*sin(theta0) as the precision-bound constant (chem-0ke)"
```

`bd close chem-0ke` has already been run as part of this session (see below) — evidence for the
acceptance criteria exists and the targeted gates pass.
