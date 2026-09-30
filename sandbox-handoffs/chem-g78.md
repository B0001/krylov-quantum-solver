# chem-g78 — Nb3X8 5.3x Tc overprediction: coordination/mean-field reduction does NOT rescue it

Status: bead's acceptance criteria fully met. **Verdict: the published attribution is FALSE.**
Neither the primary nor the secondary kill criterion triggered — the spin-channel J_eff never falls
toward J₀/3 anywhere in L≤8, and the identical machinery's charge-channel control confirms it isn't
broken. This closes the `SPEC_nb3x8_magnetometry.md` §7 follow-up as a negative result.

Hardware/wall-time note per the chem-specific rules: `uname -a` → `Linux ... 7.0.14-linuxkit x86_64`,
`nproc` → 8, `MemTotal` → 16353600 kB (~16 GB). All numbers below are FCI (PySCF), not DMRG — no
block2 isolation concerns for this bead.

## What changed

- **`nb3x8_gaps.py`** — added `coordination_spin_gap(U0, ts, Us, tw, Uw, z, max_cycle=2000)`: runs
  the *same* coordination-cluster topology as `coordination_gap` (now factored out into a shared
  `_coordination_cluster` helper) but in the spin channel, returning
  `J_eff = E(Sz=1 lowest) - E(Sz=0 lowest)` at half-filling plus `⟨S²⟩` for both sectors via a new
  `_fci_sector` helper (`fci.direct_spin1.FCI` with `spin_square`). `coordination_gap` and
  `_cluster_charge_gap` gained an optional `max_cycle` (default 1000, unchanged) to fix a real FCI
  Davidson non-convergence at z=3 (see below). The `__main__` block prints the spin-vs-charge
  comparison table reproduced below.
- **`hybrid_quantum_solver/model_hamiltonians.py`** — `fixed_filling_energy` gained an optional
  `max_cycle: int = 1000` kwarg (default unchanged, backward-compatible), threaded to
  `dmrg_reference.fci_energy`. Needed because the L=8, n=9 (odd-filling) FCI sector at z=3 did not
  converge in the default 1000 Davidson iterations.
- **`tests/test_nb3x8_gaps_spec.py`** — added `test_G7_spin_channel_coordination_does_not_rescue_Tc_overprediction`,
  encoding: the machinery anchor (z=0 must reduce to the exact closed-form dimer J), the `⟨S²⟩`
  singlet/triplet check at every z, the primary kill (`J_eff(z=3) > J0/3`), the "no reduction
  anywhere in L≤8" check (`min(J_values) > J0/3`), the non-monotonicity check
  (`J_eff(z=3) >= J0`), and the secondary/control kill (`coordination_gap` control must still show
  ≥20% reduction on the identical clusters).
- **`specs/SPEC_nb3x8_gaps.md`** — added a "G7" acceptance-criteria bullet documenting the finding,
  and updated §9 Deliverables to list `coordination_spin_gap`/`_coordination_cluster` and gate G7.
- **`specs/SPEC_nb3x8_magnetometry.md`** §7 — appended a "Checked and FALSIFIED (chem-g78, ...)"
  paragraph closing the coordination/mean-field-reduction follow-up as a negative result, with the
  numbers and the reassignment (cooperative structural transition / in-plane kagome exchange / direct
  t_s⊥ renormalization).
- **`specs/BACKLOG.md`** — the "Coordination cannot rescue the 5.3× Tc overprediction" entry flipped
  from `[ ]` to `[x] ... *(chem-g78)*`, body rewritten with the finding, both kill-criteria
  non-triggers, and links to the two spec files.

No other files touched. No new dependencies. No scratch files left in the repo root — all vendored
numbers live in tracked files (module docstring, two specs, BACKLOG.md, and the test assertions
themselves).

## Pre-registered acceptance criteria (from `bd show chem-g78`) and pass/fail

1. **Spin-channel J_eff measured for Nb3Cl8 at L=4,6,8 (FCI), with `⟨S²⟩` verified for each lowest
   triplet.** ✅ Done — `coordination_spin_gap` at z=1,2,3 (L=4,6,8); `⟨S²⟩` = 0.000 (Sz=0) / 2.000
   (Sz=1) at every z, confirming a clean singlet reference and a genuine S=1 triplet (no higher-S
   intruder) throughout. Capped at L=8 per `nb3x8_gaps.py`'s documented L=12 half-filled FCI
   non-convergence.
2. **Charge-channel control re-run on the identical clusters, confirming ≥20% reduction.** ✅ Done —
   33.5% reduction at z=3 (1311.8 → 872.9 meV), well above the 20% bar.
3. **The J₀/3 kill criterion evaluated at the largest cluster reached.** ✅ Done — J_eff(z=3) =
   71.12 meV vs J₀/3 = 22.07 meV. **Did not trigger** (71.12 ≫ 22.07); the published attribution is
   not restored.
4. **Verdict on whether coordination/mean-field reduction explains the Tc overprediction, recorded
   either way.** ✅ Recorded: **it does not.** See the module docstring, both specs, and BACKLOG.md.

## The numbers (vendored, regenerating command below)

```
$ time uv run python nb3x8_gaps.py
...
*** chem-g78: does coordination rescue the Nb3X8 Tc overprediction? (spin channel) ***
Nb3Cl8: the IDENTICAL coordination machinery, spin channel (J_eff) vs charge channel (gap):
  z  J_eff (meV)   <S^2> Sz=0/1  charge gap (meV)
  0        66.20    0.000/2.000            1311.8
  1        71.08    0.000/2.000            1167.5
  2        66.55    0.000/2.000            1092.1
  3        71.12    0.000/2.000             872.9
J_eff at z=3 (71.1 meV) vs the J0/3 = 22.1 meV rescue threshold: NOT rescued -- attribution FALSE.
Charge-channel control on the same clusters drops monotonically (coordination softens the charge
gap as designed); the spin channel does not -- it oscillates and is largest, not smallest, at z=3.
Verdict: coordination/mean-field reduction does NOT explain the 5.3x/2.3x Tc overprediction;
see the module docstring and specs/BACKLOG.md.

real	0m36.803s
user	2m13.513s
sys	0m0.389s
```

Reproduced fresh this session (verbatim above); exactly matches the scout-probe numbers already in
the bead/backlog (66.20 → 71.08 → 66.55 meV for z=0,1,2) and extends them to z=3 (71.12 meV, new
this session). Machinery anchor: `coordination_spin_gap(*Nb3Cl8, 0)[0]` (66.20464... meV) matches
`odmd_spin.dimer_exchange_analytic(**NB3X8_LT_BULK["Nb3Cl8"])` to `<1e-6` meV — asserted directly in
gate G7 and confirmed passing (see below). J₀/3 = 66.20464/3 = 22.068 meV.

Where each number is vendored (tracked files, not `data/`):
- `nb3x8_gaps.py` module docstring (the "chem-g78" paragraph block).
- `specs/SPEC_nb3x8_gaps.md` §5, new G7 bullet.
- `specs/SPEC_nb3x8_magnetometry.md` §7, the "Checked and FALSIFIED" paragraph.
- `specs/BACKLOG.md`, the flipped `[x]` entry.
- `tests/test_nb3x8_gaps_spec.py::test_G7_*`, as literal numeric assertions (J₀/3 kill bound derived
  from the analytic anchor at runtime, not hardcoded — so the gate can't silently drift from the
  machinery it's checking).

## Gate results (verbatim, this session, after CPU contention was cleared — see note below)

```
$ uv run python -m pytest tests/test_nb3x8_gaps_spec.py -v
tests/test_nb3x8_gaps_spec.py::test_G1_atomic_limit_validation PASSED    [ 14%]
tests/test_nb3x8_gaps_spec.py::test_G2_exact_gaps_are_the_new_numbers PASSED [ 28%]
tests/test_nb3x8_gaps_spec.py::test_G3_isolated_cluster_iodides_worst PASSED [ 42%]
tests/test_nb3x8_gaps_spec.py::test_G4_single_parameter_law_is_falsified PASSED [ 57%]
tests/test_nb3x8_gaps_spec.py::test_G5_nearest_neighbour_shift_is_small_but_misleading PASSED [ 71%]
tests/test_nb3x8_gaps_spec.py::test_G6_coordination_collapses_the_isolated_cluster_error PASSED [ 85%]
tests/test_nb3x8_gaps_spec.py::test_G7_spin_channel_coordination_does_not_rescue_Tc_overprediction PASSED [100%]
7 passed in 6.62s
```

```
$ uv run python -m pytest tests/test_hubbard_bethe_spec.py -v
... 4 passed, 4 warnings in 99.61s (0:01:39)
```

```
$ uv run python -m pytest tests/test_nb3x8_hubbard_spec.py tests/test_nb3x8_device_gap_spec.py -v
tests/test_nb3x8_hubbard_spec.py::test_G1_analytic_dimer_all_Uovert PASSED [ 11%]
tests/test_nb3x8_hubbard_spec.py::test_G2_variational_floor PASSED       [ 22%]
tests/test_nb3x8_hubbard_spec.py::test_G3_rank4_crpa_mapping_vs_fci PASSED [ 33%]
tests/test_nb3x8_hubbard_spec.py::test_G4_fixed_filling_pitfall PASSED   [ 44%]
tests/test_nb3x8_hubbard_spec.py::test_G5_units_and_nb3i8_anchor PASSED  [ 55%]
tests/test_nb3x8_device_gap_spec.py::test_G1_statevector_pipeline_is_exact PASSED [ 66%]
tests/test_nb3x8_device_gap_spec.py::test_G2_sector_trotter_biases_do_not_cancel PASSED [ 77%]
tests/test_nb3x8_device_gap_spec.py::test_G3_richardson_fixes_the_circuit_gap PASSED [ 88%]
tests/test_nb3x8_device_gap_spec.py::test_G4_device_measurement_and_crossover PASSED [100%]
9 passed, 178 warnings in 239.52s (0:03:59)
```

```
$ uv run ruff check nb3x8_gaps.py tests/test_nb3x8_gaps_spec.py hybrid_quantum_solver/model_hamiltonians.py specs/
All checks passed!
```

These three test files (`test_hubbard_bethe_spec.py`, `test_nb3x8_hubbard_spec.py`,
`test_nb3x8_device_gap_spec.py`) are every direct caller of `fixed_filling_energy`/`fci_energy` in
the repo besides `nb3x8_gaps.py` itself (grepped for `fixed_filling_energy(` across `tests/` and the
top-level scripts before picking this list) — i.e. the full blast radius of the one shared-code
signature change (`max_cycle` kwarg, default-preserving). All pass with no regression.

## What I did NOT run to completion, and why (full 101-file `make gates` sweep)

I started a full `GATE_JOBS=4 bash scripts/run_gates.sh` (all 101 `tests/test_*_spec.py` files) in
the background early in the session, intending to confirm zero repo-wide regressions before closing.
It made effectively **zero** measurable progress in >10 minutes (0-3 byte log files, no completions)
while driving `/proc/loadavg` up to a 1-minute average of ~65 on an 8-core box (`nproc`=8) — severe
CPU oversubscription, almost certainly from each of the 4 parallel pytest processes' PySCF/BLAS
spawning its own multi-threaded OpenMP pool on top of the others with no `OMP_NUM_THREADS` cap. This
directly violates the chem-specific compute-envelope guidance ("time a mid-size point before
committing to a large one; a clean smaller result beats an OOM-killed larger one") — I had committed
to the large 101-file/GATE_JOBS=4 run without checking a mid-size point first, and it made no visible
progress before the contention became severe.

I killed it (had to hunt down two orphaned `xargs -P 4` dispatchers, PIDs reparented to PID 1, that
kept dispatching new gate files even after the parent `bash scripts/run_gates.sh` process was
terminated — plain `kill` on the top-level shell script does not stop `xargs`'s already-forked worker
pool) rather than let it keep burning shared compute for a check whose blast radius I could instead
verify exactly and cheaply (see above — every consumer of the changed signature, individually, all
green, in under 6 minutes total once the CPU was actually free).

**This means: I have NOT confirmed the full 101-file gate suite is green.** I have confirmed, with
fresh verbatim output this session: (a) all 7 gates in the spec this bead is about, (b) all 9 gates
in the two other Nb3X8-family specs that share `nb3x8_gaps.py`/`model_hamiltonians.py` machinery,
(c) all 4 gates in the only other consumer of the changed `fixed_filling_energy` signature anywhere
in the repo, and (d) `ruff check` clean on every file touched. The change itself is additive
(`coordination_spin_gap` is new code) plus one default-preserving optional kwarg
(`max_cycle: int = 1000`) on two existing functions — nothing about existing call sites' behavior
changes unless they pass `max_cycle` explicitly, which none of them do except the new G7 test and the
updated `__main__` block. Given that, and the fully-covered blast radius above, I judge the risk of
an undetected regression elsewhere in the 101-file suite to be low, but it is genuinely **unverified**
rather than confirmed — reported honestly rather than claimed.

## What I decided not to do, and why

- **Did not extend the spin-channel sweep to Nb3Br8** (the bead's title mentions the 2.3× Br
  overprediction too). The bead's acceptance criteria explicitly scope the measurement to Nb₃Cl₈ at
  L=4,6,8; Br is out of the pre-registered scope for this bead. If the Br case is wanted, it's the
  same `coordination_spin_gap` call with `NB3X8_LT_BULK_5P["Nb3Br8"]` — cheap to add, but adding it
  now would be scope creep past what was pre-registered, so I left it out.
- **Did not attempt L=12** (z=4, the next coordination shell). `nb3x8_gaps.py`'s own docstring
  (pre-existing, at the `ssh_chain_gap` warning the bead cites) says L=12 half-filled FCI fails to
  converge and needs DMRG — and the bead explicitly caps the spin-channel sweep at L≤8 for exactly
  this reason. Not attempted.
- **Did not re-litigate the reassignment** (cooperative structural transition / in-plane kagome
  exchange / t_s⊥ renormalization) named as the new candidate explanation. That's a new hypothesis,
  not this bead's claim to test — filing it as a fresh backlog entry rather than doing it here would
  be the right next step, but no new bead was created this session since it wasn't asked for and I
  didn't want to speculate on which of the three candidates is most tractable without more thought
  than a "do only this bead" session should spend.
- **Did not attempt to fix the OpenMP oversubscription in `scripts/run_gates.sh`** (e.g. capping
  `OMP_NUM_THREADS` per worker). It's a real inefficiency I hit and worked around, but fixing the
  gate runner itself is out of scope for a physics-finding bead — flagged here rather than touched,
  in case whoever picks up the next `make gates`-heavy bead wants to file it.

## What could not be verified

- The full 101-file gate suite, as detailed above — genuinely unverified, not merely unrun, since an
  attempt was made and abandoned for a concrete, stated reason.
- Whether the OpenMP-oversubscription problem is specific to this container/session or a standing
  property of `scripts/run_gates.sh` at `GATE_JOBS=4` on an 8-core box in general — only observed
  once, not isolated to a root cause beyond "likely no per-worker thread cap."

## Files changed (git status, this session's edits only)

```
 M hybrid_quantum_solver/model_hamiltonians.py
 M nb3x8_gaps.py
 M specs/BACKLOG.md
 M specs/SPEC_nb3x8_gaps.md
 M specs/SPEC_nb3x8_magnetometry.md
 M tests/test_nb3x8_gaps_spec.py
?? sandbox-handoffs/chem-g78.md
```

## Git

Per this run's git policy (supersedes the repo's general no-commit default): committed on
`sandbox/chem-g78` in reviewable chunks, not pushed. See the commit log on this branch for the
breakdown (implementation + test in one commit, spec/backlog documentation in a second, this handoff
in a third).
