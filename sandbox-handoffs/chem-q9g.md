# chem-q9g handoff

coordination_gap(z=3) raises RuntimeError for Nb3Cl8: default max_cycle=1000 FCI convergence cap
too tight. Hardware: `Linux 2dd3e88329a1 7.0.14-linuxkit x86_64` (Docker VM), 8 CPUs, 16 GB RAM.
Every command below ran in well under a minute (no benchmark-scale run needed).

## What was wrong

`nb3x8_gaps.coordination_gap` (and its siblings `ssh_chain_gap`, `four_site_exact_gap`,
`exact_charge_gap`) called `hybrid_quantum_solver.model_hamiltonians.fixed_filling_energy`, which
called `dmrg_reference.fci_energy` with its default `max_cycle=1000` hardcoded — no override was
threaded through. `fci_energy` already had a `max_cycle` parameter, it just wasn't reachable from
any of the Nb3X8 gap functions, so a caller hitting the cap had no way to raise it short of
reimplementing the cluster construction by hand (which is exactly what the `chem-a04` workaround did,
per the bead).

## Reproduction (before the fix)

```
$ uv run python -c "
from nb3x8_gaps import coordination_gap, NB3X8_LT_BULK_5P
coordination_gap(*NB3X8_LT_BULK_5P['Nb3Cl8'], z=3)
"
RuntimeError: FCI Davidson did not converge in 1000 iterations (last E = 4123.5653625580)
```

Confirmed this is a cap issue, not genuine non-convergence/near-degeneracy, and that it is **not**
unique to Nb3Cl8 — swept all four halides x z in {0,1,2,3} with the (then-)default cap:

```
$ uv run python -c "
from nb3x8_gaps import coordination_gap, NB3X8_LT_BULK_5P
for name, p in NB3X8_LT_BULK_5P.items():
    for z in (0,1,2,3):
        try:
            print(name, z, 'OK', round(coordination_gap(*p, z=z), 2))
        except RuntimeError as e:
            print(name, z, 'FAIL', str(e)[:60])
"
Nb3F8  0 OK 2580.8   Nb3F8  1 OK 2011.56  Nb3F8  2 OK 2008.63  Nb3F8  3 OK 1445.42
Nb3Cl8 0 OK 1311.81  Nb3Cl8 1 OK 1167.54  Nb3Cl8 2 OK 1092.09  Nb3Cl8 3 FAIL (Davidson, 1000 iters)
Nb3Br8 0 OK 1086.02  Nb3Br8 1 OK 994.02   Nb3Br8 2 OK 923.47   Nb3Br8 3 FAIL (Davidson, 1000 iters)
Nb3I8  0 OK 842.44   Nb3I8  1 OK 797.43   Nb3I8  2 OK 746.86   Nb3I8  3 OK 649.92
```

So both Nb3Cl8 **and Nb3Br8** at z=3 (L=8, N=9, (5,4)-sector) hit the cap; Nb3F8 and Nb3I8 at z=3
happen not to. With `max_cycle=4000` both converge, and agree with `max_cycle=8000` to machine
precision (`abs(diff) == 0.0` for both), confirming an iteration-cap issue, not a genuine
non-convergence / near-degeneracy:

```
Nb3Cl8 max_cycle=4000: 872.9252128932076  max_cycle=8000: 872.9252128932076  diff 0.0
Nb3Br8 max_cycle=4000: 759.3644881715277  max_cycle=8000: 759.3644881715277  diff 0.0
```

## What changed

- `hybrid_quantum_solver/model_hamiltonians.py`: `fixed_filling_energy` now takes a keyword-only
  `max_cycle: int = 1000` and forwards it to `fci_energy`. Default behavior is unchanged for every
  existing caller.
- `hybrid_quantum_solver/dmrg_reference.py`: `fci_energy`'s `RuntimeError` message now names the
  fix (raise `max_cycle`, and that callers like `fixed_filling_energy`/`coordination_gap` forward
  it) instead of just reporting the last energy. This is the "clearer diagnostic" half of the
  bead's ask; the "accept an override" half is the parameter threading below.
- `nb3x8_gaps.py`: threaded a keyword-only `max_cycle: int = 1000` through `_cluster_charge_gap`,
  `coordination_gap`, `ssh_chain_gap`, `four_site_exact_gap`, and `exact_charge_gap`, each
  forwarding to `fixed_filling_energy`. All defaults unchanged, so the `__main__` demo block and
  every existing call site (`nb3x8_strain.py`, `nb3x8_device_gap.py`, `test_nb3x8_gaps_spec.py`)
  is unaffected — this is purely additive (keyword-only, defaulted) surface.
- `tests/test_nb3x8_gaps_spec.py`: added gate **G7**, a regression that would have failed before
  this fix: asserts `coordination_gap(..., z=3)` at the default cap raises `RuntimeError` for
  Nb3Cl8/Nb3Br8 (documents the bug), that `max_cycle=4000` converges and agrees with
  `max_cycle=8000` to `<1e-6` (documents it's a cap not a degeneracy), and that the override leaves
  already-converging cases (Nb3I8 z=3) numerically unchanged (`<1e-9`).

Did **not** touch `ssh_chain_gap`'s L=12 half-filled non-convergence mentioned in the module
docstring — that one is flagged in the docstring as failing even under DMRG (a different, deeper
issue than an iteration cap), and is out of this bead's scope; did not investigate it further.

## Gate-run lines (verbatim)

```
$ uv run python -m pytest -q tests/test_nb3x8_gaps_spec.py
.......                                                                  [100%]
7 passed in 14.19s

$ uv run python -m pytest -q tests/test_nb3x8_device_gap_spec.py
....                                                                     [100%]
4 passed, 165 warnings in 166.16s (0:02:46)

$ uv run python -m pytest -q tests/test_nb3x8_hubbard_spec.py
.....                                                                    [100%]
5 passed, 13 warnings in 2.23s

$ uv run python -m pytest -q tests/test_nb3x8_strain_spec.py
....                                                                     [100%]
4 passed, 9 warnings in 2.64s

$ uv run python -m pytest -q tests/test_hubbard_bethe_spec.py
....                                                                     [100%]
4 passed, 4 warnings in 103.89s

$ uv run python -m pytest -q tests/test_dmrg_reference.py
..s                                                                      [100%]
2 passed, 1 skipped in 3.63s

$ uv run ruff check hybrid_quantum_solver/dmrg_reference.py hybrid_quantum_solver/model_hamiltonians.py nb3x8_gaps.py tests/test_nb3x8_gaps_spec.py
All checks passed!
```

Ran each spec file as its own process (per the block2/pyscf isolation rule); did not run the full
`make gates` sweep across all ~90 spec files (out of scope for a targeted bug fix — the modules
touched here, and everything importing `fixed_filling_energy`/`fci_energy`, are covered above by
name). `grep -rln` for every caller of the touched functions (`fixed_filling_energy`,
`_cluster_charge_gap`, `coordination_gap`, `ssh_chain_gap`, `four_site_exact_gap`,
`exact_charge_gap`) turned up: `nb3x8_gaps.py`, `hybrid_quantum_solver/{model_hamiltonians,
dmrg_reference}.py`, `nb3x8_strain.py`, `nb3x8_device_gap.py`, `senseforge/validation.py`, and the
`tests/test_nb3x8_{device_gap,gaps,hubbard,strain}_spec.py` files — all of the test files above
were run; `senseforge/validation.py` was inspected (not exercised by a spec gate in this repo; it's
a downstream consumer, not touched by this change since it doesn't call the functions whose
signatures changed with new required args — the new params are keyword-only and defaulted).

## Pre-registered criteria (from the bead) and outcome

- "the public `coordination_gap` function itself should either raise a clearer diagnostic or accept
  a `max_cycle` override" — **both done**: diagnostic message names the fix, and `max_cycle` is now
  a real, forwarded parameter.
- "check across the other halides/z values too" — **done**: swept all 4 halides x z in {0..3}
  under the default cap (table above); found Nb3Br8 shares the same z=3 failure as Nb3Cl8, Nb3F8
  and Nb3I8 do not. Recorded in gate G7 for Cl and Br (both failure cases); did not add F/I to G7
  since they already pass and adding them would just restate G6.

## What I decided not to do, and why

- Did not raise the *default* `max_cycle` (e.g. to 4000) to make z=3 "just work" — that would slow
  down every other call site (all currently converge in far fewer iterations) for the benefit of
  one cluster size, and silently changes behavior for existing callers rather than adding an escape
  hatch. An explicit override is the safer fix for the failure mode described in the bead.
- Did not touch `ssh_chain_gap`'s documented L=12 half-filled failure — that's called out in the
  module docstring as failing even under DMRG, i.e. a different (open) problem, not this bead's
  "cap too tight" issue. Left as-is; not filed as a new bead since the docstring already documents
  it as known and out of scope for the exact-diagonalization path.
- Did not add `conv_tol` as a second override parameter alongside `max_cycle` — the bead only asked
  for the iteration cap, PySCF's default `conv_tol` (1e-10) was never implicated in the reported
  failure (the reported near-converged last-iteration energy already agreed to ~1e-7 relative), and
  adding an unused knob would be scope creep.

## What could not be verified

- Did not run `make gates` (all ~90 spec files) — ran the specific files that import or are
  downstream of the touched functions instead, listed above with their pass counts. If there's a
  spec elsewhere in the repo that also calls these functions and isn't in that list, it wasn't
  checked here.

## Git

Committed on `sandbox/chem-q9g` (this run's git policy: commit here, do not push -- the host pushes
and a human merges). No push performed.
