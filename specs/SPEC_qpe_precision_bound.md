# SPEC: QPE's precision constant is not a measurement — it is a per-system closed form

**Status:** CLOSED (chem-0ke) — CONFIRMED. G1–G3 green (`tests/test_qpe_precision_bound_spec.py`,
7 passed); measured numbers in §10.

---

## 1. Goal

`SPEC_qpe_readout_laws` G2 measured the argmax-Fejer-kernel QPE precision bound's constant as
`~2.175` on H2 CAS(2,2) and explicitly disclaimed a derivation ("not derived analytically", "not
universal"). This spec tests a specific closed-form candidate: `run_qpe` decodes
`E = lambda*cos(2*pi*phi)`, so a half-bin dyadic phase-estimation error (`|delta-phi| <= 2^-(t+1)`)
pushed through the arccos-to-cos Jacobian gives a first-order bound
`|E_est - E_0| <= pi*sin(theta_0) * lambda/2^t`, where `theta_0 = arccos(E_0/lambda)` is computable
per-system from the exact spectrum alone, before ever running a QPE sweep. False if the measured
max ratio `err(t)/(lambda/2^t)` ever exceeds `pi*sin(theta_0)` on any tested system, or if the bound
is loose by more than 20% (i.e. `max_ratio < 0.8 * pi*sin(theta_0)`) on **every single** system
tested (a real bound, but the wrong mechanism).

## 2. Background and honest framing

- This is a direct follow-up to `SPEC_qpe_readout_laws` G2/R1: that spec measured a single number
  (2.175) on a single system (H2 CAS(2,2)) and explicitly said "not derived analytically" and "not
  universal". `specs/BACKLOG.md` (Fault-tolerant stack section) proposed the arccos-Jacobian
  derivation below as the closer.
- **What you can claim if the gates pass:** `pi*sin(theta_0)` is a real, per-system-predictable
  upper envelope for the argmax point-estimate error's `err(t)/(lambda/2^t)` ratio, computable from
  the exact ground energy and the block-encoding 1-norm alone — no QPE sweep required to know the
  bound in advance.
- **What you cannot claim:** that the bound is *tight* at every `t` (it is an upper envelope; a
  short sweep can sit well under it by chance — see Caveat R1), or that the mechanism is fully
  characterized (the derivation is a first-order geometric argument, not a rigorous worst-case proof
  over all possible bin alignments).
- **Reference:** the exact ground energy and exact diagonalization of the qubit Hamiltonian
  (`np.linalg.eigh`), same reference `qpe_walk_readout.py` and `SPEC_qpe_readout_laws` already use.
  `theta_0` and the predicted bound are derived from this exact spectrum directly — never from a
  QPE run — so the check does not let the producer grade its own homework.

## 3. Approach — pre-registration first

**Before running any sweep**, compute per system: `lambda = sum |c_l|` (Pauli 1-norm),
`E_0` = exact ground electronic energy, `theta_0 = arccos(E_0/lambda)`, and the predicted bound
`pi*sin(theta_0)`. These are written down here, in this file, as the committed prediction:

| System        | qubits | lambda   | E_0 (elec.) | theta_0 (rad) | **predicted bound = pi·sin(theta_0)** |
|---------------|-------:|---------:|------------:|--------------:|---------------------------------------:|
| H2 CAS(2,2)   |      4 | 2.699278 |   -1.852388 |       2.327122 |                               **2.285077** |
| LiH CAS(2,2)  |      4 | 1.595641 |   -1.058117 |       2.295789 |                               **2.351497** |
| N2 CAS(3,4)   |      6 | 8.422876 |   -5.446405 |       2.273942 |                               **2.396444** |

(Geometries: H2 `H 0 0 0; H 0 0 0.74`; LiH `Li 0 0 0; H 0 0 1.6`; N2 `N 0 0 0; N 0 0 1.10` with
`symmetry="D2h"` to pin the HOMO/LUMO degenerate-pair gauge (see `SPEC_lambda_ladder_honest_caveat`
R2/chem-1yr) — all STO-3G, all built the same way `qpe_walk_readout.py`'s own `__main__` and
`taper_qubits.py`'s reference harness do. CAS(n,m) notation here follows this repo's existing
convention from those files: `mcscf.CASCI(mf, n, m)`, i.e. `n` orbitals / `m` electrons.)

Reproduce this table with `uv run python scripts/spec_qpe_precision_bound.py` — pinned as a
tracked script so the prediction is regenerable and auditable, not just prose.

Then: sweep `t = 4..14` on all three systems with the EXACT ground state as trial (as in
`SPEC_qpe_readout_laws` G1/G2), record `max_t err(t)/(lambda/2^t)` per system, and compare against
the table above.

## 4. Public interface

No new library code. Reuses `qpe_walk_readout.run_qpe`, `qubitization_blueprint.build_qubit_hamiltonian`,
`qubitization_blueprint.pauli_decompose` unchanged. One tracked script,
`scripts/spec_qpe_precision_bound.py`, prints the pre-registration table above (no
QPE sweep inside it — only the exact-spectrum-derived prediction).

## 5. Acceptance criteria (validation gates)

- **G1 — never exceeded.** On every system (H2 CAS(2,2), LiH CAS(2,2), N2 CAS(3,4)), for every
  `t = 4..14`: `err(t) <= pi*sin(theta_0) * lambda/2^t` (equivalently
  `max_t ratio(t) <= pi*sin(theta_0)`). A single violation on any system, at any `t`, kills the
  derivation (wrong mechanism).
- **G2 — tight on at least one system.** `max_t ratio(t) >= 0.8 * pi*sin(theta_0)` on at least one
  of the three systems (i.e. not loose by more than 20% everywhere). If the bound is a real bound
  but never within 20% on *any* system, the mechanism does not explain the observed constant even
  though it happens to be a valid upper bound.
- **G3 — reproducibility of the pre-registration.** The table's `theta_0` / bound values, recomputed
  from the exact spectrum at test time, match the committed table to `1e-6` (protects against the
  table being hand-edited out of sync with the code that generates it).

> Definition of done: **G1 and G2 both hold.** G3 is a bookkeeping gate that keeps the committed
> prediction honest, not a physics claim.

## 6. Implementation plan (test-first)

1. Commit `scripts/spec_qpe_precision_bound.py`, run it, paste its output into §3's
   table (done above — see the actual run log in the PR/handoff).
2. Write `tests/test_qpe_precision_bound_spec.py` encoding G1-G3.
3. No changes to `qpe_walk_readout.py` or `qubitization_blueprint.py`.
4. `uv run pytest -q tests/test_qpe_precision_bound_spec.py`, ruff clean.

## 7. Out of scope

- A rigorous worst-case proof of the arccos-Jacobian argument over all bin alignments (this spec
  tests the *prediction*, not the derivation's mathematical rigor).
- Systems beyond H2/LiH/N2 CAS(3,4) at <=6 qubits (`pauli_decompose` is exponential in qubit count).
- Explaining `SPEC_qpe_readout_laws` G4's separate, still-open `f(t)` non-monotonicity finding (a
  different quantity — success-probability ratio, not precision-error ratio).

## 8. Caveats and risks

- **R1 — envelope, not a tight-at-every-t prediction.** `pi*sin(theta_0)` bounds the worst-case
  ratio over the dyadic phase grid; at any single `t` the realized ratio can sit well under it if
  that `t`'s phase bin happens to align favorably. This is why the sweep must cover enough `t`
  values (4..14, eleven points) to have a real chance of hitting near-worst-case alignment before
  drawing a conclusion — a 2- or 3-point sweep could look like a miss by chance alone.
- **R2 — first-order geometric argument, not a formal proof.** The derivation combines a half-bin
  dyadic error with the arccos-to-cos Jacobian `d(cos theta)/d theta = -sin(theta)` at `theta_0`;
  it is a plausibility argument for *why* the constant is `O(pi*sin(theta_0))`, not a rigorous bound
  on the Fejer-kernel argmax estimator's worst case. G1 is the actual falsifier.
- **R3 — three systems, all <=6 qubits, all STO-3G, all near equilibrium geometry.** Generalization
  to larger active spaces, different bases, or stretched/dissociating geometries (where `theta_0`
  moves and near-degenerate excited states can shift the argmax) is not tested here.

## 9. Deliverables

- `scripts/spec_qpe_precision_bound.py` — pre-registration script (no library changes).
- `tests/test_qpe_precision_bound_spec.py` — gates G1-G3.
- This file's §3 table — the committed, before-the-fact prediction.

## 10. Results

Run `uv run python tests/test_qpe_precision_bound_spec.py` to regenerate the measured side. The
sandbox's numbers reproduce to all 6 decimals on macOS 27 / Apple M3. The pre-registration script
also reproduces the §3 table exactly.

| System       | predicted bound | measured max ratio (t = 4..14) | at t | max/bound |
|--------------|----------------:|-------------------------------:|-----:|----------:|
| H2 CAS(2,2)  |        2.285077 |                       2.175167 |   11 |     95.2% |
| LiH CAS(2,2) |        2.351497 |                       2.284570 |   14 |     97.2% |
| N2 CAS(3,4)  |        2.396444 |                       2.303179 |   14 |     96.1% |

G1 is never exceeded, and G2 is tight to within 5% on all three systems. **CONFIRMED.**

**Landing finding: the first-order envelope does not port as a strict bound.** chem-177's
identity-free walk operator (`run_qpe(recenter=True)`, `SPEC_ft_identity_shift`, landed together
with this spec) moves θ₀ toward π/2. On recentered H₂ at t = 4 the measured ratio is 2.6966, above
`π·sin θ₀'` = 2.6212. That value sits under the exact nearest-bin bound
`2^t·max±|cos(θ₀ ± π/2^t) − cos θ₀|` = 2.7739, which tends to `π·sin θ₀` as t grows. This is the
finite-bin term that R2's "first-order argument, not a proof" leaves out. It does not touch G1–G3,
which were pre-registered on the identity-included walk; there the raw sweep stays under
`π·sin θ₀` through t = 20 on all four of SPEC_ft_identity_shift's systems. "Portable" holds as a
large-t envelope, not as a bound at every t. Regenerate with
`uv run python tests/test_ft_identity_shift_spec.py`. Not gated; see the PR follow-ups.
