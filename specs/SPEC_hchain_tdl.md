# SPEC: Hydrogen-chain DMRG with bond-dimension + thermodynamic-limit extrapolation

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates below are agreed.

---

## 1. Goal

Compute the ground-state energy **per atom in the thermodynamic limit**, `e_∞`, of the
minimal-basis hydrogen chain H_n, by:

1. **Bond-dimension extrapolation** `D → ∞` at each finite `n` (DMRG truncation control), and
2. **Finite-size extrapolation** `n → ∞` of the per-atom energy.

This turns the existing one-off `benchmark_dmrg_large.py` (two bond dimensions, no extrapolation)
into a *controlled* strongly-correlated calculation with quantified uncertainty.

## 2. Background and honest framing

The H_n chain is the canonical *ab initio* model of the Mott metal–insulator transition and the
reference system of the Simons Collaboration benchmark (Motta et al., *PRX* 7, 031059, 2017).

**What we can claim:** a controlled, cross-validated `e_∞` for the **minimal-basis model** at a
fixed geometry, with honest error bars from both extrapolations.

**What we cannot claim (stated up front, not discovered later):**
- This is **not** the complete-basis-set (physical) energy — minimal basis only (see §7).
- This **reproduces** known benchmark physics; it does not extend the frontier.
- The quantum (statevector) solver plays **no role** — these `n` are far past its reach. This is a
  classical DMRG study; that is the honest, scalable lever (see REFACTOR_PLAN.md, Phase 4).

## 3. Scientific approach

**Geometry/basis (fixed model):** H_n, uniform spacing `R` (default `R = 1.8` Bohr ≈ 0.953 Å,
near the cohesive minimum), basis `STO-6G` → `n` spatial orbitals, `n` electrons, singlet.

**(a) Bond-dimension extrapolation.** DMRG energy `E(D)` converges to the exact (in-model) energy
as the discarded weight `δ(D) → 0`. The rigorous prescription (White/Chan) is a **linear fit of
`E(D)` vs the discarded weight `δ(D)`**, extrapolated to `δ = 0`:

    E(D) ≈ E(∞) + c · δ(D)

over a schedule of increasing `D` (e.g. `[200, 400, 800, 1600]`). The extrapolated `E(∞)` must lie
**below** the largest-`D` computed energy and (where checkable) **at/above** exact FCI.

> **Technical requirement / risk (R1):** block2 must expose the per-sweep **discarded weight**.
> If it cannot be obtained reliably, the fallback is `E` vs `1/D` extrapolation, which is cruder;
> the spec is met either way but the gate tolerances (below) assume the truncation-error fit.

**(b) Thermodynamic-limit extrapolation.** For an open chain the per-atom energy scales as

    E(n)/n ≈ e_∞ + a/n + O(1/n²)

(a 1/n surface/boundary correction). Fit `E(∞,n)/n` vs `1/n` over `n ∈ {…}` to obtain `e_∞` and its
uncertainty.

## 4. Public interface

Two additions, no changes to the validated solver:

```
hybrid_quantum_solver/dmrg_reference.py
    dmrg_energy_extrapolated(h1, eri, n_elec, e_core, bond_dims=(200,400,800,1600), **kw)
        -> ExtrapResult(energy, stderr, per_D=[(D, E, discarded_weight)], method="dweight"|"invD")

benchmark_hchain_tdl.py          # the study driver (resumable, like benchmark_nbn.py)
    -> data/hchain_tdl.csv        # one row per (n): E(D)-points, E_extrap, stderr, FCI, |err|
    -> prints e_∞ ± stderr and the n→∞ fit
```

`ExtrapResult.energy` is `E(D→∞)`; `stderr` is the fit standard error. The TDL fit lives in the
driver and reports `e_∞`, its standard error, and the fit residual.

## 5. Acceptance criteria (validation gates)

Each is an automated check in `tests/test_hchain_tdl_spec.py` (test-first):

- **G1 — bond-dim extrapolation is *sound* (lands at exact FCI).** For `n = 10` and `n = 12`
  (FCI tractable) with well-converged bond dims `(80, 160, 300)`, `|E_extrap − E_FCI| < 2e-4 Ha`.
  > **Refinement (found during implementation):** the original clause "extrapolation *improves*
  > on the best single-`D` energy" is **unsatisfiable for any FCI-validatable system** — at sizes
  > small enough to have an exact reference, the largest-`D` DMRG energy is itself near-exact, so
  > extrapolation cannot beat it and (with coarse points) can even *overshoot below* FCI. The
  > improvement only materialises in the **under-converged large-`n` regime**, where no FCI exists
  > to check. So G1 validates *soundness* (no spurious bias), not improvement. The `2e-4` tolerance
  > (vs `1e-4`) reflects the bond dims that keep the test fast; it tightens with larger `D`.
- **G2 — variational sanity.** `E_extrap ≤ E(D_max)` and (where FCI known) `E_extrap ≥ E_FCI − 1e-6`.
- **G3 — monotone truncation convergence.** `E(D)` is non-increasing in `D`, and `δ(D)` is
  decreasing, for every `n` tested.
- **G4 — TDL fit stability.** The `n → ∞` estimate `e_∞` changes by `< 1 mHa/atom` when the
  largest-`n` point is dropped from the fit (a leave-one-out self-consistency check).
- **G5 — reproducibility.** Re-running a fixed `(n, D-schedule, seed)` reproduces `E_extrap` to
  `< 1e-8 Ha`; the driver is resumable (skip completed `n`, like `benchmark_nbn.py`).

> Gate **G1** is the definition of "done" for part (a); **G4** for part (b). A literature
> cross-check of `e_∞` against the Simons benchmark value at the chosen `R` is a **stretch goal**,
> not a gate (basis/geometry differences make it indicative, not pass/fail).

## 6. Implementation plan (test-first)

1. Write `tests/test_hchain_tdl_spec.py` encoding G1–G5 (initially failing).
2. Implement `dmrg_energy_extrapolated` (truncation-error fit; `1/D` fallback per R1).
3. Implement `benchmark_hchain_tdl.py` (resumable driver + TDL fit).
4. Iterate until G1–G5 are green. Run block2 tests in their own process (the OpenMP-isolation
   pattern already in `run_in_chem.sh`).

## 7. Out of scope

- **Complete-basis-set limit** (cc-pVxZ extrapolation) — minimal basis only; the result is the
  *model* energy, not the physical one.
- **Periodic boundary conditions** — open chains only (hence the `1/n` surface term).
- **Bond-length/equation-of-state scan** — single fixed `R` (a follow-up spec could add the curve).
- **The quantum/statevector solver** — out of reach at these `n` by construction.
- Any **novelty claim** — this is a reproduction/validation study.

## 8. Caveats and risks

- **R1 (discarded weight access):** see §3(a). Mitigation: `1/D` fallback.
- **Minimal basis** makes `e_∞` a model number; do not present it as the physical H-chain energy.
- **Open-boundary surface effects** can bias the TDL fit at small `n`; mitigated by G4 and by
  starting the fit at `n ≥ 10`.
- **DMRG cost** grows with `n` and `D`; the largest `(n, D)` points are the runtime pole. The driver
  is resumable so partial runs are never lost.

## 9. Deliverables

- `hybrid_quantum_solver/dmrg_reference.py` — `dmrg_energy_extrapolated` (+ `ExtrapResult`).
- `benchmark_hchain_tdl.py` — resumable study driver → `data/hchain_tdl.csv`.
- `tests/test_hchain_tdl_spec.py` — gates G1–G5.
- A short results summary (the `e_∞ ± stderr` table) in the PR description, with the §2/§7 caveats.
