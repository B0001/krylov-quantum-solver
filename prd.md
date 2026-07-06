# Product Requirements Document: A Discovery Engine for Correlated-Electron Chemistry

## 1. Purpose & Objective

* **Core Function:** This repo is a *discovery engine*, not just a solver. A validated hybrid
  quantum–classical toolchain (PySCF active-space integrals → vetted Jordan–Wigner qubit
  Hamiltonian → real-time quantum Krylov / ODMD spectroscopy, cross-checked against FCI/CASCI,
  DMRG, analytic limits, and experiment) exists to **produce scientific findings that are both
  valid and profitable**:
  * **Valid** = every number carries a reference (FCI, DMRG, experiment, or analytic limit) and a
    test that could kill it. This is the repo's founding culture ("falsifiable honesty," see
    `REFACTOR_PLAN.md`) and it is non-negotiable.
  * **Profitable** = the finding *moves the field*: a number not in the literature, a boundary
    where a method or law breaks, a mechanism, or a step toward experiment — something whose being
    wrong would embarrass us, and whose being right is worth citing.
* **Target Audience:** Quantum-chemistry / quantum-algorithms researchers and correlated-materials
  theorists evaluating near-term and fault-tolerant methods on strongly-correlated molecular and
  cluster systems.
* **The shift this PRD encodes:** the validated core path and the near-term/FT method rungs are
  largely built and green. The remaining value is not more machinery — it is **pointing the
  validated machinery at questions whose answers are new.** Prefer discovery over plumbing.

## 2. What Counts as a Discovery (the two-axis bar)

Every candidate is scored on both axes *before* it becomes a spec. Record the score in the spec.

* **Validity gate (mandatory, binary):** Is there a ground truth (FCI / DMRG / experiment /
  analytic limit) and a cheap test that would falsify the claim? No reference → not a spec. This
  is unchanged from the existing SDD loop and every existing gate satisfies it.
* **Profitability gate (mandatory, graded):** rank the candidate on:
  1. **Novelty** — is the headline number/finding absent from the literature? (Best: a value a
     published paper *could have* reported but didn't, e.g. the Nb₃X₈ interlayer exchange
     constants, exciton-binding collapse, and χ(T) tables that `arXiv:2501.10320` never
     tabulated.) A *reproduction* is valid and useful — but it must be **labeled** as such and
     justified by what it unlocks (a validated primitive, a cross-check), never dressed as new.
  2. **Transfer** — does the finding generalize beyond the toy that produced it? A *mechanism* or
     a *law* (the visibility shot-cost law, Trotter-bias-in-the-gap, damping-immune eigenphases)
     outranks a single converged energy.
  3. **Consequence** — does being wrong cost something? A finding that survives a test that *could
     have* killed it is worth more than one that never risked falsification.
  4. **Reach toward experiment or strong correlation** — see §4. Numbers that approach a measured
     quantity, or that live where classical exact methods are genuinely hard, are the profitable
     frontier; a "quantum advantage" claim where an exact classical method is trivial is
     explicitly *anti-profitable* (§5).
* **Honest self-scoring is part of the deliverable.** If a candidate is a clean reproduction with
  no advantage at the reachable scale, say so in the spec's caveats — that honesty is what makes
  the *next* novel claim credible.

## 3. Release Definition (Falsifiable Honesty → Falsifiable Discovery)

* **Culture:** Rebuilt after a scientific audit found the original physics core broken. Every
  claimed result has a reference and a test that could kill it. A gate that proves unsatisfiable
  is a **finding**: revise the spec and record why (see `SPEC_hchain_tdl.md` and the `[-]` killed
  entries in `BACKLOG.md`).
* **Core Loop (`specs/README.md`):** hypothesis in `specs/BACKLOG.md` (a *claim + a cheap
  disproving check*) → `SPEC_<slug>.md` (goal, interface, acceptance gates, out-of-scope, honest
  caveats, **and the two-axis score from §2**) → failing `tests/test_<slug>_spec.py` gate →
  minimum code reusing validated primitives → `make gates` until green.
* **Success Criteria:** `make gates` and `make test` are green — the qubit Hamiltonian reproduces
  FCI, the Krylov solver converges to FCI while respecting its variational floor, shot noise stays
  bounded and improves with shots, and every spec gate passes. **And**: each new spec that claims
  novelty names the literature it would embarrass and the reference that bounds it.

## 4. The Discovery Frontier (where the profitable numbers are)

These are the directions this PRD prioritizes. They reuse validated primitives; the value is the
*question*, not new plumbing.

* **Real correlated materials — the Nb₃X₈ thread is the flagship.** The bilayer downfolds to an
  exactly-diagonalizable generalized Hubbard dimer, and the validated ODMD/ED/DMRG stack has
  already produced genuinely new, gated numbers for the family (charge gaps, interlayer exchange
  J, exciton binding, χ(T), photoemission spectra) — each with the honest caveat that it is an
  *isolated cluster*, an upper bound on the broadened solid. **Profitable next steps:** push these
  from isolated-cluster toward the real solid (coordination/TDL corrections that connect to the
  measured ~600–650 meV gap and the 90 K transition), and toward the *strongly-correlated* sectors
  (e.g. the iodide, already shown to be beyond the Heisenberg regime) where the method earns its keep.
* **Toward experiment, not just internal consistency.** The most valuable validity anchor is a
  measured number. Candidates like "Be₂ toward experiment" (CBS+CV well depth toward 929.7 cm⁻¹)
  and any Nb₃X₈ observable comparable to Sheckelton 2017 / Haraguchi 2017 magnetometry convert an
  internal cross-check into a falsifiable prediction against the lab.
* **Strong correlation where FCI is intractable by physics, not just by count.** The NbN CAS(14,14)
  finding is a warning: "FCI-intractable by determinant count" ≠ "strongly correlated." A
  profitable multireference benchmark must be hard because of *physics* (near-degeneracy, genuine
  multireference character named by sector), with DMRG at real bond dimension as the reference.
* **Methodological laws and boundaries that transfer.** The visibility shot-cost law, Trotter-bias
  survival in energy *gaps*, and depolarizing-immune eigenphases are examples: a calibrated,
  predictive rule extracted from a validated observable, usable as an experiment-planning tool.
  These outrank one-off energies.

## 5. Explicitly Out-of-Scope / Anti-Goals (crucial for AI guidance)

* **NO** numeric result without a ground-truth reference and a falsifying test.
* **NO** novelty claim that is actually a reproduction — reproductions are welcome but must be
  labeled `(reproduction)` with the paper cited and the reference that could embarrass them.
* **NO** "quantum advantage" (or any advantage) claim where a classical exact/reference method is
  trivial at the reachable scale (e.g. the 4-qubit Nb₃X₈ dimer trace, H₂-scale QITE). State
  "no advantage at this scale" plainly when true.
* **NO** treating the CIF materials path as a periodic-solid calculation — it builds a *finite
  molecular cluster* with no periodic boundary conditions. Cluster gaps/observables are upper
  bounds on the broadened solid; transition-metal systems are a research target, not a validated
  materials result, until a TDL/coordination step connects them to the bulk or to experiment.
* **NO** importing, "fixing", or building on the quarantined broken core
  (`orchestrate_hybrid_pipeline.py`, `quantum_sampler.py`) — regression fixtures only.

## 6. Technical Architecture Constraints (do not regress)

* **The validated live path (`hybrid_quantum_solver/`):** `chemistry_gateway` extracts active-space
  integrals → `molecular_hamiltonian` builds the JW `SparsePauliOp` + HF reference →
  `quantum_krylov_solver` estimates the ground state and is bounded below by the true ground state
  → `pipeline.run_geometry` returns a `PipelineResult`. **Correctness invariants:** qubit
  Hamiltonian == FCI; Krylov converges to FCI and never dips below the variational floor;
  `dmrg_reference.reference_energy` shares the integral convention with
  `build_hamiltonian_from_integrals`.
* **Environment:** Everything runs in the `chem` conda env — always prefix shell commands with
  `conda run -n chem`. `ft_resource_estimator.py` (openfermion) runs in a separate `chem-ft` env.
* **Process isolation:** `block2`/DMRG segfaults if loaded into a process that already imported
  `pyscf` or `qiskit-aer` — DMRG and every `test_*_spec.py` gate run in their own process (the
  Makefile / `run_in_chem.sh` enforce this). Preserve it when adding tests.
* **Dependency pins:** `scipy>=1.8,<1.16` (1.16+ crash the qiskit-nature import chain). DataFrames
  use **polars**, not pandas.
* **Long jobs:** hour-long DMRG/benchmark runs die on session teardown if launched from here — the
  user runs those in their own terminal; block2 scratch goes to `.dmrg_tmp/` (git-ignored).

## 7. Definition of Done for a Discovery Spec

A spec is *done* (`[x]` in `BACKLOG.md`) when:
1. `make gates` is green and the gate could have failed (no vacuous passes — cf. the QCIVET /
   silent-exactness bugs found and fixed in the ODMD thread).
2. The headline claim is stated with its reference, its two-axis score (§2), and its honest
   boundary (what it is *not* — isolated cluster, no device noise, no advantage at scale, etc.).
3. The `BACKLOG.md` entry records the **finding**, including any boundary the work discovered —
   the mismatch is the science, not a failure.
