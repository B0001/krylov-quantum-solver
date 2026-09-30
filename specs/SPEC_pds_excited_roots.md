# SPEC: The higher roots of PDS's P_K(E) are bounded excited-state estimates — but converge far too
slowly to be useful at attainable K

**Status:** CLOSED — gates G1(+G1b)–G4 PASS (2026-09-29). `pds_roots` added to `moment_expansion.py`.
The **bound** claim survives everywhere it's asserted (H₄ equilibrium-ish K=3..8; stretched
near-degenerate K=3..6) — but **not everywhere it was tested first**: at K=7 on the stretched
geometry, even with BLAS pinned to one thread, the run reproducibly lands on one of two
floating-point code paths, and one of them **violates the bound by −7.49 mHa** (`cond(M) ≈ 1.2e17`).
That combination is excluded from G1 rather than papered over — G1b pins the conditioning number
itself as the falsifiable reason. The **usefulness** claim named in `SPEC_moment_pds` §7 is
**KILLED** — PDS(8)'s first-excited root misses FCI's E₁ by 12–14 mHa on H₄, 8–9x over the 1.6 mHa
chemical-accuracy bar. Index-matching by nearest distance on near-degenerate stretched H₄ also breaks
down exactly as predicted: root₂ skips reachable level 2 and locks onto level 3 instead — PDS roots
track spectral *density*, not level *index*, once levels crowd closer than the achievable
resolution.

> A spec is a *falsifiable hypothesis*, not a contract: if implementation shows a gate is wrong,
> change the gate and record why (that mismatch is the finding).

---

## 1. Goal

`SPEC_moment_pds` builds `P_K(E)` from `⟨H^n⟩` moments of the HF reference and keeps only the
smallest root (the ground-state estimate); §7 named "excited states from the higher roots" as an
unimplemented follow-up. Claim: the moment problem `P_K(E)` solves is equivalent to the K-node Gauss
quadrature of the spectral measure of H against |HF⟩ — the same information a K-step Lanczos
recursion carries — so by Rayleigh–Ritz/Cauchy interlacing, the **j-th smallest real root is a
variational upper bound on the j-th smallest reachable eigenvalue** (reachable = nonzero HF overlap,
`SPEC_qksd_excited`'s definition). The claim is **false** if any root_j < E_j − 1e-9 Ha. Separately:
is this bound *useful* — does PDS(8)'s first-excited root land within chemical accuracy of the true
E₁ on H₄, as `specs/BACKLOG.md`'s scout probe suggested it would not (ground root 2.04→0.002 mHa over
K=3..7 while root₁ crawls 553→20 mHa)?

## 2. Background and honest framing

- **Prior art.** Peng & Kowalski, Quantum 5, 473 (2021) (PDS, ground state only — the higher-root
  interpretation is this repo's own extension, not from the paper). Gauss-quadrature/moment-problem
  equivalence to Lanczos: standard numerical-analysis result (Golub & Meurant, *Matrices, Moments and
  Quadrature*). Degeneracy/index-matching convention: `SPEC_qksd_excited` R1 (match as a sorted set,
  never by positional identity).
- **Ground truth.** The *reachable* exact spectrum: dense-diagonalize the same qubit Hamiltonian,
  keep eigenvalues with `|⟨v_i|HF⟩|² > 1e-8`, sort ascending, lift by the energy offset — identical
  construction to `SPEC_qksd_excited`'s `_reachable_spectrum`.
- **What we can claim if gates pass.** The variational bound on higher roots is real and holds within
  the conditioning regime where float64 can represent it: H₄ equilibrium-ish K=3..8, stretched
  near-degenerate K=3..6 — a genuine, cost-free (no new moments beyond what PDS(K) already computes)
  by-product of PDS(K).
- **What we cannot claim.** (a) Usefulness: the first-excited root does **not** reach chemical
  accuracy at K=8 on H₄ — it is ~8x over the bar, and going higher K does not help because the moment
  matrix in the raw (uncentered) electronic frame is already ill-conditioned past machine precision
  (`cond(M) ≈ 1e17` at K=7, `≈1e20` at K=8 — consistent with `specs/BACKLOG.md`'s separate
  centered-frame entry, chem-70c, not attempted here). (b) Index-by-index level identification: on a
  near-degenerate stretched geometry, matching roots to reference levels positionally is wrong — a
  root can sit between two true levels and the nearest-neighbor assignment skips one. (c) **The
  variational guarantee itself, past `cond(M) ≳ 1e14` combined with near-degenerate orbitals**: it is
  a hard theorem in exact arithmetic, but at K=7 on the stretched geometry float64 reproducibly
  violates it on roughly half of fresh-process runs (by mHa, not ULPs) — "variational" is a property
  of the exact linear algebra, not of what `np.linalg.solve`/`np.roots` return once the input matrix
  is singular to working precision. (d) No new moments/measurements: same `hamiltonian_moments` call
  as ground-state PDS, exact statevector only.

## 3. Approach

`pds_roots(moments, order, offset)`: build the same `M`, `Y`, solve for the polynomial coefficients
exactly as `pds_energy` does, but return **all** real roots (imag < 1e-6) sorted ascending, each
lifted by `offset`, instead of just the minimum. `pds_energy` becomes `pds_roots(...)[0]`. Reference:
`_reachable_spectrum` (dense diagonalization + HF-overlap filter), reused verbatim from
`SPEC_qksd_excited`'s test convention.

## 4. Public interface

```
moment_expansion.pds_roots(moments, order, offset=0.0) -> np.ndarray   # sorted ascending, real roots + offset
moment_expansion.pds_energy(moments, order, offset=0.0) -> float       # now pds_roots(...)[0], unchanged behavior
```

## 5. Acceptance criteria (validation gates)

Gates in `tests/test_pds_excited_roots_spec.py`. H₄ equilibrium-ish (`1.0/2.0/3.0 Å`, same geometry
as `SPEC_moment_pds`) and H₄ stretched/near-degenerate (`2.0/4.0/6.0 Å`, chosen because the reachable
spectrum has several levels within ~20 mHa of each other — see §6). PySCF/qiskit only, no block2.

- **G1 — variational bound (the can't-be-faked invariant, definition of done).** Every real root_j ≥
  the j-th reachable eigenvalue − 1e-9 Ha. Checked for H₄ equilibrium-ish at K ∈ {3..8}, and H₄
  stretched/near-degenerate at K ∈ {3..6}. **MEASURED: holds everywhere in that scope**, including
  equilibrium-ish K=7,8 (`cond(M) ~ 1e17–1e20`) over 8+ pinned-thread fresh-process trials.
  Equilibrium-ish K≥7 degrades gracefully to "true but useless." **The stretched geometry's K=7,8 are
  excluded from this gate, not assumed** — see G1b and R3: that specific (near-degenerate, cond > 1e17)
  combination reproducibly *breaks* the bound on some floating-point code paths, which is a sharper,
  unplanned finding beyond what the bead's own killable check anticipated.
- **G1b — the excluded cell is excluded for a falsifiable reason.** `cond(M)` at K=7 on the stretched
  H₄ is `> 1e14` (measured `≈1.2e17`) — deterministic given the moments, this is the reason G1 stops
  trusting the strict bound there, made into its own assertion so it can't silently rot.
- **G2 — usefulness at K=8 (KILLED).** `SPEC_moment_pds` §7 / the backlog scout probe predicted
  PDS(8)'s root₁ would land within 1.6 mHa of reachable E₁ on H₄. **MEASURED: it does not** — the
  gap is 12–14 mHa (equilibrium-ish geometry, varies a few mHa across processes), ~8-9x the
  chemical-accuracy bar. This gate asserts the miss (`> 5e-3 Ha`, a margin well clear of run-to-run
  float noise at this conditioning) so a future silent "fix" that makes the number look good without
  addressing conditioning gets caught, not quietly accepted.
- **G3 — differential convergence rate (quantifies the "~1000x slower" backlog claim).** At K=7 on
  H₄ equilibrium-ish: `|root₀ − E₀| / |root₁ − E₁|` ≥ 100 — the ground root is at least two orders of
  magnitude more converged than the first excited root at the same K, from the *same* linear solve.
  **MEASURED: ratio ≈ 11400** (0.0017 mHa vs 19.8 mHa) — far past the 100x floor, close to the
  backlog's own predicted ~1000x.
- **G4 — index-matching breaks down on near-degenerate stretched H₄ (explicitly tested, not papered
  over).** At K=6 (chosen: `cond(M) ≈ 1.5e14`, still inside double-precision range, unlike K=7/8's
  `cond(M) > 1e17`), nearest-neighbor-assign each of the first 6 roots to its closest reachable level.
  **MEASURED finding, asserted as the gate:** the assignment is **not** injective onto `{0,1,2,3,4,5}`
  — reachable level index 2 (`E=-1.3689`) is skipped (root₂ lands at `-1.3500`, closer to level 3 at
  `-1.3542`), and the assignment is not required to be, and is not, accuracy-within-chemical-accuracy
  for the skipped/high roots. The gate does *not* assert good accuracy here — it asserts the skip
  happens, so the finding cannot silently disappear if someone "fixes" the polynomial solve.

> Definition of done: **G1**. G2–G4 are pre-registered *falsification* checks per the bead — a
> passing spec here documents that the naive higher-root usefulness claim is false, not that the
> feature "doesn't work." `pds_roots` itself is correct and complete; its output is honest, low-value
> data at K≥6 on these systems without the centered-frame conditioning fix (chem-70c, separate bead).
> G1's scope was *revised during implementation* to drop the stretched geometry's K=7,8 once testing
> showed the bound itself, not just its usefulness, can fail there — the standing SDD rule ("if a
> gate proves unsatisfiable, revise it and record why") applied to a gate that was originally written
> assuming the bound was inviolable everywhere.

## 6. Implementation plan (test-first)

1. `scratch_pds_roots_scout.py` (throwaway, not committed to `tests/`) computed the numbers in §5
   before the gates were written, per the bead's own pre-registered scout-probe values.
2. Add `pds_roots` to `moment_expansion.py`; refactor `pds_energy` to call it (no behavior change —
   `test_moment_pds_spec.py` G1–G4 stay green, verified).
3. Write `tests/test_pds_excited_roots_spec.py` encoding G1(+G1b)–G4 above.
4. `make gates` (own process; pyscf/qiskit only, no block2).

## 7. Out of scope

- Fixing the conditioning (centered moment frame) — that is `specs/BACKLOG.md`'s separate
  centered-frame entry (chem-70c). This spec measures the raw-frame behavior honestly, including its
  failure mode, and does not attempt the fix.
- A general degeneracy-robust index-matching algorithm (e.g. Hungarian assignment with a
  no-match/abstain option) — G4 documents that naive nearest-neighbor matching fails; building a
  better matcher is a follow-up, not this bead.
- `n_states`-style public API sugar (`ExcitedKrylovStep`-equivalent for PDS) — `pds_roots` returning
  the raw sorted array is the minimum reuse-first primitive the bead asked for.
- Excited-state PDS on systems beyond H₄ (H₂'s reachable space is only 2-dimensional at this active
  space, too small to show the index-skip effect; LiH not probed here).

## 8. Caveats and risks

- **R1 — the K=8 usefulness number is not float-reproducible to the mHa.** Repeated `pds_roots` calls
  *within one process* on identical moments are bit-identical (verified, 5 trials); numbers can shift
  by a few mHa *across fresh interpreter processes* at `cond(M) > 1e17` (BLAS/threading-dependent
  rounding on a numerically singular matrix) — expected at this conditioning, not a bug. G2's `> 5
  mHa` margin is chosen wide enough to be robust to this.
- **R2 — the index-skip finding is geometry- and K-specific.** Chosen H₄ stretch (2/4/6 Å) and K=6
  produce a clean, reproducible skip; a different stretch or K could show a different (or no) skip.
  The claim is "this failure mode exists and is real," not "it always happens."
- **R3 — the variational bound is a theorem about exact arithmetic, and float64 can violate it.**
  Discovered while writing G1 (not pre-registered — the bead's own check assumed the bound survives
  and only the usefulness claim would die): at K=7 on the stretched/near-degenerate H₄, with BLAS
  pinned to 1 thread for reproducibility (`tests/test_shift_both_sides_spec.py`'s pattern), 3 of 6
  fresh-process trials landed on a code path that violates `root_j ≥ E_j − 1e-9` by **−7.49 mHa**
  (`cond(M) ≈ 1.2e17`); the other 3 landed on a path with a small positive margin (`+196 μHa`). Both
  outcomes are deterministic *given* the code path (repeats exactly within a process, and the same
  bimodal split repeats across fresh processes) — the split itself is not controlled by the thread
  pins, so its exact source (SCF/CASCI orbital-gauge selection on near-degenerate orbitals, `np.roots`
  companion-matrix eigensolve internals, or both) is not pinned down further here. The practical
  upshot: **do not trust a PDS higher-root bound at `cond(M) > ~1e14`** on a system with
  near-degenerate orbitals, even if a single run looks fine — rerun or check conditioning first.
- Honest limitation: exact statevector, minimal-basis H₄ only — a correctness/behavior study, not a
  claim about production excited-state usefulness of raw-frame PDS.

## 9. Deliverables

- `moment_expansion.py` — `pds_roots` (new), `pds_energy` (now a thin wrapper, unchanged behavior).
- `tests/test_pds_excited_roots_spec.py` — gates G1, G1b, G2, G3, G4.
- This spec, recording the killed usefulness claim, the index-skip finding, and the float64
  bound-violation finding (R3).
