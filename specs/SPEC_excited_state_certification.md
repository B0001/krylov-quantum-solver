# SPEC: Certified lower bounds on the first excited state — and why they cannot be self-certified

**Status:** v2.0 — IMPLEMENTED. The headline claim (a *self-certified* excited-state solver) is
**FALSIFIED**; two of the four v1.0 gates were unsatisfiable as written. v1.0 is preserved verbatim
in §8 with the measurement or argument that resolved each claim.
**Depends on:** `temple_bounds.py`, `certified_gaps.py`, `reachability.py`,
`hybrid_quantum_solver/quantum_krylov_solver.py`
**Gates:** `tests/test_excited_bounds_spec.py` (G1, G1b, G1c, G2, G2b, G3, G4, G4b — 10 tests, 98 s)
· `excited_bounds.py` · `data/excited_certification_bench.csv`

---

## 1. Goal

`temple_bounds.py` certifies E₀ from below. `certified_gaps.py` certifies the gap, but its floor for
E₁ is `ε₁ = θ₁ − σ₁`, which its own docstring calls a *premise, unverifiable*. This spec asks
whether the first excited state can be given a genuine two-sided bracket from the same real-time
Krylov data, and whether the separator that such a bound needs can be obtained **without a
classical oracle** — the "self-certified" mode that was v1.0's entire headline.

**Answers.** Yes to the bracket, with an oracle separator: it is rigorous, it contains E₀ and E₁
at every gated depth on three molecules, and it costs one extra `H|u⟩` per Ritz vector. **No to
self-certification** — not "not yet", but *never from this data*: the separator requires an
eigenvalue **count**, a Krylov subspace cannot certify a count, and a level with small HF overlap
is invisible to every quantity the subspace produces. §2 and G4b give the construction.

### Reproduction vs novelty (required up front by `specs/README.md`)

**Reproduction, and textbook:** the bounds themselves. Lehmann's optimal subspace lower bounds
(Lehmann 1949/1950; Maehly 1952; modern treatments Beattie & Goerisch, *Numer. Math.* 72, 143
(1995); Zimmermann & Mertins, *ZAMM* 75 (1995)) and their numerical use are a settled literature,
as are Temple (1928), Weinstein and Kato. Nothing here is a new inequality.

**Novel here:** only the *composition* — Lehmann bounds evaluated on a real-time QKSD basis inside
this repo's certified arc — and the *negative result* on self-certification, which is a statement
about this construction, not about the literature. v1.0's claims to the **"world's first
self-certified excited-state solver"** and to closing **"a major gap in the theoretical computer
science of spectral certification"** are withdrawn: they were false as novelty claims, and the
thing they named is the thing this spec falsified.

---

## 2. What v1.0 assumed, and what is actually true

**v1.0 §3.2 wrote down the wrong inequality for the theorem it invoked.** §2 describes Lehmann's
multi-dimensional generalisation correctly; §3.2 then writes
`E_i ≥ θ_i − σ_i²/(β − θ_i)`, which is the **single-vector Temple/Kato** bound applied per Ritz
vector. This implementation uses **Lehmann's pencil** (§3), keeps the per-vector Temple number
alongside it for comparison (`temple_lower_1`), and gates that the two differ in the direction the
theory requires (G1b: Lehmann is tighter, by +0.067 mHa on LiH at M = 4).

**v1.0 never stated the precondition the excited-state bound needs.** The ground-state Temple bound
is safe because θ₀ < E₁. For E₁ the requirement is:

> **β must satisfy θ₁ < β ≤ E₂, where E₂ is the third eigenvalue of the HF-*reachable* sector.**

Derivation: for any state |u⟩ in the reachable sector, `⟨u|(H − E₁)(H − β)|u⟩ = Σ_λ |c_λ|²(λ − E₁)(λ − β) ≥ 0`
provided no reachable eigenvalue lies strictly between E₁ and β, i.e. β ≤ E₂. Expanding and
dividing by β − θ > 0 gives the bound. Three consequences v1.0 missed:

1. The bound is valid for **any** state once β is valid — the Ritz vector does *not* have to be
   dominated by the right eigenvector. All of the risk lives in β, and none of it in the vector.
2. Because the Krylov space lies inside the reachable sector, E₁ and E₂ are the sector's
   eigenvalues; unreachable levels in between are harmless. This is the same sector restriction
   `certified_gaps.py` and `reachability.py` already carry.
3. **A degenerate E₁ admits no β at all.** θ₁ ≥ E₁ by Poincaré separation, so E₁ = E₂ forces
   θ₁ < β ≤ E₁ ≤ θ₁. The bound is then vacuous by construction, at any M, in any basis. This is
   exactly what happens to Be₂ (§6.5).

**And the system v1.0 picked for its convergence gate has no separator at all.** STO-3G H₂'s
reachable sector holds **two** levels: E₂ does not exist, E₁ is the top of the sector, and the
Krylov space saturates the sector at rank 2 from M = 2 with θ_i − E_i = ±1e-15. "Micro-Hartree at
M = 12" is not a hard target there, it is a category error (G2).

---

## 3. Approach (as built)

1. **Ritz data.** `QuantumKrylovSolver.eigenstates(M, n_states=3)` gives orthonormal Ritz vectors
   |u₀⟩, |u₁⟩, |u₂⟩ and their Ritz values. One sparse `H|u_i⟩` per vector supplies both θ_i and
   σ_i² (`temple_bounds.mean_and_variance`, reused, including its clip of the negative-rounding
   variance).
2. **Upper bounds, free.** θ_i ≥ E_i by Poincaré separation — the same interlacing fact
   `solve_excited` already documents.
3. **Lower bounds: Lehmann's pencil.** On the span of the Ritz vectors with θ_i < β, form
   `A₁ = W†(H − β)W` (negative definite by that selection — Lehmann's hypothesis) and
   `A₂ = W†(H − β)²W`, and solve `A₂ y = τ A₁ y`. Sorted descending, `β + τ` bounds the eigenvalues
   immediately below β: `β + τ₍₁₎ ≤ E₁`, `β + τ₍₂₎ ≤ E₀` **when exactly two eigenvalues lie below
   β**. At d = 1 this reduces algebraically to Temple; at d ≥ 2 it is strictly tighter.
4. **Refusal, not guessing.** If no Ritz value lies below β, or the selection leaves `A₁`
   indefinite, the function returns an empty array and the bracket carries `−inf`. A `−inf` lower
   bound is a valid, vacuous bound: callers must check `isfinite` before quoting a width.
5. **Self-certified mode** (v1.0's `β = θ₂ − σ₂`) is implemented and is **labelled an estimate**.
   It exists so the gates can measure how it fails, not because it is usable.

---

## 4. Public interface (`excited_bounds.py`)

```python
lehmann_lower_bounds(H, states, beta) -> np.ndarray     # descending; [] = refusal
bracket_from_states(H, states, beta=None, m=0, offset=0.0) -> ExcitedBracket   # matrix level
lehmann_excited_brackets(mh, krylov_dim, beta_oracle=None, solver=None) -> ExcitedBracket
bracket_ladder(mh, dims, beta_oracle=None, solver=None) -> list[ExcitedBracket]
reachable_separator(mh) -> float        # REFERENCE ONLY (dense): E_2 of the reachable sector
```

`ExcitedBracket` keeps v1.0's field names (`lower_0/upper_0/lower_1/upper_1`, `best_estimate_*`,
`gap_lower`, `gap_upper`) and adds `temple_lower_1`, `beta`, `beta_source`, `sigma1`, `rank`.

**Three deviations from v1.0's signature, all deliberate:**

- `beta_oracle` defaults to `None` = *self* mode in v1.0's design, and that default is now the
  **unsafe** one. It is kept (so the gates can exercise it) but documented as an estimate, and
  `reachable_separator` is provided so the safe call is as short as the unsafe one.
- A **matrix-level** entry point (`bracket_from_states`) was added. Without it the missed-level
  witness in G4b could not drive the module's own decision logic, and a falsification that cannot
  be run against the shipped code is not a gate.
- The function is named for Lehmann's *bounds*, not "Lehmann enclosure": the upper half of the
  bracket is Poincaré interlacing, not Lehmann's theorem.

---

## 5. Acceptance gates (as gated, with the measured numbers)

- **G1 — Simultaneous containment, oracle separator.** LiH CAS(2,5) / H₄ / N₂ CAS(6,6), M ∈
  {4,6,8,10,12,16,20,24}: L₀ ≤ E₀ ≤ U₀ and L₁ ≤ E₁ ≤ U₁ at **24/24** brackets, zero escapes. Finite
  (non-vacuous) lower bounds: LiH 8/8, H₄ 7/8, N₂ 5/8 — the refusals are exactly the depths where
  θ₁ has not yet fallen below E₂, and the gate requires ≥ 5 finite per system so it cannot be
  passed by refusing. Tightest margins E₀ − L₀ = 1.1e-11 Ha and E₁ − L₁ = 1.4e-6 Ha (H₄, M = 20)
  against run-to-run scatter 9e-15 / 9e-9. *Passes.* (v1.0 specified H₂ and LiH; H₂ is
  uncertifiable — G2 — so it was replaced by H₄ and N₂.)
- **G1b — It is Lehmann, not the per-vector Temple bound of v1.0 §3.2.** Lehmann ≥ Temple at every
  finite point; best gain **+6.7e-5 Ha** (LiH, M = 4). Compared at 1e-6 Ha because the pipeline's
  own arithmetic floor is ~1e-8 Ha (§7). *Passes.*
- **G1c — The pencil is a bound, not an estimate.** 300 seeded random Hermitian 8×8 with perturbed
  eigenvector subspaces and every separator position: worst escape **exactly 0.0**, tightest bound
  within **1.6 mHa** of its eigenvalue (so containment is not vacuous). Misstating the eigenvalue
  count by one — the G4 failure mode in miniature — produces **510** violations. *Passes.*
- **G2 — THE FIRST KILL (replaces v1.0's G2).** STO-3G H₂'s reachable sector holds exactly **2**
  levels ⇒ no E₂, no separator, `lower_1 = −inf` at **every** M; rank saturates at **2** from M = 2
  and |θ₁ − E₁| ≤ 1e-13. v1.0's "monotonically … micro-Hartree at M = 12 on STO-3G H₂" is
  unsatisfiable in both halves. *Passes as the kill.*
- **G2b — Convergence, re-gated where it is measurable.** U₁ − L₁ (mHa): LiH 88.26 → 89.01 → 90.06
  → 3.89 → 3.98 → 4.19 → **0.0256** → 0.0277; H₄ 981 → 215.7 → 113.3 → 15.77 → 6.98 → 0.752 →
  **0.0017** → 0.0033 over M = 4…24. Closes **> 3 orders** from M = 6 to M = 20 on both — and is
  **not monotone**: LiH widens at five of the seven steps (M = 6, 8, 12, 16, 24), H4 at M = 24. *Passes (as a ceiling plus an explicit
  non-monotonicity assertion).*
- **G3 — v1.0's "< 5×" is true on one system of three.** Overhead = certified gap-bracket width ÷
  |Ritz gap − exact gap|: LiH **1.23–1.30×**, H₄ **2.30–5.42×**, N₂ **2.58–11.68×**. Gated as
  `< 1.5×` on LiH, `> 5×` on H₄ and N₂ (the falsification), `< 12×` overall ceiling, and `> 1.0×`
  everywhere. *Passes.*
- **G4 — THE SECOND KILL.** Self mode (β = θ₂ − σ₂) escapes on **real molecules**: LiH M = 4, 6, 8
  (E₁ under-bounded by **2.02, 1.97, 1.89 mHa**; E₀ also escapes by 2.0 µHa at M = 4) and N₂
  CAS(6,6) M = **16** (by **17.4 mHa** — well past v1.0's "all M ≥ 8"). Every escape has
  β_self > E₂, and the converse fails: N₂ at M = 4 has an unsound β and survives anyway, so a
  passing self-mode run is not evidence. *Passes as the kill.*
- **G4b — Why it cannot be repaired.** Constructed witness: levels {0, 0.5, 1, 2, 3} with amplitude
  **1e-4** on the 0.5 level (population 1e-8, **above** the certified arc's 1e-10 reachability cut,
  so E₁ = 0.5 genuinely belongs to the sector). At M = 4 the Krylov space returns Ritz values
  [0, 1, 2, 3] — the level is simply absent — with **σ₁ = 6.8e-5** reporting a converged subspace.
  Self mode picks β = 1.9999 and certifies **E₁ ≥ 1.000000**: an overshoot of **0.500 Ha**. With
  the true separator β = E₂ = 1.0 the precondition survives only by 3e-9 (θ₁ *is* the level the
  subspace mistook for E₁), the pencil is correspondingly ill-conditioned, and the bound collapses
  onto E₁ ± 1e-5 — it stops lying and stops saying anything. *Passes.*

---

## 6. The findings

1. **The excited-state bracket works, given a separator.** One extra `H|u⟩` per Ritz vector buys a
   rigorous two-sided bracket on E₁ that closes by three orders of magnitude with depth and
   contains E₀ and E₁ at every gated point across three molecules (G1, G2b). This is the piece
   `certified_gaps.py` was missing: it replaces that module's unverifiable `θ₁ − σ₁` floor with a
   bound whose premise is explicit and checkable.
2. **Certification costs O(1) in width, and v1.0's "5×" was a guess.** The certified interval is
   1.23–1.30× the (unknowable) raw Ritz gap error on LiH, 2.30–5.42× on H₄, 2.58–11.68× on N₂. The
   *reason* the overhead is small — bracket width and Ritz error both scale as σ²/gap — also says it
   is not system-independent, which is exactly how a single number like "5×" went wrong.
3. **Self-certification is impossible from Krylov data, and the reason is a counting argument.**
   Lehmann's theorem needs to know how many eigenvalues lie below β. `β = θ₂ − σ₂` is an estimate:
   Kato's interval guarantees only that *some* eigenvalue lies within σ₂ of θ₂, never that it is
   E₂. When the subspace has not resolved level 2, three eigenvalues sit below β instead of two,
   and the pencil's top root — a true bound on E₂ — is returned as a bound on E₁ (G4).
4. **The impossibility is structural, not a depth problem.** A reachable level with small HF
   overlap is invisible to *every* quantity the subspace produces, so no function of that data can
   supply the count. The witness in G4b hides a level at amplitude 1e-4 and the "certified" bound
   is 0.5 Ha wrong while every residual reports convergence (G4b). Same failure class as
   `SPEC_subspace_floor_resolvability`'s ~1e-4-amplitude level near the cluster boundary and
   `SPEC_reachability_tolerance`'s two thresholds selecting different ground states — one rung up:
   there it corrupted a floor and a target, here it corrupts a *certificate*.
5. **Be₂ — v1.0's showcase — is out of scope for a structural reason worth recording.** At
   CAS(4,8)/cc-pVDZ (R = 2.45 Å) the lowest excited singlets are a **degenerate π pair**,
   E₁ = E₂ = −29.085033 Ha, and §2.3 shows a degenerate E₁ admits no separator at any M or basis.
   Independently, those states are not HF-reachable by real-time evolution: the Krylov rank
   saturates at **3** and θ₁ sits **147 mHa above** E₁^CASCI at every M ≤ 14, i.e. the optically
   bright transition v1.0 asked to certify is not in the subspace at all. Be₂ is kept in the
   benchmark as a refusal row rather than deleted.
6. **Where the bound refuses is informative.** With an oracle β the method returns −inf exactly
   when θ₁ has not yet dropped below E₂ (N₂ at M ≤ 8, H₄ at M = 4). That refusal is the honest
   behaviour the self-certified mode lacks: given a *correct* separator the construction never
   lies, it only declines.

---

## 7. Honest scope and caveats

- **Sector-restricted**, like every certificate in this arc: E₀/E₁/E₂ are the lowest HF-reachable
  levels (`reachability.reachable_eigenpairs`, the `1e-10` population cut whose defects
  `SPEC_reachability_tolerance` documents). QKSD's own scope.
- **The precondition is θ₁ < β ≤ E₂ over that sector, and nothing weaker.** Degenerate E₁ ⇒ no β
  exists. β from an oracle means an independent *lower* bound on E₂ (CASCI, DMRG, an experimental
  level); a *point estimate* of E₂ is not a separator.
- **Self mode is not a certificate.** It is shipped labelled as an estimate only so the gates can
  measure its failure. Do not quote `beta_source == "self"` numbers as bounds.
- **`−inf` is a deliberate refusal**, valid but vacuous — check `isfinite` before quoting a width,
  and note a refusal is *not* evidence that the bracket would have been wrong.
- **Arithmetic floor ~1e-8 Ha.** The Ritz vectors of a thresholded canonical orthogonalisation are
  orthonormal only to ~1e-9 (measured on H₄), the pencil divides that defect by β − θ₁ ≈ 0.3 Ha,
  and the default `dt` comes from ARPACK with a random start vector, so L₁ moves by up to ~1e-8 Ha
  between identical runs. The bound is rigorous in exact arithmetic; this implementation is not a
  certificate below ~1e-8 Ha, and a vanishing β − θ₁ makes that worse (G4b's witness: β − θ₁ = 3e-9
  leaks 1e-6 Ha). Gated containment margins exceed the floor by 3+ orders, which is why G1 stands.
- **Exact statevector.** The hardware shot cost of ⟨H²⟩ (a ~λ²-sized Pauli expansion) is not
  modelled and shot noise on H/S is not propagated into the bound. `certified_noise.py` is the
  pattern for doing that; it is not done here.
- **Does not plug into `certkit_bridge.py`.** certkit v0.2.0's only claim kind is
  `lambda_min_enclosure` and its Temple rule (jn1.1) takes one vector plus a β. A two-vector pencil
  claim about the **second** eigenvalue has no rule there, so an emitted certificate would be
  unchecked by construction — which is precisely what the CI gate exists to prevent. Adding a
  Lehmann rule to certkit is the follow-up; until then these bounds stay out of `certkit_out/`.
- **Three molecules, one geometry each**, plus a 5-level synthetic witness. The counting argument
  behind G4/G4b is system-independent; the widths and overheads in G2b/G3 are not.

---

## 8. v1.0, and what resolved each claim

> **§1 (v1.0):** *"This closes a major gap in the theoretical computer science of spectral
> certification by providing the **world's first self-certified excited-state solver** for
> QMA-complete Hamiltonians."*

Withdrawn on both halves. Lehmann/Temple/Weinstein bounds and their numerical application are a
settled literature (§1, "reproduction vs novelty"); and the self-certified mode the sentence names
is falsified by G4/G4b. What is left is a composition and a negative result.

> **§2 (v1.0):** *"the first excited state E₁ … is left entirely uncertified, with zero bounds
> available in the published literature."*

False as stated: Lehmann's 1949/1950 bounds *are* bounds on excited states, which is why §2 could
cite them. The accurate version is narrower — E₁ was uncertified **in this repo**, and
`certified_gaps.py` said so in its own docstring.

> **§3.2 (v1.0):** *"If a separator interval [a,b] is known to contain exactly k exact eigenvalues
> … `E_0 ≥ θ_0 − σ_0²/(β − θ_0)`, `E_1 ≥ θ_1 − σ_1²/(β − θ_1)`"*

The formula is per-vector Temple/Kato, not Lehmann (§2). Implemented both; gated that Lehmann is
the tighter one (G1b, +6.7e-5 Ha on LiH). The precondition the prose skipped — θ₁ < β ≤ E₂ over the
reachable sector, and no separator at all for degenerate E₁ — is derived in §2 and is what makes
Be₂ uncertifiable (finding 5).

> **G2 (v1.0):** *"The excited state bracket width (U₁ − L₁) must close cleanly and monotonically
> as M increases, reaching micro-Hartree scales at M=12 on STO-3G H2."*

Unsatisfiable twice over. STO-3G H₂'s reachable sector has **2** levels, so no separator exists and
U₁ − L₁ = ∞ at every M; and the Krylov space saturates that sector at rank 2 from M = 2 with Ritz
values exact to 1e-15, so M = 12 carries no more information than M = 2. Re-gated on LiH/H₄ (G2b),
where the width closes by > 3 orders — and **not monotonically**: it widens at five of seven steps
on LiH and at the last step on H₄.

> **G3 (v1.0):** *"The certified gap bracket width must be less than 5× the raw uncertified Ritz
> gap error, proving that Lehmann certification is computationally efficient."*

Measured: **1.23–1.30×** (LiH), **2.30–5.42×** (H₄), **2.58–11.68×** (N₂ CAS(6,6)). True on one
system of three. Re-gated as the measured ceiling (< 12×) with the falsification pinned explicitly
(> 5× on H₄ and N₂). The surviving statement is stronger than it looks — the overhead is O(1) at
all, which variance-based intervals usually are not.

> **G4 (v1.0):** *"The self-certified mode (estimating β from θ₂ − σ₂) must be mathematically
> rigorous and pass all containment gates for all dimensions M ≥ 8 on H2, establishing a clear,
> machine-checked domain of validity."*

**False**, and the "domain of validity" cannot exist. It fails on real molecules inside the claimed
domain (N₂ CAS(6,6) at M = 16, E₁ under-bounded by 17.4 mHa; LiH at M = 4–8 by ~2 mHa), and the
construction in G4b shows why no depth threshold can rescue it: a reachable level with 1e-4
amplitude is invisible to the subspace, so the eigenvalue count β depends on is not a function of
the data. The gate that replaced it demonstrates the failure and pins the sound condition
(β ≤ E₂, supplied from outside) instead.

> **§6.3 / §7 (v1.0):** *"Test and verify on the Be₂ dimer to output a certified optical gap
> curve."*

Not achievable, for two independent structural reasons (finding 5): E₁ = E₂ at CAS(4,8) leaves no
separator, and the bright π states are not HF-reachable by real-time evolution. The deliverable
`data/excited_certification_bench.csv` therefore certifies LiH / H₄ / N₂ and keeps Be₂ as a
documented refusal row — the spec was revised rather than the system quietly swapped
(`specs/README.md` step 5).
