# SPEC: Certified two-sided energy brackets from the Krylov solve (Temple/Weinstein)

**Status:** IMPLEMENTED — gates G1–G4 green (`tests/test_temple_bracket_spec.py`).

---

## 1. Goal

Every energy this repo produces is a variational *upper* bound (QKSD Ritz values, PDS, DMRG);
nothing certifies from below. Claim: the QKSD ground Ritz eigenstate |Ψ₀(M)⟩ the solver already
builds yields, at the cost of **one extra expectation value** ⟨Ψ₀|H²|Ψ₀⟩, a rigorous **lower**
bound via Temple's inequality — so each Krylov solve carries a certified two-sided bracket
`[E_Temple, E_Ritz]` that provably contains the exact (sector) ground energy and whose width →
0 as M grows. Falsifiable at machine precision: a single M on a single system where the exact
energy escapes the bracket kills the claim.

## 2. Background and honest framing

- Classical bound theory: Temple (Proc. R. Soc. A 119, 276, 1928) — `E₀ ≥ ⟨H⟩ − σ²/(ε − ⟨H⟩)`
  for any `ε ≤ E₁` with `⟨H⟩ < ε`; Weinstein (1934) — `E₀ ≥ ⟨H⟩ − σ` when ⟨H⟩ lies closer to E₀
  than E₁. Modern revival for correlated chemistry: Pollak & Martinazzo (e.g. JCTC 15, 1498,
  2019; JCP 152, 244110, 2020). **Reproduction of known theory** — the composition with QKSD
  Ritz data (and the quantified certification overhead / self-consistency boundary) is what is
  new *for this repo*.
- **What we can claim if gates pass:** micro-Hartree-wide certified brackets around FCI from the
  standard Krylov solve + one ⟨H²⟩ evaluation, with the certification overhead (bracket width vs
  the uncertified Ritz error) measured; and an oracle-free variant whose validity region is
  mapped, not assumed.
- **What we cannot claim:** certification is **sector-restricted** — states built from |HF⟩ have
  exactly zero weight outside its symmetry sector, so E₀/E₁ mean the lowest/second reachable
  levels (QKSD itself sees nothing else). The Temple premise `ε ≤ E₁` is only *rigorous* with an
  oracle gap; the self-consistent mode `ε = θ₁ − σ₁` cannot verify its own premise (G4 maps
  where it holds empirically). Exact statevector; the shot cost of measuring ⟨H²⟩ on hardware
  (its Pauli expansion is ~λ²-sized) is NOT modeled — this is a classical-postprocessing
  certificate for the simulated pipeline.

## 3. Approach

From `QuantumKrylovSolver.eigenstates(M, n_states=2)` (validated in `SPEC_qksd_excited.md` /
`SPEC_qksd_properties.md`): electronic Ritz values θ₀, θ₁ and states Ψ₀, Ψ₁. One sparse matvec
each gives `σᵢ² = ⟨Ψᵢ|H²|Ψᵢ⟩ − θᵢ²`. Then

- upper = θ₀ (variational, already pinned by the test suite),
- Temple lower = `θ₀ − σ₀²/(ε − θ₀)` if `ε > θ₀`, else −∞ (vacuous, still valid),
- Weinstein lower = `θ₀ − σ₀`,
- ε: **oracle mode** — the exact reachable E₁ (dense diagonalization; validation only);
  **self mode** — `ε = θ₁ − σ₁`, a Weinstein-style E₁ estimate from the same Krylov data.

**Reference:** exact reachable spectrum by dense diagonalization (as `SPEC_qksd_excited.md`).

## 4. Public interface

```
temple_bounds.EnergyBracket           # dataclass: m, upper, lower, weinstein_lower, width,
                                      #   variance, eps, eps_source  (total energies, Ha)
temple_bounds.krylov_bracket(mh, m, eps=None, solver=None) -> EnergyBracket
    # eps: exact E1 as a TOTAL energy (oracle mode) | None -> self-consistent mode
    # solver: pass a shared QuantumKrylovSolver to reuse the cached Krylov basis
temple_bounds.bracket_ladder(mh, dims, eps=None) -> list[EnergyBracket]
```

Top-level module `temple_bounds.py` (a method rung, like `msd.py`/`odmd.py`).

## 5. Acceptance criteria (validation gates)

All in `tests/test_temple_bracket_spec.py`; systems H2, H4 chain, LiH, N2 CAS(6,6); depths
M ∈ {2, 4, 6, 8, 12, 16, 20, 24}. Noiseless (the claim is about rigor, not sampling).

- **G1 — containment, no exceptions (DEFINITION OF DONE).** For every system × M, oracle-mode
  `lower ≤ E_FCI ≤ upper` (tolerance 1e-9 Ha). One escape kills the spec.
- **G2 — the bracket closes.** width(M=16) < width(M=4) on every system, and
  width(M=16) < 1e-5 Ha on H4 and N2 (measured ~1e-8 and 3.8e-6), width(M=24) < 1e-4 Ha on LiH
  (measured 3.3e-5 — LiH's small sector gap of 133 mHa makes it the slow case).
- **G3 — certification overhead is small.** At mid-convergence (H4 M=8, N2 M=12) the certified
  width is < 5× the uncertified Ritz error (measured ~2.7× both) — a rigorous bracket costs
  under 5× the raw error, not orders of magnitude.
- **G4 — the self-consistent mode's validity region (the boundary).** (a) With ε = θ₁ − σ₁:
  `lower ≤ E_FCI` holds for every system at every M ≥ 6, and width(M=24) < 1e-4 Ha on all four.
  (b) At M = 4 the premise fails on H4 and N2 (ε > exact E₁ — the Krylov space has not resolved
  the excited state yet), so small-M self-certification is *not* rigorous: the recorded finding.

## 6. Implementation plan (test-first)

1. `tests/test_temple_bracket_spec.py` encoding G1–G4 (RED — `temple_bounds` missing).
2. `temple_bounds.py`: minimum code on top of `eigenstates()` — no new subspace machinery.
3. `make gates` to green.

## 7. Out of scope

- Lehmann/Kato optimal intervals and the Pollak–Martinazzo self-consistent ladder (sharper, more
  machinery — a follow-up).
- Shot-noise / hardware measurement cost of ⟨H²⟩ (λ²-sized Pauli expansion).
- Lower bounds for excited states; brackets under Trotterized evolution.
- Cross-sector certification (a different reference state per sector would be needed).

## 8. Caveats and risks

- **R1 — premise sensitivity:** Temple silently degrades (not breaks) if ε overshoots E₁ only
  slightly; G4(b) pins where that happens so the self mode is never quoted below M=6.
- The bracket certifies the *reachable-sector* ground energy: identical scope to QKSD itself,
  but state it when quoting.
- Vacuous −∞ lowers (ε ≤ θ₀, seen on LiH M=8–12 self mode) are valid but useless — callers must
  check `width` is finite, not just trust `lower`.

## 9. Deliverables

- `temple_bounds.py` — `EnergyBracket`, `krylov_bracket`, `bracket_ladder` (+ `__main__` table).
- `tests/test_temple_bracket_spec.py` — gates G1–G4.
- `BACKLOG.md` entry with measured widths/overheads and the validity boundary.
