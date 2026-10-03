# SPEC: Where ODMD repairs — and where it silently breaks — the Temple bracket's oracle-free premise

**Status:** IMPLEMENTED — gates pre-registered 2026-10-02 (commit 8eb7978) before the gated run.
G-a and G-c pass as pre-registered; two G-b sub-predictions were **KILLED** by the run and revised
to pin what was measured — (b2) escapes are not confined to K ≤ 8, (b4) the flag's coverage is
not set by M (§5 notes, §10). Backlog hypothesis: specs/BACKLOG.md "Method rungs", *"A
non-variational estimator cannot supply the Temple premise — ODMD repairs it at shallow M, then
silently breaks it"* (chem-z3h).

---

## 1. Goal

`SPEC_temple_bracket`'s self mode feeds ε = θ₁ − σ₁ into Temple's inequality as a stand-in for E₁,
and its G4(b) records that the premise ε ≤ E₁ fails at M = 4. The backlog claim: ODMD's E₁
estimate (`odmd.py`), read off the first row of the overlap matrix S the Krylov solve already
measured, repairs that premise at shallow M with **no extra measurement** — but ODMD is
non-variational, so no depth guarantees ε ≤ E₁. Deliverable: a **measured K-vs-validity map** —
per system, per base Krylov depth M and ODMD signal depth K: does the substituted premise hold,
does the Temple bracket still contain E₀, and can anything short of an oracle tell. Falsifiable:
G-a kills the repair as vacuous where self mode is not actually violated; G-b fails if an ODMD
overshoot is never located or is absorbed rather than recorded INVALID; G-c fails if the free
(K ≤ M) region does repair the premise.

## 2. Background and honest framing

- Temple (1928): E₀ ≥ ⟨H⟩ − σ²/(ε − ⟨H⟩) for any ε ≤ E₁ with ⟨H⟩ < ε. `SPEC_temple_bracket` R1:
  a slightly-overshooting ε *degrades* the bound before it *breaks* it — so "premise violated" and
  "bracket wrong" are different events, and only the second is a false certificate.
- ODMD (`SPEC_odmd.md`) is explicitly not variational (its G2). `SPEC_odmd_uq` G4: resampling one
  signal cannot see model bias — so its error bars cannot flag an overshoot either.
- **"No extra measurement" holds only for K ≤ M.** The real-time Krylov S is Toeplitz,
  S_ij = s_{j−i} with s_k = ⟨φ₀|e^{−ikΔtH}|φ₀⟩: an M-dim solve has measured s₀…s_{M−1} and nothing
  more. An ODMD signal of length K > M costs K − M extra overlap measurements.
- **Provenance — what was seen before this spec was written.** `odmd_temple_repair_map.py`
  (snapshot 5c5c6db) and a research pass for this batch unit reported (scout, not gated):
  self mode at M = 4 overshoots E₁ by +15 mHa (H₄) / +22 mHa (N₂) while its bracket still contains
  E₀ (by 0.45 / 1.64 mHa); ODMD at K = 4 overshoots by +950 / +1530 mHa and the bracket escapes
  (H₄ at M ∈ {2, 4, 6}; N₂ at K ∈ {4, 5}, largest escape +1.259 mHa); premise violated without
  escape at H₄ K = 6–14, 28 and N₂ K = 6–10, 18–28 (M ∈ {2, 4}); repair at M = 4 first at K = 16
  (H₄) / K = 12 (N₂); LiH CAS(2,5) self mode not violated at M = 4 (−110 mHa); full-space LiH
  violated (+466 mHa) and escapes at K = 4, 6, 8; ε > θ₁(M) fires on the M = 4 escapes but not at
  M = 2. The backlog scout's H₄ numbers (ε = −1.59181, ODMD K=16 −1.653511, K=20 −1.649481 Ha)
  do **not** reproduce at the gated 0.9 Å H₄ chain. The gates below encode those observations as
  predictions with thresholds fixed here; anything not in this list (the flag at M ∈ {6, 8}, LiH's
  escape depths in M, its repair/degrade sets) is first measured by the gated run.
- **What we can claim if the gates pass:** a measured map on five STO-3G system variants (noiseless)
  of where the ODMD substitution keeps, degrades, or breaks the Temple certificate; and the measured
  coverage of one oracle-free, provably sound (but one-sided) premise-failure flag.
- **What we cannot claim:** a certificate of any kind from the substitution; anything about K < 8
  as ODMD-the-method (below `SPEC_odmd`'s validated depths — this is the substitution's failure, not
  ODMD's); the centered frame (ODMD runs here in the solver's *uncentered* frame with `solver.dt`,
  so K-thresholds do not transfer to `SPEC_odmd`'s τ = π/W); shot noise; K or M off the grid.

## 3. Approach

- **Reference (oracle, validation only):** exact HF-reachable E₀, E₁ (total Ha) from
  `reachability.reachable_eigenpairs(mh, tol=1e-8)`.
- One shared `QuantumKrylovSolver(mh)` per system (default Δt = π/width, threshold 1e-10).
- **Self mode:** `krylov_bracket(mh, M)` for M ∈ {2, 4, 6, 8}.
- **ODMD ε(K):** second retained DMD eigenphase of `odmd_spectrum(S[0, :K], solver.dt)` (+ offset)
  from the solver's own S at depth K. Independent of M. `None` if the fit resolves < 2 modes —
  recorded, not dropped.
- **Substituted bracket:** `krylov_bracket(mh, M, eps=ε(K))` at every grid cell (M, K).
- **Classification**, TOL = 1e-9 Ha (`SPEC_temple_bracket` G1's containment tolerance) throughout:
  - *premise violated* ⇔ ε > E₁ + TOL;
  - *escape* (the only INVALID verdict — a false certificate) ⇔ lower > E₀ + TOL;
  - *degraded* ⇔ premise violated but no escape;
  - *free* ⇔ K ≤ M;
  - *repaired* (at an M where self mode's premise is violated) ⇔ ODMD premise holds and the lower
    bound is finite.
- **The oracle-free flag (the "bracket invalid" decision).** Ritz values bound eigenvalues from
  above (Courant–Fischer / Hylleraas–Undheim–MacDonald, within the reachable sector): θ₁(M) ≥ E₁.
  So **ε > θ₁(M) proves ε > E₁** with no oracle, and the Temple lower bound is then not a
  certificate. One-sided by construction: ε ∈ (E₁, θ₁] passes undetected, and self mode
  (ε = θ₁ − σ₁ ≤ θ₁) can never trip it — which is `SPEC_temple_bracket`'s "self mode cannot verify
  its own premise" restated. Added to `temple_bounds.EnergyBracket` as a defaulted field + property;
  every existing caller is unchanged.

## 4. Public interface

```
temple_bounds.EnergyBracket.theta1            # NEW defaulted field (total Ha): 2nd Ritz value,
                                              #   inf when the subspace is rank-1
temple_bounds.EnergyBracket.premise_refuted   # NEW property: eps > theta1 + 1e-9 (oracle-free,
                                              #   sound, one-sided)
odmd_temple_repair_map.odmd_eps1(mh, solver, k) -> float | None
odmd_temple_repair_map.k_vs_validity_map(key) -> dict  # e0, e1, n_reachable (oracle);
                                              #   self {M: EnergyBracket}; odmd {K: eps | None};
                                              #   cells {(M, K): EnergyBracket}
odmd_temple_repair_map.summarize(map) -> dict # the recorded form (§10)
uv run python odmd_temple_repair_map.py       # full table + summary for all five systems
```

## 5. Acceptance criteria (validation gates)

All in `tests/test_odmd_temple_repair_map_spec.py`. Systems (STO-3G): `h2` H₂ 0.74 Å; `h4` H₄
chain 0.9 Å; `lih_cas` LiH 1.6 Å CAS(2,5) (`certified_gaps`' variant); `lih` LiH 1.6 Å full space
(`SPEC_temple_bracket`'s variant); `n2` N₂ 1.1 Å CAS(6,6). M ∈ {2, 4, 6, 8};
K ∈ {4, 5, 6, 7, 8, 10, 12, …, 28}. Noiseless.

- **G-a — is the repair non-vacuous? (self mode at M = 4).** The backlog claim is "violated on all
  three reference systems". Pre-registered prediction, which **kills that claim as stated**: the
  self-mode premise is violated (ε − E₁ > TOL) on `h4`, `n2`, `lih`; NOT on `lih_cas`
  (ε − E₁ < −TOL); and vacuous on `h2` — its reachable sector has ≤ 4 levels, so M = 4 exhausts it
  and ε ≤ E₁ + TOL. Count: 2 of 3 under `SPEC_odmd`'s set (H₂, H₄, N₂) and under `certified_gaps`'
  set (H₄, LiH CAS(2,5), N₂); 3 of 4 under `SPEC_temple_bracket`'s (H₂, H₄, LiH, N₂).
- **G-b — the K-sweep locates the overshoot, and the bracket is recorded INVALID, not absorbed
  (DEFINITION OF DONE).**
  - (b1) On `h4`, `n2`, `lih`: at least one K overshoots (ε(K) > E₁ + TOL); on `h4` and `n2` the
    K = 4 estimate overshoots by > 100 mHa.
  - (b2) On each of `h4`, `n2`, `lih` the map records ≥ 1 containment escape, its worst escape
    > 1e-6 Ha (1000× TOL — not float noise); every escape on the grid sits at K ≤ 8; the largest
    escape over all systems is > 1 mHa.
    *Revised after the gated run:* the K ≤ 8 clause is **KILLED** by full-space LiH, which escapes
    at M = 6, K ∈ {12, 14, 16} (+0.106 / +0.118 / +0.031 mHa) — ODMD still overshoots E₁ there by
    250–430 mHa (27 reachable levels in a 133 mHa gap: its 2nd DMD mode is unresolved through
    K = 18). The gate now pins both halves: H₄, N₂ and LiH CAS(2,5) escape only at K ≤ 8;
    full-space LiH escapes at some K > 8.
  - (b3) Premise violation ≠ broken bracket: on `h4` and `n2` there are *degraded* cells (premise
    violated, E₀ still contained) — so INVALID must be decided by containment, not by the premise.
  - (b4) The oracle-free flag: (i) **sound** — on every cell of all five systems, `premise_refuted`
    ⇒ ε > E₁ + TOL, and it never fires in self mode; (ii) it fires on every escape with M ≥ 4;
    (iii) it fires on no escape at M = 2 — the recorded coverage gap (θ₁ of a 2-dim subspace sits
    above the overshoot).
    *Revised after the gated run:* (ii) and (iii) are **KILLED both ways** — the flag misses
    full-space LiH's escapes at (M, K) = (4, 4), (4, 5), (6, 16) and catches LiH CAS(2,5)'s at
    (2, 4), (2, 5). Coverage is decided by how far ε lands above the *Ritz* E₁, ε − θ₁(M), not by
    M. The gate now pins (i) soundness and that on every system with escapes the flag catches
    some and misses some (37 of 45 caught, §10).
- **G-c — the free region does not repair.** On `h4` and `n2` at M = 4 (the self-mode failure
  depth): the only free ODMD depth (K = 4) escapes, while the self-mode bracket at M = 4 still
  contains E₀ — the free substitution makes the certificate *worse*; and the smallest repairing K
  at M = 4 exceeds 2M = 8 — repair costs more extra overlap measurements than the solve itself made.
- **Map deliverable.** `summarize(k_vs_validity_map(key))` for all five systems, regenerated by
  `uv run python odmd_temple_repair_map.py`, recorded in §10 and in the BACKLOG entry.

## 6. Implementation plan (test-first)

1. Make `k_vs_validity_map` return its data (it only printed); add TOL to every comparison; add
   both LiH variants; reference via `reachable_eigenpairs` (macOS ZHEEVD-safe).
2. `tests/test_odmd_temple_repair_map_spec.py` encoding G-a/G-b/G-c (RED: `premise_refuted` missing).
3. `EnergyBracket.theta1` + `premise_refuted` (defaulted, backward-compatible).
4. Run once; record §10; any gate reality refuses is revised here with the reason.

## 7. Out of scope

- Making the substitution safe (e.g. a depth-convergence check on ε(K), Lehmann/Kato rigorous E₁
  floors) — a follow-up; this spec maps the failure, it does not fix it.
- Centered-frame ODMD (`SPEC_odmd`'s τ = π/W), shot noise, Trotter error, other geometries.
- Changing `certified_gaps` (its lower gap certificate has the same premise; not touched here).

## 8. Caveats and risks

- **R1 — brittle exact sets.** Gates use margins (> 100 mHa, > 1e-6 Ha, K ≤ 8) and structure, not
  exact cell sets; §10 records the exact sets as measured on one machine (Apple M3, macOS
  Accelerate + scipy ZHEEVR fallback). A cell within ~1e-9 Ha of a boundary can flip elsewhere.
- **R2 — the flag is one-sided.** It catches only gross overshoot (ε above the *Ritz* E₁). A small
  overshoot — the degraded cells, and every M = 2 escape — is invisible to it; detecting those
  still needs the oracle. Never read "flag silent" as "premise holds".
- Sector-restricted (reachable E₀/E₁), like everything Temple-based in this repo.

## 9. Deliverables

- `odmd_temple_repair_map.py` — `odmd_eps1`, `k_vs_validity_map`, `summarize`, `__main__` table.
- `temple_bounds.py` — `EnergyBracket.theta1`, `EnergyBracket.premise_refuted`.
- `tests/test_odmd_temple_repair_map_spec.py` — gates G-a, G-b, G-c.
- `BACKLOG.md` chem-z3h entry closed in place with the verdict; `SPEC_temple_bracket.md` G4(b)
  points here for the full-space LiH violation.

## 10. Results (gated run)

`uv run python odmd_temple_repair_map.py` (Apple M3, 2 BLAS threads, 118 s wall; the gate file
takes ~135 s). Energies in mHa; cells are (M, K); "caught" = escapes with `premise_refuted`.

| system | reach. levels / gap | self ε−E₁ @M=4 | ODMD overshoot K | escapes (M, K) | worst escape | caught | repair K @M=4 |
|---|---|---|---|---|---|---|---|
| h2 | 2 / 1620 | −3e-5 (exhausted) | — | — | — | 0/0 | n/a |
| h4 | 12 / 625 | **+14.77** | 4, 6–8, 10–14, 28 (K=5: no fit) | (2,4) (4,4) (6,4) | +0.147 (2,4) | 2/3 | 16–26 |
| lih_cas | 7 / 133 | −110.42 | 4–8, 10–18 | K = 4–7 at every M (16) | +0.589 (2,4) | 14/16 | n/a |
| lih | 27 / 133 | **+465.75** | 4–8, 10–18, 26 | M=4: K 4–8; M=6: K 4–8, 12, 14, 16; M=8: K 4–8 (18) | +0.203 (6,7) | 15/18 | 20–24, 28 |
| n2 | 18 / 708 | **+22.19** | 4–8, 10, 18–28 | K = 4, 5 at every M (8) | **+1.259** (2,4) | 6/8 | 12–16 |

- **G-a — the claim "violated on all three" is killed:** violated on H₄, N₂, full-space LiH; not
  on LiH CAS(2,5) (−110 mHa); vacuous on H₂ (2 reachable levels, exhausted by M = 2).
- **The free diagonal makes things worse (G-c).** At K = M = 4 self mode overshoots E₁ by
  +14.77 / +22.19 mHa (H₄ / N₂) yet its bracket still contains E₀ (by 0.446 / 1.640 mHa); the free
  ODMD estimate overshoots by +950.16 / +1529.71 mHa and the bracket escapes (+0.070 / +0.245 mHa).
  Repair at M = 4 first appears at K = 16 (H₄), 12 (N₂), 20 (full-space LiH): 8–16 overlaps
  beyond the 4 the solve measured.
- **Not monotone in K, as the backlog caveat predicted** — but on H₄/N₂ it degrades rather than
  breaks: H₄ is repaired at K = 16–26 and violated again at K = 28 (+0.74 mHa); N₂ repaired at
  K = 12–16 (ε 118–332 mHa *below* E₁ — valid but looser), violated again at K = 18–28
  (+14.8…+35.4 mHa); no escape at any of those cells. Full-space LiH does break at deep K (M = 6,
  K = 12–16). On LiH CAS(2,5), where self mode needed no repair, the substitution only does harm.
- **The oracle-free flag:** sound on every cell of all five systems (never fires without a real
  overshoot; never in self mode), and catches 37 of 45 escapes. The 8 misses: H₄ (2,4); N₂ (2,4),
  (2,5); LiH CAS(2,5) (2,6), (2,7); LiH (4,4), (4,5), (6,16). On H₄ and N₂ the worst escape on
  the grid is a miss — so a silent flag never certifies a bracket; full detection still needs the
  oracle (consistent with `SPEC_odmd_uq` G4: nothing single-signal sees model bias).
- **Scout numbers that did not reproduce:** the backlog's H₄ ε = −1.59181 / E₁ = −1.649658 /
  K=20 overshoot of 0.176 mHa. At the gated 0.9 Å chain E₁ = −1.555053 Ha and K = 20 sits
  3.81 mHa *below* E₁ at every M; the non-monotonicity shows up at K = 28 instead.
