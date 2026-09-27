# SPEC: NbN CAS(14,14) low-spin sector — is the committed reference in the right spin sector?

**Status:** DONE (bead chem-dc7, 2026-09-25); **Tables 1-3, G4, G5 corrected (bead chem-g1i,
2026-09-26)** — chem-dc7's numbers described a molecule with the wrong Nb-N bond length. Sector
question, on the real geometry: the CAS ground is still **S = 1**, not the committed S = 3 (still
not the (7,7) singlet either) — that qualitative call survives. Only the *sizes* of the gaps, and
the G4 hardness verdict, change: the S1-S3 gap is **1.537 mHa**, not 15.1 mHa, and the "hard
multireference benchmark" framing (G4, CONFIRMed on the wrong geometry) is **reversed** — on the
real geometry the (7,7) sector converges by D≈400-1200 like every other sector here. **chem-czw
(2026-09-27, §0b):** re-derived all four sectors on the OTHER locally-stable UHF minimum
chem-owr/chem-y6w found; decision is to KEEP the vendored chk — its orbitals give a lower
(better-converged) CASCI(14,14) energy in every sector, by 1.44-4.84 mHa.

---

## 0. Correction (chem-g1i, 2026-09-26) — read this before §3's old tables

**Root cause (found by chem-bbi, `SPEC_nbn_dmrg_reference.md`):** the CIF/chk this spec's chem-dc7
reconstruction targeted (`nbn_low_spin.py cif`, MP's reported d(Nb-N) = 2.25 Å) is not what the
*committed* `data/nbn_scf.chk` contains. The real, vendored files (`specs/nbn_scf_reference.chk`,
`specs/nbn_mp-2634.cif` — present on a persistent volume, gitignored so never committed, restored
by `nbn_dmrg_reference.ensure_reference_data()`) encode a Nb-N separation of **3.7255 Å**, not
2.2464 Å (MP's own value) — `benchmark_nbn.py`'s `ground_state_mf` builds `atom_str` from ASE's
raw `atom.position`, and the real CIF stores the N atom's fractional coordinate unwrapped (as
`(1,0,0)` rather than `(0,0,0)`, its periodic-image equivalent), so the naive Cartesian distance
is a full lattice vector too long. **Every number chem-dc7 recorded (Tables 1-3 below, unless
marked otherwise) describes this different, 2.25 Å molecule — not the one `data/nbn_scf.chk`
actually contains.** They are kept in §3 for the record, clearly marked, not deleted.

**The real numbers** (all via `nbn_dmrg_reference.load_nbn_cas()` / `run_schedule()`, which
restore the committed chk directly — sidesteps the CIF-reading bug entirely, the same route
`SPEC_nbn_dmrg_reference.md` G1-G3b use):

| S | nelec | E (Ha) | ⟨S²⟩ | Δ vs S=1 (mHa) |
|---|---|---|---|---|
| 1 | (8,6) | **−110.04756504307636** (exact FCI, unconstrained Ms=0 sector) | 2.0000095761529892 | 0 (ground) |
| 2 | (9,5) | −110.04698592997889 (exact FCI) | 6.0000045405141735 | +0.579 |
| 3 | (10,4) — **committed sector** | −110.04602841723886 (exact FCI) / −110.04602823118734 (DMRG A′) | 12.000000005607337 | +1.537 |
| 0 | (7,7) | −110.04249952166984 (DMRG headline, D≤1200, converged) | — | +5.066 |

Full provenance (wall time, host, per-D discarded weights) in §3's new Table 1 and
`results/nbn_low_spin/runs.jsonl` (2026-09-26T15:40-16:24, host `9bb17f03582d`, x86_64, 8 cores,
16 GB, except the S=1 exact-FCI row, run by chem-bbi on host `84689dcda54e`).

**G4's pre-registered hardness check also reverses.** At cheap D (≤300, the CI-gate schedules),
neither the pre-registered KILL nor CONFIRM condition fires on the real geometry (a genuine,
reproduced inconclusive result — not a threshold nudge). The headline schedules (D≤1200) resolve
it: the (7,7) sector converges (dw(1200) ≈ 3×10⁻¹⁰, two independent schedules agree to 0.3 nHa) —
**it is soft, like every other sector gated here.** The "hard multireference benchmark" framing
chem-dc7 CONFIRMed was a property of the wrong-bond-length molecule, not the real one. Full numbers
in §3's new Table 2; `tests/test_nbn_low_spin_spec.py::test_G4_*` revised accordingly (§6).

**A second, separate anomaly surfaced while correcting Table 3, DIAGNOSED 2026-09-27 (chem-owr;
consequence filed as chem-y6w, out of scope here):** re-running the reconstruction-scan pipeline
(`nbn_low_spin.py scan`, a *fresh* SCF via a synthetic single-cell CIF, independent of the
committed chk) at d = 3.7255 Å — the real molecule's actual bond length — does **not** reproduce
the committed/vendored-chk energy: E(10,4) = −110.041189 (scan) vs −110.046028 (real chk), a
4.84 mHa discrepancy at nominally the same geometry. The scan pipeline's own geometry construction
is verified bug-free (§3 Table 3: no-MIC/MIC distances agree to machine precision at every d
tried), and using the *real* vendored CIF directly (not the synthetic reconstruction) at the exact
committed geometry still reproduces the discrepancy — so it is not a geometry, basis, ECP, or MIC
bug at all. **Root cause: the committed (10,4)/S=3 UHF solution sits essentially exactly at a
`pyscf` internal-stability boundary, and whether `benchmark_nbn._tight_scf`'s stability-following
loop reports the solution stable (keeping it) or unstable (following it 0.548 mHa lower in SCF
energy, 4.84 mHa lower in downstream CAS(14,14) FCI(10,4)) depends on default multi-threaded BLAS
floating-point non-associativity, not on any geometry difference:**

- Fresh SCF on `specs/nbn_mp-2634.cif` (real CIF, no chk restore, default threads=8 on this
  container): 4 independent trials → 3 found the instability and landed on the lower solution
  (e_tot ≈ −110.02813375 Ha), 1 reported "stable" and reproduced the committed
  −110.02758538272639 Ha (matches `specs/nbn_scf_reference.chk`'s −110.0275853827266 Ha to 9
  decimals).
- The same command with `OMP_NUM_THREADS=1`: 2/2 trials deterministically and bit-identically
  (−110.02813374576796 Ha both times) find the instability and land on the lower solution.
- Both the committed solution and the lower one independently pass an explicit
  `mf.newton().stability(return_status=True)` check (`stat_i=True`) when checked directly — i.e.
  both are genuinely, separately locally-stable UHF stationary points; this is true multi-minima
  degeneracy in the SCF landscape of a near-dissociated (3.7255 Å, cf. the ~2.25 Å MP-reported
  equilibrium bond) transition-metal open-shell diatomic, not a false instability report.
- Reproducing commands: see `sandbox-handoffs/chem-owr.md`.

This means the currently-vendored chk is not a *reliably reproducible* output of the current code
— a deterministic (single-threaded) rerun of the same geometry finds a different, lower, equally
"stable" solution every time. It does not change any number already reported above (all go through
the vendored-chk path, which is a fixed, frozen artifact, not a rerun), but it means **Table 3's
distance scan is not a validated stand-in for the real geometry's absolute energies at any point**,
and it raises a question — filed as **chem-y6w** — about whether the vendored chk itself is the
reference this repo wants to keep.

**chem-y6w's decision (2026-09-27) — option (b): keep the vendored chk, caveat every dependent,
fix the regeneration path's determinism.** Re-vendoring from the deterministic lower-energy
minimum would require re-deriving every downstream number this session's predecessors produced —
`SPEC_nbn_dmrg_reference.md`'s G1-G3b *and* this spec's own Tables 1-2/S1-S3 gap/G4/G5, all of
which use the SAME UHF orbitals (loaded once from the chk, then re-used across every `nelec=`
sector via `CASCI.get_h1eff()`/`get_h2eff()` in `load_nbn_cas`) — meaning switching minima is not
a one-number fix, it changes the active-space orbitals for every sector in both specs at once.
`CLAUDE.md`'s own compute envelope note (one worker, ~8 CPUs/16 GB, "time a mid-size point before
committing to a large one") and its long-running-jobs note (hour-plus jobs should run in the
user's own terminal, not this container) both argue against attempting that full re-derivation
inside one bug-fix session — chem-bbi's and chem-g1i's own numbers alone cost ~50+ minutes of
wall-clock for the current (higher-energy) minimum; redoing them for the other minimum is the same
order of cost again. **Decision:** keep `specs/nbn_scf_reference.chk` / `data/nbn_scf.chk` as the
pinned reference, explicitly labelled metastable (every table/gate below that depends on it is
marked); separately, **fix what is actually a bug regardless of which minimum is "the" reference**
— that a fresh SCF regeneration's outcome depends on ambient thread count at all. Root cause and
fix: `mf.stability()`'s Davidson solver runs through pyscf's own OpenMP-threaded contraction
routines, and `benchmark_nbn._tight_scf` now pins `pyscf.lib.num_threads(1)` for its duration
(restoring the ambient count after) — verified bit-identical across ambient thread counts
1/2/4/8 and across 5 repeated trials at the previously-flaky ambient=8 (`specs/
SPEC_nbn_scf_determinism.md`, `tests/test_nbn_scf_determinism_spec.py` G1/G2). **This does not
make a fresh regeneration reproduce the vendored chk** — the deterministic answer is the *other*,
lower-energy minimum on the Linux container but the vendored one on macOS arm64 (BLAS-build
dependent; `SPEC_nbn_scf_determinism.md` §8, 2026-09-27); the vendored chk remains pinned as a frozen binary file,
materialized by `ensure_reference_data()`, never regenerated by the normal test/CI path. A
follow-up bead, **chem-czw**, tracks re-deriving both specs' numbers on the lower-energy minimum
if that is ever judged worth the compute; not attempted here.

## 0b. chem-czw (2026-09-27) — re-derived all four sectors on the OTHER minimum; decision: keep the vendored chk

chem-y6w deferred re-deriving Tables 1-2/G4/G5 (and `SPEC_nbn_dmrg_reference.md` G1-G3b) on the
deterministic lower-energy UHF minimum on cost grounds alone. This bead did the re-derivation
(all four spin sectors, same orbitals-shared-across-sectors structure as Table 1) to answer the
question chem-y6w left open: **should the lower-energy minimum replace the vendored chk as the
primary reference?**

**Method:** regenerated the lower-energy minimum's chk via `pyscf.lib.num_threads(1)`-pinned
`benchmark_nbn.ground_state_mf` (the fix from `SPEC_nbn_scf_determinism.md`), confirmed it lands
on the same minimum chem-owr/chem-y6w characterized (E_SCF = −110.02813374576796 Ha, reproduced to
5.7e-14 Ha on this container — cross-container float noise, not a third minimum), then ran the
*same* battery Table 1 used — exact FCI for S=1/2/3, two independent DMRG schedules for S=0 — on
these orbitals via `nbn_low_spin.py fci --chk ...` / `dmrg --chk ...` (the `fci` subcommand's
`--chk` option, previously `dmrg`-only, was added for this — mirrors the existing pattern 1:1).

**Table 4 — spin ladder, LOWER-ENERGY minimum's orbitals (chem-czw, 2026-09-27; host
`e1fd4a2d0570`, 8-core x86_64, 16 GB, no GPU). Cross-check against Table 1, not a replacement.**

| S | nelec | E (Ha) | ⟨S²⟩ | Δ vs S=1 (mHa) | wall | Δ vs Table 1 (same sector, mHa) |
|---|---|---|---|---|---|---|
| 1 | (8,6), exact FCI | −110.04577676472408 | 2.0000222708995032 | 0 (ground) | 259.7 s / 8 cores | **+1.788** (Table 1 lower) |
| 2 | (9,5), exact FCI | −110.04546990757419 | 6.000… | +0.307 | 192.2 s / 8 cores | **+1.516** (Table 1 lower) |
| 3 | (10,4), exact FCI — committed sector's nelec | −110.04119080464493 | 12.000000000270804 | +4.586 | 5.6 s / 8 cores | **+4.838** (Table 1 lower; matches chem-owr's originally-cited "4.84 mHa" almost exactly) |
| 0 | (7,7), DMRG A perD 400/800/1200 | −110.0410570692685 | — | +4.720 | 901.6 s / 4 threads | **+1.440** (Table 1 lower) |
| 0 | (7,7), DMRG B ramp 300/600/1200 | −110.04105709905623 | — | +4.720 | 444.0 s / 4 threads (\|A−B\|=2.98e-8 Ha, agrees) | — |

Raw records: `results/nbn_low_spin/runs.jsonl` (rows with `"chk": "results/nbn_lower_minimum/
nbn_scf_lower.chk"`), regenerable via `scripts/gen_nbn_lower_minimum_chk.py`.

**Findings:**

1. **Ordering is orbital-choice-independent.** S=1 < S=2 < S=3 < S=0 on the lower minimum's
   orbitals too — the chem-dc7/chem-g1i "committed sector is not the CAS ground" call, and the
   softness of the (7,7) sector (two independent DMRG schedules agree to 2.98e-8 Ha, same order as
   Table 1's 2.9e-10 Ha), both survive on the other minimum. Nothing here undermines the repo's
   existing sector-ordering or hardness conclusions.
2. **The vendored/committed chk's orbitals give a strictly LOWER (more favorable) CASCI(14,14)
   energy than the lower-SCF minimum's orbitals, in every sector checked — by 1.44 to 4.84 mHa.**
   This is the opposite direction from the raw whole-molecule SCF energy (the vendored chk is
   0.548 mHa *higher* there). "Lower whole-molecule HF energy" and "better CAS(14,14) active-space
   orbitals for this multireference problem" are not the same claim — CASCI is not variational
   over orbital choice (only over the CI expansion within a fixed orbital set); comparing CASCI
   totals across two different UHF references is a standard but informal empirical quality
   criterion, not a rigorous bound (a rigorous comparison would need CASSCF orbital relaxation,
   out of scope per §7). Taken as such, it is a real, reproducible, cite-worthy finding: the
   vendored chk's orbitals are empirically the better CAS(14,14) reference of the two, in addition
   to chem-y6w's compute-cost argument for not switching.

**Decision: keep `specs/nbn_scf_reference.chk` / `data/nbn_scf.chk` as the primary reference.**
This upgrades chem-y6w's deferral (cost alone) to a substantive reason (finding 2 above) — the
committed chk is not merely cheaper to keep, it demonstrably gives a better-converged CAS(14,14)
description in every sector tested. Table 1, `SPEC_nbn_dmrg_reference.md` G1-G3b, and every other
number anchored to the vendored chk are correct as recorded and are NOT overwritten. This decision
is pinned as a regression gate: `tests/test_nbn_czw_lower_minimum_spec.py` (own process; ~12s,
plus a fixture-gated ~10-20s live regeneration sub-test) — fails if a future change makes the
lower-SCF minimum's orbitals equal or beat the vendored chk's CASCI energy in any of the four
sectors, or if the sector ordering stops holding on the other minimum.

**What was not done:** DMRG on the lower minimum's (10,4)/(9,5)/(8,6) sectors (redundant — exact
FCI already covers them, same as Table 1's approach); CASSCF orbital relaxation on either minimum
(would settle the "better reference" question rigorously but is out of scope, §7); re-deriving
`SPEC_nbn_dmrg_reference.md`'s own G1-G3b gates on the lower minimum (its (10,4) exact-FCI number
above, −110.04119080464493 Ha, is the number those gates would need; not threaded back into that
spec's own gate file since the decision is to keep the vendored chk, not switch).

## 1. Goal

`SPEC_nbn_dmrg_reference.md` gated only the high-spin nelec=(10,4) sector (block2 SU(2), S = 3)
and recorded an ungated low-spin E(7,7) "3.5 mHa above". Claim under test (specs/BACKLOG.md,
"Is NbN's flagged 'hard multireference benchmark' actually hard?"): the (7,7) sector is genuinely
hard, and the sector ordering survives at converged bond dimension.

**Pre-registered criteria (fixed before any number was seen; not revised):**

- **KILL** the hard-benchmark framing if dw(D=300) < 1e-7 **and** per-D spread < 1e-5 Ha.
- **CONFIRM** it if dw > 1e-5 **or** |E_A′ − E_B′| > 0.1 mHa.
- Report the converged sector ordering; if the committed sector is not the ground, correct the
  committed reference, don't just caveat it.
- The gate skips cleanly on a fresh clone.

Applied to the real geometry (§0), neither KILL nor CONFIRM fires at the pre-registered D=300
checkpoint — an honest gap in the pre-registered design, not anticipated when the criteria were
written, closed by going to converged D (§3 Table 2, §6).

## 2. Honest framing — read this before the numbers

- **The original geometry is lost — corrected: it was not.** (chem-bbi, `SPEC_nbn_dmrg_reference.md`,
  2026-09-26) `data/nbn_scf.chk` and its source CIF were gitignored, not deleted; they persisted on
  this container's `data/` volume and are now vendored at `specs/nbn_scf_reference.chk` /
  `specs/nbn_mp-2634.cif`. `nbn_dmrg_reference.ensure_reference_data()` restores them into `data/`
  on a fresh clone.
- **The 2.25 Å reconstruction does NOT reproduce the committed number — now explained, not
  mysterious.** At d = 2.25 Å the reconstruction pipeline gives E(10,4) = −110.056110 (exact FCI),
  not −110.046028 — because the *real* committed molecule sits at 3.7255 Å, not 2.25 Å (§0). Every
  number in the original §3 tables below is for the 2.25 Å reconstruction, a real but different
  diatomic; it is not the committed one and none of it should be read as describing
  `data/nbn_scf.chk`.
- **To make the sector conclusion independent of the lost — corrected: mis-targeted — geometry,**
  the ordering was checked at nine distances 2.0–3.0 Å (§3, table 3) via the *same* synthetic,
  self-consistent reconstruction pipeline (not the committed chk). It holds at every one, and
  (chem-g1i) also at the real d = 3.7255 Å via that pipeline — but see §0's second anomaly: that
  pipeline's own absolute energies do not agree with the committed chk's, at any distance tried.
  The ordering-robustness claim (S ≥ 2 beats S = 3 across a family of hypothetical bond lengths)
  stands on its own terms; it is not evidence about the real geometry's absolute energies, which
  come only from Table 1 (real, vendored-chk numbers).
- **"(7,7)" means two different things here.** block2 runs `SymmetryTypes.SU2` with
  spin = na − nb, so DMRG nelec=(7,7) is the **S = 0 singlet**. PySCF FCI in nelec=(7,7) is the
  **Ms = 0 determinant sector**, which contains an Ms = 0 component of *every* S (the CASCI
  integrals are spin-restricted, so H is spin-free). Its unconstrained lowest root is therefore
  the global CAS ground of any spin — the single calculation that answers "right sector?".

## 3. Results

Raw record: `results/nbn_low_spin/runs.jsonl` (one JSON line per run, per-D weights included).

### Table 1 — spin ladder, REAL geometry (chem-g1i, 2026-09-26; host `9bb17f03582d`: 8-core
x86_64, 16 GB, no GPU, except the marked exact-FCI S=1 row, run by chem-bbi on host
`84689dcda54e`, also 8 cores). All rows via `nbn_dmrg_reference.load_nbn_cas`/`run_schedule`,
which restore the committed `data/nbn_scf.chk` (or its vendored copy) directly.

| S | how | E (Ha) | ⟨S²⟩ | Δ vs S=1 (mHa) | wall |
|---|---|---|---|---|---|
| 1 | FCI, Ms=0 sector, unconstrained (1.18e7 dets) | **−110.04756504307636** | 2.0000095761529892 | 0 | 1010.6 s / 8 cores |
| 1 | DMRG A′ nelec=(8,6), cheap (D≤300) | −110.04756608985913 | — | −0.00105 (overshoot) | 144.0 s / 2 threads (chem-bbi) |
| 2 | FCI, (9,5) sector (4.0e6 dets) | −110.04698592997889 | 6.0000045405141735 | +0.579 | 205.4 s / 8 cores |
| 3 | FCI, (10,4) sector (1.0e6 dets) — **the committed sector** | −110.04602841723886 | 12.000000005607337 | +1.537 | 13.2 s / 8 cores |
| 3 | DMRG A′ nelec=(10,4), cheap — matches historical −110.046028 to 6 decimals | −110.04602823118734 | — | +1.537 | ~205 s / 2 threads (chem-bbi) |
| 0 | DMRG A perD 400/800/1200, nelec=(7,7) | −110.04249952166984 | — | +5.066 | 1029.4 s / 4 threads |
| 0 | DMRG B ramp 300/600/1200, nelec=(7,7) | −110.04249952196103 | — | +5.066 | 530.0 s / 4 threads (\|A−B\|=2.9e-10 Ha) |
| 0 | DMRG A′ perD 100/200/300 (CI-gate schedule) | −110.04250013542114 | — | +5.066 | 238.1 s / 2 threads |
| 0 | DMRG B′ ramp 80/160/300 (CI-gate schedule) | −110.04249595430977 | — | +5.065 | 124.3 s / 2 threads |

Ordering S=1 < S=2 < S=3 < S=0 (S=1 the CAS ground) is unchanged from the wrong-geometry table;
only the gaps shrank by roughly an order of magnitude (15.1→1.537 mHa for S=3, 39.3→5.066 mHa for
S=0).

### Table 2 — the pre-registered test, (7,7) singlet, REAL geometry (chem-g1i, 2026-09-26).

**Cheap schedules (D ≤ 300, the CI-gate dims):**

| schedule | D ladder | dw at D=300 | per-D spread (100/80→300) | E (Ha) |
|---|---|---|---|---|
| A′ perD | 100/200/300 | **1.0533e-7** | 3.14e-5 Ha | −110.04250013542114 |
| B′ ramp | 80/160/300 | **1.0362e-7** | 6.07e-5 Ha | −110.04249595430977 |

\|E_A′ − E_B′\| = 4.181e-6 Ha (4.18 µHa) — far under CONFIRM's 1e-4 Ha floor. dw(300) sits ~5%
*outside* the 1e-7 KILL floor on both schedules; per-D spread sits 3-6× outside the 1e-5 Ha KILL
floor on the high side. **Neither pre-registered condition fires** — genuinely inconclusive at
D≤300, not a forced or nudged call.

**Headline schedules (D ≤ 1200) — the definitive answer:**

| schedule | D ladder | dw at D=1200 | E (Ha) | wall |
|---|---|---|---|---|
| A perD | 400/800/1200 | **3.061e-10** | −110.04249952166984 | 1029.4 s / 4 threads |
| B ramp | 300/600/1200 | **2.733e-10** | −110.04249952196103 | 530.0 s / 4 threads |

\|E_A − E_B\| = 2.91e-10 Ha (0.3 nHa) — five orders of magnitude tighter than the pre-registered
CONFIRM floor. dw already drops to 3.86e-8 by D=400 alone (an order of magnitude under the KILL
floor). **Conclusion: on the real geometry the (7,7) sector is soft, like the committed (10,4)
sector — the "hard multireference benchmark" framing chem-dc7 CONFIRMed does not survive.** The
cheap-D check landed inconclusive here only because D=300 sits close to this sector's particular
convergence knee on this molecule, not because the sector is actually hard.

**Table 2b (chem-dc7's own-natural-orbitals cross-check) is NOT redone on the real geometry** —
left as a gap, called out in §7. The wrong-geometry numbers below (§3-old) do not carry over; no
claim is made about them here either way.

### Table 3 — reconstruction-independence scan (self-consistent synthetic family, NOT the
committed geometry — see §0's second anomaly).

Unchanged from chem-dc7 (2.0-3.0 Å, table below, kept for the record) **plus** one new point at
the real molecule's actual bond length, added by chem-g1i:

| d (Å) | SCF 2S | E(10,4) S≥3 | E(9,5) S≥2 | E(S≥2) − E(S=3) (mHa) | lower of the two |
|---|---|---|---|---|---|
| 3.7255 (chem-g1i, = the real Nb-N separation) | 6 | −110.041189 | −110.045469 | −4.28 | S=2 |

Ordering (S≥2 below S=3) holds here too — consistent with Table 1's real-geometry finding that
S=1/S=2 both beat S=3. **But the absolute energies do not match Table 1's real-chk numbers at this
same distance** (E(10,4): −110.041189 here vs −110.046028 from the real chk — 4.84 mHa apart).
The synthetic scan's own geometry construction is verified bug-free (a direct
`write_cif`/`get_distance` check shows no-MIC and MIC distances agree to machine precision at
every d tried, including 2.25 and 3.7255 Å — this is a *different* code path issue than the
CIF-unwrapping bug that misidentified the committed geometry in the first place). **Diagnosed
2026-09-27 (chem-owr, §0):** not a geometry/basis/ECP bug — the (10,4)/S=3 UHF solution has (at
least) two genuinely locally-stable minima 0.548 mHa apart in SCF energy (4.84 mHa apart after
CASCI(14,14)/FCI), and default multi-threaded BLAS non-associativity makes which one
`benchmark_nbn._tight_scf`'s stability-following loop lands on non-deterministic run to run; a
deterministic single-threaded rerun reliably finds the lower one, not the committed one. Use
Table 1, not this table, for any claim about the real molecule's absolute energies; use this table
only for the ordering-robustness claim over a family of hypothetical geometries. Whether the
vendored chk should be replaced by the lower, single-thread-reproducible solution is filed as
**chem-y6w**, not resolved here.

<details>
<summary>Original (chem-dc7, 2026-09-25) Table 1-3 — WRONG-BOND-LENGTH RECONSTRUCTION, d(Nb-N)
nominally "2.25 Å" but see §0: the real separation this same reconstruction pipeline implies at
the committed molecule's own SCF-preferred spin is 3.7255 Å, not 2.25 Å. Kept for the record, not
as a reference for the committed molecule.</summary>

**Table 1 (old) — spin ladder at d = 2.25 Å, exact FCI (the reference) vs DMRG.**

| S | how | E (Ha) | ⟨S²⟩ | Δ vs S=1 (mHa) | wall |
|---|---|---|---|---|---|
| 1 | FCI, Ms=0 sector, unconstrained (1.18e7 dets) | −110.0712455 | 2.000 | 0 | 19 min |
| 1 | DMRG A′ nelec=(8,6), dweight-extrapolated | −110.0712315 ± 0.045 mHa | — | +0.014 | 6 min |
| 1 | DMRG B′ nelec=(8,6) | −110.0715881 ± 0.29 mHa | — | −0.34 | 2 min |
| 2 | FCI, (9,5) sector, unconstrained (4.0e6 dets) | −110.0695221 | 6.000 | +1.72 | 3 min |
| 3 | FCI, (10,4) sector (1.0e6 dets) — the committed sector (name only — wrong molecule) | −110.0561097 | 12.000 | +15.14 | 1 min |
| 0 | DMRG A perD 400/800/1200, nelec=(7,7) | −110.0319893 ± 0.002 mHa | — | +39.26 | 30 min |
| 0 | DMRG B ramp 300/600/1200, nelec=(7,7) | −110.0319942 ± 0.008 mHa | — | +39.25 | 33 min |
| 0 | FCI, (7,7) with S²-penalty pinned to S=0 | did not converge (2 attempts, 300 Davidson cycles each; best upper bound −110.031418) | — | — | 82 + 27 min |

**Table 2 (old) — the pre-registered test, (7,7) singlet, cheap schedules (D ≤ 300).**

| schedule | D ladder | dw at D=300 | per-D spread | E (Ha) |
|---|---|---|---|---|
| A′ perD | 100/200/300 | 4.8e-5 | 2.26 mHa | −110.0321533 |
| B′ ramp | 80/160/300 | 4.8e-5 | 3.87 mHa | −110.0320896 |

\|E_A′ − E_B′\| = 0.064 mHa (< 0.1 mHa, so that arm does not fire); dw(300) = 4.8e-5 > 1e-5 ⇒
CONFIRM (on this molecule). The KILL condition misses by ~500× in dw and ~200× in spread.

**Table 2b (old) — cross-check in the singlet's OWN natural orbitals** (block2 SU(2) D=500 1-RDM,
NOs sorted by occupation). Occupations: 1.997, 1.996, 1.996, 1.959, 1.592, 1.455, 1.028, 1.025,
0.503, 0.263, 0.140, 0.026, 0.012, 0.008.

| schedule | basis | dw at D=300 | E(D=300) | E extrapolated |
|---|---|---|---|---|
| A′ | SCF (septet UHF) | 4.8e-5 | −110.0318512 | −110.0321533 |
| A′ | singlet NOs | 7.2e-5 | −110.0318387 | −110.0322844 |
| B′ | singlet NOs | 7.1e-5 | −110.0318302 | −110.0324142 |

**Table 3 (old) — lowest S≥2 (FCI (9,5)) vs S=3 (FCI (10,4)) across d.**

| d (Å) | SCF 2S | E(10,4) S≥3 | E(9,5) S≥2 | E(S≥2) − E(S=3) (mHa) | lower of the two |
|---|---|---|---|---|---|
| 2.00 | 6 | −110.048525 | −110.087822 | -39.3 | S=2 |
| 2.10 | 6 | −110.056333 | −110.082257 | -25.9 | S=2 |
| 2.20 | 6 | −110.057331 | −110.074056 | -16.7 | S=2 |
| 2.30 | 6 | −110.030642 | −110.041952 | -11.3 | S=2 |
| 2.40 | 6 | −110.023289 | −110.029678 | -6.4 | S=2 |
| 2.50 | 6 | −110.017942 | −110.016381 | +1.6 | S=3 (but S=1 is 28 mHa lower, DMRG) |
| 2.60 | 6 | −110.022624 | −110.025552 | -2.9 | S=2 |
| 2.80 | 6 | −110.027910 | −110.041455 | -13.5 | S=2 |
| 3.00 | 6 | −110.040977 | −110.042432 | -1.5 | S=2 |
| 2.25 (main) | 6 | −110.056110 | −110.069522 | −13.4 | S=2 (S=1 lower still, FCI) |

</details>

## 4. What this means

1. **The committed reference names the wrong sector.** CASCI(14,14) prefers S = 1 over the
   committed S = 3 at every geometry checked, real or reconstructed. On the real, vendored
   geometry the gap is **1.537 mHa** (not the 15.1 mHa recorded against the wrong-bond-length
   reconstruction). `SPEC_nbn_dmrg_reference.md` is corrected (chem-bbi) to carry this number.
2. **The singlet gap also shrinks, but the qualitative call is unchanged: not a close competitor.**
   Real geometry: S=0 is 5.066 mHa above S=1 (was 39.3 mHa on the wrong molecule) — still well
   above S=2/S=3, not a near-degenerate state.
3. **The hard-benchmark framing does NOT survive on the real geometry (REVERSED from chem-dc7).**
   The (7,7) sector converges by D≈400-1200 (dw → 3×10⁻¹⁰), just like the committed (10,4)
   sector. chem-dc7's CONFIRM verdict was a property of the wrong-bond-length molecule; NbN
   CAS(14,14) should be retired as a strong-correlation target in every sector checked so far, not
   just the high-spin one `SPEC_nbn_dmrg_reference.md` originally gated.
4. **"FCI-intractable" is still not true of this CAS on this hardware.** PySCF's direct FCI solved
   every fixed-Sz sector tried (up to 1.18e7 determinants) in well under 20 minutes on 8 cores.

## 5. Public interface

```
nbn_dmrg_reference.load_nbn_cas(..., nelec=None)      # nelec=(na, nb) forces an active sector
nbn_dmrg_reference.run_schedule(name, n_threads, nelec=None)
nbn_low_spin.py cif                                     # reconstructed CIF + regenerated chk (wrong bond length; see §0)
nbn_low_spin.py fci  --nelec NA NB [--twos 2S] [--save-no]
nbn_low_spin.py dmrg --schedule S --nelec NA NB [--orbitals scf|no] [--chk PATH]
nbn_low_spin.py scan --d D1 D2 ...                      # (10,4) vs (9,5) exact FCI per distance (synthetic reconstruction, NOT the committed chk — §0)
```

## 6. Acceptance gates — `tests/test_nbn_low_spin_spec.py` (own process; ~10-11 min on 8 cores)

- **G4 — pre-registered hardness, real geometry (chem-g1i, revised).** Asserts the actual,
  reproduced cheap-D (≤300) outcome: neither the pre-registered KILL nor CONFIRM condition fires,
  and pins the near-boundary dw/agreement numbers as a regression check (order-of-magnitude bands
  around the measured values — not re-derived "safe" thresholds). The *definitive* verdict (soft,
  reversing chem-dc7's CONFIRM) comes from the headline D≤1200 runs, recorded in §3 Table 2 but
  not re-run every CI cycle (~9-17 min each) — same convention as the (10,4) sector's own headline
  record in `SPEC_nbn_dmrg_reference.md`. **Metastable-chk caveat (chem-y6w, §0):** every number
  this gate pins comes from orbitals anchored to the vendored, metastable (10,4)/S=3 UHF solution.
- **G5 — the committed sector is not the CAS ground (chem-g1i, recalibrated).** DMRG A′ S=1 below
  A′ S=3 by > 1 mHa (measured: 1.538 mHa — margin chosen with ~35% headroom below the measured
  gap and ~1000× above the ~1 µHa DMRG/FCI cross-check noise floor), and the singlet above the
  septet. Green. **Metastable-chk caveat (chem-y6w, §0):** the 1.538 mHa gap itself is only
  measured for the vendored solution; chem-owr found the *other* locally-stable minimum shifts the
  (10,4) sector's own CASCI energy by 4.84 mHa (3× the gap size) — this gate does not (and, per
  chem-y6w's decision, does not currently plan to) check whether S=1 still beats S=3 on that other
  minimum.
- Both skip (not fail) without `data/nbn_scf.chk` (now restored on a fresh clone by
  `ensure_reference_data()`, wired into `load_nbn_cas`) or block2.
- **G6 — chem-czw, the vendored chk beats the other minimum's orbitals (own file,
  `tests/test_nbn_czw_lower_minimum_spec.py`, ~12s + one ~10-20s fixture-gated live regen).**
  Pins Table 4 above: the sector ordering on the lower-SCF-energy minimum's orbitals, the two
  independent DMRG schedules' agreement, and — the decision-bearing check — that the vendored
  chk's CASCI(14,14) energy is strictly lower than the other minimum's in all four sectors.

**Superseded, not fixed:** `test_nbn_dmrg_reference_spec.py::G3` (dw(300) < 1e-7 for the (10,4)
sector) was recorded failing against chem-dc7's wrong-bond-length reconstruction; it passes again
on the real geometry (chem-bbi) — see that spec.

## 7. Out of scope

- Recovering why the real, vendored CIF's N-atom fractional coordinate is stored unwrapped in the
  first place (a materialsproject.org export convention, ASE's reader, or something else) —
  `SPEC_nbn_dmrg_reference.md` records the bug's effect; its origin is not diagnosed.
- **chem-owr (filed 2026-09-26, diagnosed 2026-09-27):** the `nbn_low_spin.py scan` reconstruction
  pipeline does not reproduce the committed/vendored-chk energy even at the matching real bond
  length (3.7255 Å) — a 4.84 mHa discrepancy. Root cause found: SCF-solution multiplicity (two
  genuinely locally-stable UHF minima 0.548 mHa apart) plus BLAS-threading non-determinism in
  pyscf's stability-following loop, not a geometry/basis/ECP bug (§0). **chem-y6w (decided
  2026-09-27):** keep the vendored chk as reference (re-deriving Tables 1-2/G4/G5 and
  `SPEC_nbn_dmrg_reference.md` G1-G3b on the other minimum is out of scope, filed separately); the
  regeneration path's non-determinism itself is fixed (`specs/SPEC_nbn_scf_determinism.md`) — a
  fresh regeneration is now reproducible, just not equal to the vendored chk.
- **Re-deriving this spec's Tables 1-2/S1-S3 gap/G4/G5 on the deterministic lower-energy UHF
  minimum instead of the vendored one — DONE (chem-czw, 2026-09-27, §0b/Table 4).** All four
  sectors re-derived; decision is to keep the vendored chk (§0b). Still out of scope, and not
  attempted: re-deriving `SPEC_nbn_dmrg_reference.md`'s own G1-G3b gate *file* on the other
  minimum (chem-czw's Table 4 has the (10,4) number that spec's gates would need, but the gates
  themselves are not switched over, since the decision is not to switch references) and CASSCF
  orbital relaxation to settle the "better reference" question rigorously (see §0b).
- Redoing Table 2b (the singlet's own-natural-orbitals cross-check) on the real geometry — flagged
  as a gap in §3, not filled.
- Remeasuring Table 3's full nine-distance scan's absolute energies against the real chk (only the
  qualitative ordering claim is reused; the scan pipeline's own energies are shown, chem-owr
  notwithstanding, not to be equated with Table 1's).
- CASSCF orbital relaxation per spin state (unchanged from chem-dc7: CASCI orbitals come from the
  UHF(S=3) solution throughout).
- Materials claims of any kind (finite diatomic, LANL2DZ ECP, 6-31G* N).

## 8. Caveats

- CASCI energies jump by up to ~25 mHa between neighbouring distances in the *synthetic scan*
  pipeline (chem-dc7): which UHF orbitals land in the 14-orbital window changes with d. Table 1's
  real-geometry numbers do not have this issue (fixed geometry, one SCF, restored not regenerated).
- DMRG cheap-D dweight extrapolation overshoots the exact FCI reference by ~1 µHa for S=1 and (on
  the wrong-geometry reconstruction) by up to 0.34 mHa for a different cheap schedule — small next
  to every margin used in G4/G5, but a reminder that cheap dweight numbers are extrapolated, not
  raw sweep energies (see `SPEC_nbn_dmrg_reference.md` R2).
- `data/` being gitignored means any developer's own pre-existing checkpoint under that name takes
  precedence over the vendored copy (`ensure_reference_data()` never overwrites); "passes here" is
  not yet confirmed on a literal fresh `git clone` end-to-end (chem-bbi's R3).
- **The vendored (10,4)/S=3 chk itself is metastable, not uniquely reproducible (chem-owr,
  chem-y6w).** `benchmark_nbn._tight_scf`'s stability-following Newton/SOSCF loop lands on one of
  (at least) two genuinely locally-stable UHF solutions 0.548 mHa apart in SCF energy (4.84 mHa
  apart after CASCI(14,14)/FCI in the (10,4) sector) depending on multi-threaded-BLAS floating-point
  non-associativity in pyscf's Davidson-based stability check; a deterministic
  (`OMP_NUM_THREADS=1`) rerun reproducibly finds the *lower* one, not the one currently vendored at
  `specs/nbn_scf_reference.chk`. Every number in Tables 1-2 above, and `SPEC_nbn_dmrg_reference.md`
  G1-G3b, is anchored to the higher (committed) solution; none of them are wrong given that chk,
  but the chk's own status as "the" SCF ground state for this sector is now in question.
