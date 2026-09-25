# SPEC: NbN CAS(14,14) low-spin sector — is the committed reference in the right spin sector?

**Status:** DONE (bead chem-dc7, 2026-09-25). Pre-registered outcome: **CONFIRM** (the (7,7)
singlet is not soft). Sector question: the committed reference names the **wrong sector** — the
CAS ground is **S = 1**, not the committed S = 3 and not the (7,7) singlet. G4 green, G5 green;
`SPEC_nbn_dmrg_reference.md` G3 **fails** on the regenerated checkpoint (recorded, not loosened).

---

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

## 2. Honest framing — read this before the numbers

- **The original geometry is lost.** `data/nb_structures/NbN_mp-2634.cif` and `data/nbn_scf.chk`
  are gitignored and were never committed; materialsproject.org is unreachable from the session
  that ran this. mp-2634 is WC-type NbN (P-6m2, Nb (0,0,0), N (1/3,2/3,1/2)), so the 2-atom
  "cluster" is a **diatomic** whose only parameter is d(Nb–N); MP's own description gives
  "All Nb–N bond lengths are 2.25 Å". `nbn_low_spin.py cif` writes that reconstruction and
  regenerates the checkpoint through the committed `benchmark_nbn.ground_state_mf` path.
- **The reconstruction does NOT reproduce the committed number.** At d = 2.25 Å the same pipeline
  gives E(10,4) = −110.056110 (exact FCI), not −110.046028 — 10.1 mHa lower. A fine scan over the
  whole rounding window 2.240–2.255 Å stays in −110.05594…−110.05643 (`kind: fingerprint` rows). So either the original CIF
  was not at MP's reported d, or something else (pyscf version, SCF solution) differed. Every
  number below is for the reconstruction; the committed −110.046028 is **not reproducible** from
  anything in the repo.
- **To make the sector conclusion independent of the lost geometry,** the ordering is checked at
  nine distances 2.0–3.0 Å (§3, table 3). It holds at every one.
- **"(7,7)" means two different things here.** block2 runs `SymmetryTypes.SU2` with
  spin = na − nb, so DMRG nelec=(7,7) is the **S = 0 singlet**. PySCF FCI in nelec=(7,7) is the
  **Ms = 0 determinant sector**, which contains an Ms = 0 component of *every* S (the CASCI
  integrals are spin-restricted, so H is spin-free). Its unconstrained lowest root is therefore
  the global CAS ground of any spin — the single calculation that answers "right sector?".

## 3. Results (host `vm`: 4-core x86_64, 15 GB, no GPU; wall times in `runs.jsonl`)

Raw record: `results/nbn_low_spin/runs.jsonl` (one JSON line per run, per-D weights included).

**Table 1 — spin ladder at d = 2.25 Å, exact FCI (the reference) vs DMRG.**

| S | how | E (Ha) | ⟨S²⟩ | Δ vs S=1 (mHa) | wall |
|---|---|---|---|---|---|
| 1 | FCI, Ms=0 sector, unconstrained (1.18e7 dets) | **−110.0712455** | 2.000 | 0 | 19 min |
| 1 | DMRG A′ nelec=(8,6), dweight-extrapolated | −110.0712315 ± 0.045 mHa | — | +0.014 | 6 min |
| 1 | DMRG B′ nelec=(8,6) | −110.0715881 ± 0.29 mHa | — | −0.34 | 2 min |
| 2 | FCI, (9,5) sector, unconstrained (4.0e6 dets) | −110.0695221 | 6.000 | +1.72 | 3 min |
| 3 | FCI, (10,4) sector (1.0e6 dets) — **the committed sector** | −110.0561097 | 12.000 | +15.14 | 1 min |
| 0 | DMRG A perD 400/800/1200, nelec=(7,7) | −110.0319893 ± 0.002 mHa | — | +39.26 | 30 min |
| 0 | DMRG B ramp 300/600/1200, nelec=(7,7) | −110.0319942 ± 0.008 mHa | — | +39.25 | 33 min |
| 0 | FCI, (7,7) with S²-penalty pinned to S=0 | **did not converge** (2 attempts, 300 Davidson cycles each; best upper bound −110.031418) | — | — | 82 + 27 min |

**Table 2 — the pre-registered test, (7,7) singlet, cheap schedules (D ≤ 300).**

| schedule | D ladder | dw at D=300 | per-D spread | E (Ha) |
|---|---|---|---|---|
| A′ perD | 100/200/300 | **4.8e-5** | 2.26 mHa | −110.0321533 |
| B′ ramp | 80/160/300 | **4.8e-5** | 3.87 mHa | −110.0320896 |

|E_A′ − E_B′| = 0.064 mHa (< 0.1 mHa, so that arm does not fire); dw(300) = 4.8e-5 > 1e-5 ⇒
**CONFIRM**. The KILL condition misses by ~500× in dw and ~200× in spread. Compare the septet on
the same checkpoint: dw(300) = 3.0e-7 and per-D spread 0.088 mHa (A′ −110.056118, 8 µHa from exact FCI) — the
singlet is ~160× less converged at equal D. At converged D the singlet's cheap dweight
extrapolation **overshoots** by ~0.16 mHa (A′ −110.03215 vs A/B −110.03199): D ≤ 300 is not in the
asymptotic regime for this sector.

**Table 2b — cross-check in the singlet's OWN natural orbitals** (block2 SU(2) D=500 1-RDM, NOs
sorted by occupation; `nbn_low_spin.py no`, then `dmrg --orbitals no`). Occupations:
1.997, 1.996, 1.996, 1.959, **1.592, 1.455, 1.028, 1.025, 0.503, 0.263, 0.140**, 0.026, 0.012,
0.008 — two near-singly-occupied orbitals and five more far from 0/2: an open-shell,
genuinely multireference singlet, not an artefact of borrowing septet orbitals.

| schedule | basis | dw at D=300 | E(D=300) | E extrapolated |
|---|---|---|---|---|
| A′ | SCF (septet UHF) | 4.8e-5 | −110.0318512 | −110.0321533 |
| A′ | singlet NOs | **7.2e-5** | −110.0318387 | −110.0322844 |
| B′ | singlet NOs | **7.1e-5** | −110.0318302 | −110.0324142 |

Energy invariance under the in-CAS rotation holds to 12 µHa at D=300 (truncation-level). The
discarded weight does *not* drop in the sector's own NOs (occupation ordering is not an
entanglement-optimised site ordering — not tried here), and |E_A′ − E_B′| = 0.13 mHa now fires the
second CONFIRM arm too. The verdict is not an orbital-choice artefact.

**Table 3 — reconstruction independence: lowest S≥2 (FCI (9,5)) vs S=3 (FCI (10,4)) across d.**

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

At d = 2.5 Å, the one point where the septet beats the quintet, DMRG A′ for S = 1 gives
−110.04594 (D=300: −110.04582), still **28 mHa below** the septet. The septet is the CAS ground at
**none** of the nine distances. (The SCF spin scan picks 2S = 6 at every d ≥ 2.0 Å — the SCF-level
spin ranking and the CAS-level ranking disagree; that disagreement is the finding.)

## 4. What this means

1. **The committed reference names the wrong sector.** The spin-scanned UHF prefers S = 3, and
   `load_nbn_cas` inherits that as nelec=(10,4); CASCI(14,14) prefers S = 1 by 15.1 mHa at
   d = 2.25 Å (and S = 1 or S = 2 below S = 3 at every distance scanned). The quantity
   `SPEC_nbn_dmrg_reference.md` certified is the energy of an **excited spin state**. That spec
   is corrected (not caveated) to name S = 1 as the ground and carry the exact FCI number.
2. **The "3.5 mHa singlet gap" does not reproduce.** On the reconstruction the singlet is
   24.1 mHa *above* the septet and 39.3 mHa above the true ground. The singlet is not a
   near-degenerate competitor; the triplet and quintet are.
3. **The hard-benchmark framing survives** (CONFIRM): the singlet needs D ≈ 1200 for
   dw < 1e-6. But the right follow-up target is the **S = 1 ground** (dw(300) = 1.8e-5, also
   above 1e-5), not the singlet.
4. **"FCI-intractable" is not true of this CAS on this hardware.** PySCF's direct FCI solved the
   1.18e7-determinant Ms = 0 sector in 19 min on 4 cores / < 3 GB. The repo's 5e6 cutoff is
   operational (as `SPEC_nbn_dmrg_reference.md` §8 says), and it hid an exact answer that
   overturns the reference.

## 5. Public interface

```
nbn_dmrg_reference.load_nbn_cas(..., nelec=None)      # nelec=(na, nb) forces an active sector
nbn_dmrg_reference.run_schedule(name, n_threads, nelec=None)
nbn_low_spin.py cif                                     # reconstructed CIF + regenerated chk
nbn_low_spin.py fci  --nelec NA NB [--twos 2S] [--save-no]
nbn_low_spin.py dmrg --schedule S --nelec NA NB [--orbitals scf|no] [--chk PATH]
nbn_low_spin.py scan --d D1 D2 ...                      # (10,4) vs (9,5) exact FCI per distance
```

## 6. Acceptance gates — `tests/test_nbn_low_spin_spec.py` (own process; ~8 min on 4 cores)

- **G4 — pre-registered hardness (DEFINITION OF DONE).** (7,7) A′/B′: the KILL condition does not
  hold and the CONFIRM condition does. Green.
- **G5 — the committed sector is not the CAS ground.** DMRG A′ S=1 below A′ S=3 by > 5 mHa
  (exact: 15.1), and the singlet above the septet. Green.
- Both skip (not fail) without `data/nbn_scf.chk` or block2; so does
  `tests/test_nbn_dmrg_reference_spec.py` now (G1 needs only the chk).

**Recorded failure, not fixed:** `test_nbn_dmrg_reference_spec.py::G3` asserts septet
dw(300) < 1e-7; on the regenerated checkpoint it is 3.0e-7. The threshold was calibrated on the
lost geometry. It is left failing rather than loosened after seeing the number; whoever holds the
original `nbn_scf.chk` can say whether it still passes there. G3's premise — that the septet is
the reference worth pinning — is itself superseded by §4.1.

## 7. Out of scope

- Recovering the original CIF (needs MP access or the owner's copy).
- CASSCF orbital relaxation per spin state. The ordering here is CASCI on orbitals from the
  UHF(S=3) solution, which one would expect to favour S = 3 — but that is an expectation, not a
  check; state-specific CASSCF could reorder the active window and is not run here.
- Materials claims of any kind (finite diatomic, LANL2DZ ECP, 6-31G* N).

## 8. Caveats

- CASCI energies jump by up to ~25 mHa between neighbouring distances (e.g. 2.27 → 2.28 Å):
  which UHF orbitals land in the 14-orbital window changes with d. The ordering claim is made at
  each fixed d, never across d.
- DMRG B′ for S = 1 over-extrapolates by 0.34 mHa relative to exact FCI (stderr 0.29 mHa) — the
  D ≤ 300 dweight fit is not trustworthy at sub-mHa here; exact FCI is the reference.
