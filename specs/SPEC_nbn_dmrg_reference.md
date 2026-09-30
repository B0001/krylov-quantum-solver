# SPEC: DMRG-referenced transition-metal active space — NbN CAS(14,14) beyond the FCI cutoff

**Status:** IMPLEMENTED — gates G1, G2, G3, G3b green (`tests/test_nbn_dmrg_reference_spec.py`);
headline two-schedule run recorded below (2026-07-04, driver-level).

> **CORRECTED 2026-09-26 (bead chem-bbi) — the original files were not lost; the reconstruction
> targeted the wrong bond length.** chem-dc7 believed `data/nb_structures/NbN_mp-2634.cif` and
> `data/nbn_scf.chk` were unrecoverable and rebuilt a diatomic at MP's reported d(Nb–N) = 2.25 Å,
> which does not reproduce −110.046028 (see the struck 2026-09-25 note below). The real files were
> in fact present on this container's persistent `data/` volume all along (just gitignored, so
> never committed) and are now vendored at `specs/nbn_scf_reference.chk` / `specs/nbn_mp-2634.cif`
> — `nbn_dmrg_reference.ensure_reference_data()` materializes them into the `data/` paths
> `load_nbn_cas` expects on a fresh clone (wired into `load_nbn_cas`'s default-path branch, so
> every caller gets it, not just this spec's gate). Running the pipeline on these real files
> reproduces −110.046028 to 6 decimal places, closing the "unreproducible" finding.
>
> The 2.25 Å reconstruction's mismatch is now explained rather than mysterious: the real
> checkpoint's stored Cartesian coordinates put Nb–N at **3.7255 Å**, not ~2.25 Å. Root cause:
> `benchmark_nbn.py`'s `ground_state_mf` builds its `atom_str` from ASE's raw, unwrapped
> `atom.position` values, not the minimum-image-reduced separation — `atoms.get_distance(0, 1)`
> (no MIC) gives 3.7255 Å, `atoms.get_distance(0, 1, mic=True)` gives 2.2464 Å, matching MP's
> reported bond length. The reconstruction was simply solving a different, longer-bond molecule;
> its Table 1–3 numbers in `SPEC_nbn_low_spin.md` are not a probe of the real geometry.
>
> **The sector conclusion still holds, now on the real geometry, with tighter numbers.** Exact
> FCI (Ms=0 sector, all 1.18e7 determinants, 1010.6 s / 8 cores) finds the CAS ground is **S = 1**,
> not the committed S = 3: E(S=1) = −110.04756504307636 Ha (⟨S²⟩ = 2.0000095761529892), only
> **1.538 mHa** below E(S=3) = −110.04602823118735 Ha (DMRG A′, matching the committed
> −110.046028 to 6 decimals) — not the 15.1 mHa gap recorded against the wrong-bond-length
> reconstruction. S = 0 (nelec=(7,7)) is above both, at −110.04250013542108 Ha. See G3b below and
> `results/nbn_low_spin/runs.jsonl` (2026-09-26T15:07:42) for provenance. **chem-g1i (filed
> 2026-09-26) tracks that `SPEC_nbn_low_spin.md`'s own Tables 1–3 and its G5 gate (`> 5 mHa`
> margin, calibrated against the 15.1 mHa figure) need the same correction** — out of this bead's
> scope, since that spec belongs to the already-closed chem-dc7.
>
> ~~CORRECTED 2026-09-25 (bead chem-dc7, `SPEC_nbn_low_spin.md`) — the certified energy is not the
> CAS ground state.~~ The (10,4) sector this spec certifies is block2 SU(2) **S = 3**, inherited
> from the UHF spin scan — that part still holds. ~~Exact FCI of the CAS(14,14) Ms = 0 sector (all
> 1.18e7 determinants, 19 min on 4 cores) finds the ground is S = 1, 15.1 mHa below S = 3; S = 2 is
> also below it. Corrected reference (reconstructed geometry d(Nb–N) = 2.25 Å, see below):
> E₀ = −110.0712455 Ha, S = 1, DMRG-cross-checked to 14 µHa. The septet is the CAS ground at none
> of nine Nb–N distances 2.0–3.0 Å. Also: the original CIF/chkfile were never committed, and the
> MP-bond-length reconstruction gives E(10,4) = −110.056110, not the −110.046028 below — that
> number is not reproducible from the repo.~~ The −110.046028 figure IS reproducible from the real
> geometry (chem-bbi, above); the 2.25 Å reconstruction and everything computed on it (15.1 mHa
> gap, E₀ = −110.0712455 Ha, the nine-distance scan) described a different, longer-bond molecule.

---

## 1. Goal

Close the backlog item: for a transition-metal active space **large enough that FCI is
intractable** (NbN 2-atom cluster, CAS(14,14): the half-filling Sz=0 sector is
comb(14,7)² ≈ 1.18×10⁷ determinants, beyond the repo's 5×10⁶ FCI cutoff — a black-box FCI must
handle it), DMRG gives a **converged correlation energy** — certified not by a single run but by
**two independent sweep schedules agreeing** and by bond-dimension convergence with usable
discarded weights. The claim is false if the schedules disagree at the mHa level or the
extrapolation leaves the discarded-weight regime.

(The cached ground sector is the high-spin nelec=(10,4); its *own* fixed-sector count is
comb(14,10)·comb(14,4) ≈ 1.0×10⁶. The FCI-intractability is a property of the full active-space
problem, and DMRG's structural advantage is precisely that it never enumerates the whole space —
it works one Sz sector at a time.)

## 2. Background and honest framing

- Builds on the pinned DMRG plumbing (`SPEC_singleramp.md`; `dmrg_energy_extrapolated` with
  `protocol="perD"` vs `protocol="ramp"` — genuinely independent schedules: separate converged
  runs per D vs one ramping run) and the cached spin-scanned SCF (`data/nbn_scf.chk`,
  `benchmark_nbn.py`).
- **What we can claim if gates pass:** a reference number, not a materials claim (finite
  2-atom cluster, LANL2DZ ECP, fixed geometry — the backlog item's own framing): NbN CAS(14,14)
  ~~ground energy~~ energy **E = −110.046028 Ha** of its high-spin **S = 3 excited spin state**
  (not the ground — see the correction above), reproducible to sub-µHa across independent
  schedules on the original, uncommitted geometry.
- **What we cannot claim / the recorded findings:** (i) this CAS is a **soft DMRG target** —
  the high-spin nelec=(10,4) (2S=6) sector is low-entanglement: D=400 converges to sub-nHa
  (discarded weight ~1e-9), even D≤300 sits within 1 µHa. "FCI-intractable by determinant count"
  did not mean "strongly correlated". (ii) **A near-degeneracy the SCF spin scan alone does not
  surface:** the low-spin nelec=(7,7) sector lies just **3.5 mHa above** the (10,4) ground
  (E₇,₇ = −110.042500 vs E₁₀,₄ = −110.046028, both schedule-agnostic to µHa) — so the "reference
  energy" is only meaningful once the sector is stated, and a thermally/ligand-field-perturbed
  NbN could invert them. **UPDATE (chem-dc7/chem-g1i, `SPEC_nbn_low_spin.md`; re-verified fresh on
  new hardware by chem-csf, 2026-09-29):** the (7,7) sector was run at real bond dimension — it is
  soft too (dw → 3e-10 by D=1200, two independent schedules agree to 0.3 nHa), so this CAS is not
  a hard multireference benchmark in *any* sector checked so far; see `SPEC_nbn_low_spin.md` §0/§3
  Table 2/G4 for the numbers. This does not rule out a genuinely hard sector existing in a larger
  cluster or basis, which remains unexplored.

## 3. Approach

Restore the cached ground-spin SCF, build CAS(14,14) integrals, run
`dmrg_energy_extrapolated` twice: schedule A = perD, schedule B = ramp, with different bond-dim
ladders, seeds, and scratch dirs. CI gates use cheap dims (≤300, ~2 min, still in the
discarded-weight regime); the headline D≤1200 numbers are a driver-level record (as in
`SPEC_hchain_largen2.md`).

**Headline record (driver `nbn_dmrg_reference.py`, 2026-07-04, 16 GB laptop):**

| schedule | dims | protocol | E (Ha) | stderr | method |
|---|---|---|---|---|---|
| A | 400/800/1200 | perD | −110.04602843 | ~0 | invD (weights ~1e-9–1e-13: nothing left to extrapolate) |
| B | 300/600/1200 | ramp | −110.04602846 | ~0 | invD |

|E_A − E_B| = 3×10⁻⁸ Ha — five orders below the 1 mHa gate. (The `invD` fallback here signals
*convergence*, the opposite of the failure mode that killed `SPEC_hchain_largen.md`.)

**Spin-sector map, real vendored geometry (chem-bbi, 2026-09-26 — see the correction above):**

| S | nelec | E (Ha) | Δ vs S=1 ground | source |
|---|---|---|---|---|
| 1 | (8,6) | **−110.04756504307636** | 0 (**ground**) | exact FCI, Ms=0 sector unconstrained, ⟨S²⟩ = 2.0000095761529892, 1.18e7 dets, 1010.6 s / 8 cores |
| 1 | (8,6) | −110.04756608985925 | −1.05 µHa (extrapolation overshoot, cf. §8 caveat below) | DMRG A′ cross-check, agrees with FCI to ~1 µHa (G3b) |
| 3 | (10,4) | −110.04602823118735 | +1.538 mHa | DMRG A′ — the cached SCF sector G1/G2/G3 certify, matches the historical committed −110.046028 to 6 decimals |
| 0 | (7,7) | −110.04250013542108 | +5.065 mHa | DMRG A′ |

S = 2 (nelec=(9,5)) has not been remeasured on the real geometry (the −110.0695221 figure below
is the wrong-bond-length reconstruction's number, not this geometry's) — tracked by chem-g1i
alongside the rest of `SPEC_nbn_low_spin.md`'s tables.

The gap between the CAS ground (S=1) and the committed sector (S=3) is **1.538 mHa**, not the
15.1 mHa recorded against the reconstruction below. It is still a real, signed gap (S=1 below
S=3) — the qualitative finding ("the committed sector is not the ground") survives; only its
size does not.

**Corrected spin map, WRONG-BOND-LENGTH RECONSTRUCTION (d = 2.25 Å nominal, but see the
2026-09-26 correction above: the real separation implied by this reconstruction's own pipeline
is 3.7255 Å, not 2.25 Å — `benchmark_nbn.py`'s unwrapped-position bug. Kept for the record, not
as a reference; `SPEC_nbn_low_spin.md` Table 1):**

| S | E (Ha) | Δ vs ground | source |
|---|---|---|---|
| 1 | −110.0712455 | 0 (ground of *this* molecule) | exact FCI, Ms=0 sector unconstrained, ⟨S²⟩ = 2.000 |
| 2 | −110.0695221 | +1.7 mHa | exact FCI, (9,5) sector, ⟨S²⟩ = 6.000 |
| 3 | −110.0561097 | +15.1 mHa | exact FCI, (10,4) sector |
| 0 | −110.03199 | +39.3 mHa | DMRG A/B at D ≤ 1200 agree to 5 µHa |

These four numbers describe a real diatomic — just not the one this repo committed. They are not
comparable to the real-geometry table above; do not average or interpolate between the two.

## 4. Public interface

```
nbn_dmrg_reference.load_nbn_cas(norb=14, nelec_cas=14, chk="data/nbn_scf.chk")
    -> (h1, eri, nelec, e_core)      # restored-SCF CAS integrals (no SCF re-run)
nbn_dmrg_reference.py --schedule A|B|cheap    # one schedule per process; prints E/stderr/method
```

## 5. Acceptance criteria (validation gates)

`tests/test_nbn_dmrg_reference_spec.py` — pyscf + block2 only (no qiskit imports; `make gates`
process isolation applies).

- **G1 — beyond-FCI guard + sector pin (instant).** The CAS(14,14) determinant count exceeds
  the repo's 5×10⁶ FCI cutoff, the chkfile restores without an SCF run, and the SCF sector
  is the spin-scanned nelec=(10,4). (~~ground sector~~ — it is the UHF ground spin, not the CAS
  one; corrected 2026-09-25.) A fresh clone now passes this too: `ensure_reference_data()`
  materializes the vendored `specs/nbn_scf_reference.chk` / `specs/nbn_mp-2634.cif` into the
  `data/` paths `load_nbn_cas` expects (chem-bbi, 2026-09-26).
- **G2 — two independent schedules agree (DEFINITION OF DONE, ~2 min).** Cheap-dims
  A′ (perD 100/200/300) vs B′ (ramp 80/160/300), different seeds/scratch:
  `|E_A′ − E_B′| < 0.1 mHa` (measured 0.0012) and **both** `method == "dweight"` (the regime
  guard: at these dims the discarded weights are usable, unlike the converged headline dims).
- **G3 — the softness finding, gated.** Discarded weight at D=300 < 1e-7 and the per-D energy
  spread < 0.01 mHa for the (10,4) sector — the low-entanglement character of *that sector* is
  pinned, so no one later mistakes it for a strong-correlation benchmark. (The 2026-09-25 failure
  — dw(300) = 3.0e-7 — was measured on chem-dc7's wrong-bond-length reconstruction, a different
  molecule; on the real, vendored geometry this gate passes and reproduces the original headline
  numbers exactly — see the correction note above.)
- **G3b — the true ground is pinned (chem-bbi, 2026-09-26).** The G1/G3 sector (10,4)/S=3 is not
  the CAS ground. DMRG A′ on nelec=(8,6)/S=1 must land within 0.5 mHa of the exact-FCI reference
  −110.04756504307636 Ha (measured agreement ~1 µHa) and strictly below the (10,4) sector's own
  energy (measured gap 1.538 mHa).

## 6. Implementation plan (test-first)

1. `tests/test_nbn_dmrg_reference_spec.py` encoding G1–G3 (RED — driver module missing).
2. `nbn_dmrg_reference.py` — the loader + schedule runner (promoted from the launch scripts).
3. `make gates`; record the headline A/B numbers here.

## 7. Out of scope

- A genuinely multireference TM benchmark (larger cluster / bigger basis — the low-spin (7,7)
  sector of *this* CAS was checked, chem-dc7/chem-g1i/chem-csf, and found soft too;
  `SPEC_nbn_low_spin.md` §0/§3 Table 2/G4).
- Materials claims (finite cluster, ECP, fixed geometry); periodic NbN.
- Krylov/ODMD on this system (nothing at 14 orbitals needs a quantum method — no advantage).

## 8. Caveats and risks

- **R1 (checked, chem-csf 2026-09-29):** the high-spin sector makes this easy; the low-spin (7,7)
  variant was expected to possibly need real bond dimension, but converges just as softly
  (`SPEC_nbn_low_spin.md` §3 Table 2) — `dweight` extrapolation is usable at D≤300 in both sectors,
  it just sits closer to its regime boundary for (7,7) (dw(300)≈1e-7 vs ≪1e-7 for (10,4)).
- The 5×10⁶ "FCI-intractable" line is the repo's operational cutoff, not a fundamental wall.
- **R2 (chem-bbi):** cheap DMRG A′ (dweight-extrapolated, D ≤ 300) can overshoot below the exact
  variational bound by ~1 µHa on the S=1 sector (measured; cf. `SPEC_nbn_low_spin.md` §8's 0.34
  mHa overshoot for DMRG B′ on the same sector) — small relative to G3b's 0.5 mHa tolerance, but
  a reminder that `dweight`-regime numbers are extrapolated, not raw sweep energies.
- **R3 (chem-bbi):** `data/` being gitignored means any *other* uncommitted checkpoint under that
  name on a given machine will silently take precedence over the vendored one in
  `ensure_reference_data()` (it only materializes when the target path is missing). This is by
  design (don't clobber a developer's own in-progress `data/`), but means "passes here" is not
  automatically "passes everywhere" until confirmed on a genuinely fresh clone.
- **R4 (chem-owr/chem-y6w, 2026-09-27): the vendored chk itself is metastable.** G1, G2, G3, and
  G3b above all load their CAS integrals from the SAME vendored `specs/nbn_scf_reference.chk`,
  which sits at one of (at least) two genuinely, separately-stable UHF minima in the (10,4)/S=3
  sector — 0.548 mHa apart in SCF energy, 4.84 mHa apart after CASCI(14,14)/FCI. A deterministic
  fresh regeneration (see below) finds the *other*, lower-energy minimum on the Linux container
  (on macOS arm64 it finds this one; BLAS-build dependent, `SPEC_nbn_scf_determinism.md` §8) — not
  a false-instability artifact, both minima independently pass their own `stability()` check.
  **Decision (chem-y6w):** keep this chk as the pinned reference — re-deriving G1-G3b's numbers on
  the other minimum is filed separately (**chem-czw**), since every sector shares these same
  orbitals (`load_nbn_cas`'s `CASCI.get_h1eff`/`get_h2eff`), so it is not a one-gate fix. What *is*
  fixed here: the regeneration path's own determinism was previously a coin flip (ambient-thread-
  dependent BLAS non-associativity inside pyscf's `mf.stability()`); `benchmark_nbn._tight_scf`
  now pins `pyscf.lib.num_threads(1)` for its duration, verified bit-identical across ambient
  thread counts and repeats (`specs/SPEC_nbn_scf_determinism.md`,
  `tests/test_nbn_scf_determinism_spec.py`). This does **not** make a fresh regeneration reproduce
  this chk — it makes it reproducibly find the *other* minimum instead, every time.
  **chem-czw resolution (2026-09-27):** re-derived; the vendored chk stays. All four spin sectors
  re-run on the other minimum's orbitals (`SPEC_nbn_low_spin.md` §0b, Table 4) give strictly
  *higher* (less favorable) CASCI(14,14) energies than this spec's own G1-G3b numbers, by 1.44-4.84
  mHa depending on sector — a substantive reason to keep this chk, not just the compute-cost
  argument chem-y6w gave. G1-G3b themselves are unchanged and not re-pointed at the other minimum;
  the cross-check gate lives in `tests/test_nbn_czw_lower_minimum_spec.py`, own process.

## 9. Deliverables

- `nbn_dmrg_reference.py` (`ensure_reference_data`, `load_nbn_cas`, `run_schedule`);
  `tests/test_nbn_dmrg_reference_spec.py` (G1, G2, G3, G3b).
- `specs/nbn_scf_reference.chk`, `specs/nbn_mp-2634.cif` — the real, tracked reference inputs
  (chem-bbi, 2026-09-26), vendored alongside `specs/hchain_tdl_localized_table.csv` and
  `specs/hubbard_lieb_wu_table.csv` as this repo's convention for spec-gate fixtures.
- Headline record in §3; `BACKLOG.md` item closed with the softness finding (2026-07) and
  corrected sector-ground finding (chem-bbi, 2026-09-26).
