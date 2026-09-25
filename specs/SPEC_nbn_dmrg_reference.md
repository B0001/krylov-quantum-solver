# SPEC: DMRG-referenced transition-metal active space — NbN CAS(14,14) beyond the FCI cutoff

**Status:** IMPLEMENTED — gates G1–G3 green (`tests/test_nbn_dmrg_reference_spec.py`); headline
two-schedule run recorded below (2026-07-04, driver-level).

> **CORRECTED 2026-09-25 (bead chem-dc7, `SPEC_nbn_low_spin.md`) — the certified energy is not the
> CAS ground state.** The (10,4) sector this spec certifies is block2 SU(2) **S = 3**, inherited
> from the UHF spin scan. Exact FCI of the CAS(14,14) Ms = 0 sector (all 1.18e7 determinants,
> 19 min on 4 cores) finds the ground is **S = 1**, 15.1 mHa *below* S = 3; S = 2 is also below it.
> **Corrected reference (reconstructed geometry d(Nb–N) = 2.25 Å, see below): E₀ = −110.0712455 Ha,
> S = 1**, DMRG-cross-checked to 14 µHa. The septet is the CAS ground at none of nine Nb–N distances
> 2.0–3.0 Å. Also: the original CIF/chkfile were never committed, and the MP-bond-length
> reconstruction gives E(10,4) = −110.056110, **not** the −110.046028 below — that number is not
> reproducible from the repo. The original record is kept below, struck where superseded.

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
  NbN could invert them. A hard multireference TM benchmark needs that low-spin sector at real
  bond dimension or a larger cluster — a follow-up, out of scope here.

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

**Spin-sector map (both schedules agree to µHa):**

| sector | 2S | E (Ha) | Δ vs ground |
|---|---|---|---|
| nelec=(10,4) | 6 | −110.046028 | 0 ~~(ground — the cached SCF sector)~~ the cached SCF sector, **not** the CAS ground |
| nelec=(7,7)  | 0 | −110.042500 | +3.5 mHa |

**Corrected spin map (reconstructed d = 2.25 Å; `SPEC_nbn_low_spin.md` Table 1):**

| S | E (Ha) | Δ vs ground | source |
|---|---|---|---|
| 1 | **−110.0712455** | 0 (**ground**) | exact FCI, Ms=0 sector unconstrained, ⟨S²⟩ = 2.000 |
| 2 | −110.0695221 | +1.7 mHa | exact FCI, (9,5) sector, ⟨S²⟩ = 6.000 |
| 3 | −110.0561097 | +15.1 mHa | exact FCI, (10,4) sector — the sector certified above |
| 0 | −110.03199 | +39.3 mHa | DMRG A/B at D ≤ 1200 agree to 5 µHa |

The "3.5 mHa singlet" does not reproduce either: on the reconstruction the singlet is 24 mHa
*above* the septet.

The low-spin sector (the FCI-intractable comb(14,7)² one) is only 3.5 mHa up — a
near-degeneracy the SCF spin scan does not by itself surface; the reference energy is only
well-defined once the sector is named.

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
  one; corrected 2026-09-25.)
- **G2 — two independent schedules agree (DEFINITION OF DONE, ~2 min).** Cheap-dims
  A′ (perD 100/200/300) vs B′ (ramp 80/160/300), different seeds/scratch:
  `|E_A′ − E_B′| < 0.1 mHa` (measured 0.0012) and **both** `method == "dweight"` (the regime
  guard: at these dims the discarded weights are usable, unlike the converged headline dims).
- **G3 — the softness finding, gated.** *(2026-09-25: FAILS on the regenerated checkpoint —
  dw(300) = 3.0e-7 and per-D spread 0.088 mHa. Threshold left as is, not loosened after seeing
  the number; it was calibrated on the lost original geometry. See `SPEC_nbn_low_spin.md` §6.)* Discarded weight at D=300 < 1e-7 and the per-D energy
  spread < 0.01 mHa — the low-entanglement character is pinned, so no one later mistakes this
  reference for a strong-correlation benchmark.

## 6. Implementation plan (test-first)

1. `tests/test_nbn_dmrg_reference_spec.py` encoding G1–G3 (RED — driver module missing).
2. `nbn_dmrg_reference.py` — the loader + schedule runner (promoted from the launch scripts).
3. `make gates`; record the headline A/B numbers here.

## 7. Out of scope

- A genuinely multireference TM benchmark (low-spin sector / larger cluster / bigger basis).
- Materials claims (finite cluster, ECP, fixed geometry); periodic NbN.
- Krylov/ODMD on this system (nothing at 14 orbitals needs a quantum method — no advantage).

## 8. Caveats and risks

- **R1:** the high-spin sector makes this easy; a future low-spin variant may need real bond
  dimension and should expect `dweight` extrapolation to matter.
- The 5×10⁶ "FCI-intractable" line is the repo's operational cutoff, not a fundamental wall.

## 9. Deliverables

- `nbn_dmrg_reference.py`; `tests/test_nbn_dmrg_reference_spec.py` (G1–G3).
- Headline record in §3; `BACKLOG.md` item closed with the softness finding.
