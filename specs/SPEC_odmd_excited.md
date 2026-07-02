# SPEC: Excited-state ODMD — the gap from the same signal, via noise-edge thresholding

**Status:** IMPLEMENTED — gates G1–G4 green (`tests/test_odmd_excited_spec.py`).

---

## 1. Goal

The *same* K-number survival amplitude that ODMD uses for the ground state
(`SPEC_odmd.md`) also carries the low-lying excited spectrum in its higher DMD eigenphases —
**no additional measurements** — *provided* the SVD truncation is an **absolute noise-edge
cutoff** `sigma_i > c * sigma_s * (sqrt(d) + sqrt(m))` (the random-matrix / Marchenko–Pastur
scale of the Hankel noise) instead of `SPEC_odmd.md`'s relative `5 sigma_s * sigma_max` floor.
The relative floor is inflated by the dominant ground mode (`p_0 ~ 0.95`) until it swallows the
excited singular value `~ p_1 * sqrt(d m)`; the absolute edge keeps any mode whose amplitude
clears the noise, giving a testable **visibility law**: mode `n` is recoverable iff
`p_n * sqrt(d m) > c * sigma_s * (sqrt d + sqrt m)` — so *signal depth K buys excited-state
visibility as ~ sqrt(K)*. Falsifiable at both ends: the eigenphases must match the FCI spectrum,
and the visibility onset must move with K exactly as predicted.

## 2. Background and honest framing

- Builds on `SPEC_odmd.md` (arXiv:2306.01858 — ground state only). Excited-state ODMD variants
  in the literature use *multiple* observables/references; here we test what the **single**
  survival amplitude supports. Absolute singular-value hard-thresholding at the noise edge is
  standard random-matrix practice (Gavish–Donoho-style); we claim no priority on either
  ingredient — the contribution is the *quantified visibility boundary* on this repo's validated
  stack, and the head-to-head against our QKSD `solve_excited` under the matched noise model of
  `SPEC_qksd_noise.md`.
- **What we can claim if gates pass:** from one 48-element overlap signal, the first excitation
  gap of H₄ to a few mHa at 10⁵ shots/element — an order of magnitude below noisy QKSD
  `solve_excited` at *its own best depth* — plus the recorded law for when an excited mode is
  visible at all.
- **What we cannot claim:** Hankel noise entries are *correlated* (antidiagonal structure), so
  the MP edge with `c = 1.2` is a calibrated heuristic, not a theorem. Exact statevector,
  idealized i.i.d. per-element shot noise, tiny systems, no variational bound on any eigenphase.
  States with `p_n` below the visibility law are simply invisible — this method cannot see dark
  or weakly-overlapped states (same physics as the rodeo peak-height and SKQD findings).

## 3. Approach

Same centered frame and signal as `SPEC_odmd.md`. New pieces:

1. **Spectrum extraction:** keep ALL near-unimodular DMD eigenvalues (not just the min phase);
   sort by energy; amplitudes `a_n` from a Vandermonde least-squares refit of the signal
   (`s_k ~ sum_n a_n lambda_n^k`) for optional filtering/diagnostics.
2. **Noise-edge truncation:** `cutoff = c * sigma_s * (sqrt d + sqrt m)`, `c = 1.2` (calibrated;
   G3/G4 would fail if miscalibrated), replacing the relative floor when a noise scale is known.
3. **References:** FCI spectrum restricted to HF-reachable states (dense diagonalization —
   exactly the reference set of `SPEC_qksd_excited.md`). Noisy comparator: `solve_excited` with
   the solver's own Hermitian per-element noise at matched `sigma_s` (generous to QKSD: its H
   elements get the same sigma, with no lambda factor).

## 4. Public interface

```
odmd.noise_edge(sigma, d, m, c=1.2) -> float
    # absolute Hankel singular-value cutoff for per-element noise sigma
odmd.odmd_spectrum(s, tau, cutoff=0.0, mod_window=0.2, amp_floor=0.0)
    -> (energies, amplitudes, rank)      # ascending centered-frame eigenphases + |a_n|
odmd.sample_odmd_spectrum(prob, shots, seed, n=None, c=1.2, amp_floor=0.0) -> np.ndarray
    # one shot-noisy spectrum draw; msd.py noise conventions, noise-edge truncation
```

`odmd_energy` / `sample_odmd_energy` keep their exact `SPEC_odmd.md` semantics (gates pinned);
both now share a private `_dmd_modes` core.

## 5. Acceptance criteria (validation gates)

All in `tests/test_odmd_excited_spec.py`; noisy gates use 100-seed medians.

- **G1 — same signal, whole low-lying spectrum (noiseless).** `|E_1 − E_1^FCI| < 1e-5 Ha` and
  `|gap − gap^FCI| < 1e-5 Ha` on H₄ (K=24) and N₂ CAS(6,6) (K=48).
- **G2 — depth is the excited-state resource.** On H₄ the K=16 signal has gap error > 1 mHa
  (measured 6.7) while K=24 is < 1e-5 Ha — excited eigenphases need deeper K than the ground
  state (the `SPEC_qksd_excited.md` finding, reproduced in signal space).
- **G3 — noisy gap beats QKSD at matched noise (DEFINITION OF DONE).** H₄, K=48, 10⁵
  shots/element, c=1.2: median |gap error| < 10 mHa (measured 5.8), the excited mode resolved in
  ≥ 95% of seeds, and the median beats noisy QKSD `solve_excited` at BOTH M=16 and M=24 by > 10×
  (measured ~31× and ~21×; QKSD is noiselessly converged at M=24, so noise — not depth — is its
  limit).
- **G4 — the threshold is the mechanism, and visibility follows the law.** Same setting:
  (a) the relative `5 sigma * sigma_max` floor of `SPEC_odmd.md` fails to resolve the excited
  mode in ≥ 90% of seeds (measured 100%); (b) at K=16 even the noise edge cannot see it
  (unresolved in ≥ 50% of seeds, measured 74% — `p_1 sqrt(dm)` is below the edge), while K=48
  resolves it in ≥ 95% — the sqrt(K) visibility onset.

## 6. Implementation plan (test-first)

1. `tests/test_odmd_excited_spec.py` encoding G1–G4 (RED — `odmd_spectrum` missing).
2. Extend `odmd.py`: factor `_dmd_modes`, add `noise_edge` / `odmd_spectrum` /
   `sample_odmd_spectrum`. Re-run `tests/test_odmd_spec.py` (pinned) after the refactor.
3. `make gates` to green.

## 7. Out of scope

- Multi-observable / multi-reference ODMD (richer signals raise `p_n`; a follow-up).
- States invisible under the visibility law (dark states; needs a better reference, not a
  better threshold).
- Higher excited states (E₂+ — same machinery, unvalidated here), Trotter circuits, hardware.

## 8. Caveats and risks

- **R1 — edge calibration:** the MP constant for *structured* (Hankel) noise is empirical;
  c too small admits noise modes below E₀, c too large re-creates the relative-floor failure.
  G3/G4's seed medians bound both failure modes at c=1.2.
- The QKSD comparison inherits `SPEC_qksd_noise.md`'s idealized noise model, applied equally to
  both methods (and favourably to QKSD).
- Never quote an eigenphase as a variational bound (`SPEC_odmd.md` G2).

## 9. Deliverables

- `odmd.py` — `noise_edge`, `odmd_spectrum`, `sample_odmd_spectrum` (+ shared `_dmd_modes`).
- `tests/test_odmd_excited_spec.py` — gates G1–G4.
- `BACKLOG.md` entry with measured numbers and the visibility law.
