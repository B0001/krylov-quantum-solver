# SPEC: ODMD — ground-state energy from the survival amplitude alone

**Status:** IMPLEMENTED — gates G1–G4 green (`tests/test_odmd_spec.py`).

---

## 1. Goal

The complex survival amplitude `s_k = <phi_0| e^{-i k tau H} |phi_0>` — the *first row* of the
overlap matrix the QKSD pipeline already measures, a 1-D time series of `K` numbers — is by itself
sufficient to recover the FCI ground-state energy, via observable dynamic mode decomposition
(ODMD): Hankel-matrix DMD with an SVD-truncated pseudoinverse. Because no Hamiltonian matrix
element is ever measured, the lambda-scaled sampling noise that dominates KQD (variance ~ the
Pauli 1-norm; see `SPEC_msd_sampling.md`) vanishes entirely, so at a **matched number of measured
elements and shots per element** the ODMD median error beats KQD. Falsifiable: if the DMD
eigenphases do not converge to FCI, or noisy ODMD does not beat noisy KQD at matched budget, the
claim dies.

## 2. Background and honest framing

- Reproduction of ODMD — Shen, Klymko, Sud, Williams-Young, de Jong & Van Beeumen,
  `arXiv:2306.01858` ("Estimating eigenenergies from quantum evolution: a dynamic mode
  decomposition approach"). Their headline is exactly the noise robustness we gate on: the DMD
  least-squares + SVD truncation tolerates overlap noise that destroys Krylov GEVP methods.
- **What we can claim if gates pass:** the validated repo stack reproduces ODMD's central results
  on dense-diagonalizable systems — convergence to FCI from overlaps alone, and a quantified
  matched-budget sampling advantage over our own validated KQD implementation (`msd.py`'s KQD arm)
  under the same per-element shot-noise model used in `SPEC_msd_sampling.md`.
- **What we cannot claim:** no hardware demonstration (exact-statevector overlaps, idealized
  i.i.d. Gaussian shot noise); no quantum advantage at this scale; **no variational bound** —
  unlike QKSD's Ritz values, DMD eigenphases can dip *below* FCI at small K (a recorded finding,
  cf. CMX(2) in `SPEC_moment_pds.md`). We use the **complex** signal (both Hadamard-test
  quadratures, as our hardware path already measures); the paper's cheaper Re-only variant is out
  of scope. Depths here (K <= 20) resolve these small systems; larger gaps/overlap structure may
  need deeper K.

## 3. Approach

Work in the centered frame `H - mu` with `mu` the midpoint of the HF-reachable spectrum and
`tau = pi / W` (identical conventions to `msd.py`, so the KQD comparison arm is the already
validated one). Then:

1. `s_k = <phi_0| e^{-i k tau (H - mu)} |phi_0>`, `k = 0..K-1` (exact statevector,
   `expm_multiply`).
2. Hankel data matrices `X[i,j] = s_{i+j}`, `X'[i,j] = s_{i+j+1}` with `d = K//2` rows.
3. Reduced DMD operator `A = U_r^H X' V_r Sigma_r^{-1}` from the SVD of `X` truncated at
   `sigma_i > delta * sigma_max` — the truncation is the noise-robustness mechanism (G4).
4. Eigenvalues of `A` are `~ e^{-i E_n tau}`; energies `E_n = -arg(lambda_n)/tau` (the centered
   frame keeps `|E_n| tau <= pi/2`, so no phase wrapping). Keep near-unimodular modes
   (`||lambda|-1| < 0.2`) and take the minimum — the ground-state estimate.

**Reference:** FCI via `mh.ground_state_energy()` (the repo's vetted reference). **Comparator:**
the KQD sampler in `msd.py` (`build_msd_problem` / `sample_ground_energy(..., "kqd", ...)`),
whose per-element noise scales are pinned by `tests/test_msd_sampling_spec.py`. Matched budget
means matched *measured element count*: ODMD at K=16 measures 16 overlap elements; KQD at n=8
measures 8 overlap + 8 Hamiltonian elements, each at the same `shots`.

## 4. Public interface

```
odmd.build_odmd_problem(mh, n) -> ODMDProblem
    # fields: n, tau, mu, s (complex survival amplitudes), dim, offset, ref
    #   (ref = exact ground energy in the centered frame; ref + offset = total FCI)
odmd.odmd_energy(s, tau, svd_threshold=1e-10, mod_window=0.2) -> (energy, rank)
    # centered-frame ground-state estimate from the raw signal; rank = kept SVD rank
odmd.sample_odmd_energy(prob, shots, seed, n=None, svd_threshold=None) -> float
    # one shot-noisy estimate; sigma_s and the default noise-aware truncation (5*sigma_s)
    # follow the msd.py Hadamard-test conventions; s_0 = 1 exactly
```

Top-level module `odmd.py` (a method rung, like `msd.py` / `rodeo.py` — not package API).

## 5. Acceptance criteria (validation gates)

All in `tests/test_odmd_spec.py` (test-first). Median over noise seeds, as in
`SPEC_msd_sampling.md`.

- **G1 — overlaps alone recover FCI.** Noiseless ODMD at K=20 satisfies
  `|E_odmd + offset - E_fci| < 1e-5 Ha` on H2, H4 chain, and N2 CAS(6,6).
- **G2 — convergence with depth + the non-variational boundary.** On N2 CAS(6,6):
  `|err(K=16)| < |err(K=8)|` and `|err(K=16)| < 1e-5 Ha`; **and** the K=8 estimate falls at least
  1e-4 Ha *below* FCI — documenting that ODMD, unlike QKSD, carries no variational floor.
- **G3 — matched-budget advantage over KQD (DEFINITION OF DONE).** N2 CAS(6,6), 100 noise seeds,
  ODMD(K=16, 16 elements) vs KQD(n=8, 16 elements): median `|E - E_fci|` for ODMD beats KQD by
  > 15x at 1e4 shots and > 4x at 1e5 shots (measured ~57x and ~10x in the feasibility probe).
- **G4 — SVD truncation is the mechanism.** At 1e4 shots the noise-aware truncation keeps the
  median error < 5 mHa, and removing it (`svd_threshold = 1e-10`) inflates the median by > 20x
  (measured ~900x).

## 6. Implementation plan (test-first)

1. `tests/test_odmd_spec.py` encoding G1–G4 (RED — `odmd` does not exist yet).
2. `odmd.py`: minimum code — reuse `expm_multiply` overlaps and the `msd.py` frame/noise
   conventions; the KQD arm is imported from `msd.py`, not reimplemented.
3. `make gates` to green.

## 7. Out of scope

- The paper's Re-only signal variant (halves the Hadamard tests; needs +/-E disambiguation).
- Excited states / gaps from the higher DMD modes — done as the follow-up
  [`SPEC_odmd_excited.md`](SPEC_odmd_excited.md) (noise-edge thresholding + the visibility law).
- Trotterized evolution circuits, device noise models, hardware runs.
- Larger-than-dense-diagonalizable systems (no reference).

## 8. Caveats and risks

- **R1 — mode selection under noise:** a spurious near-unimodular DMD eigenvalue below the true
  ground would poison the min-phase selection. Mitigation: the noise-aware SVD truncation (5x the
  per-element sigma) plus the modulus window; G3/G4's 100-seed medians would catch a failure.
- ODMD is **not variational** (G2 documents it) — never quote an ODMD number as a bound.
- The advantage is quantified against *this repo's* KQD noise model (LCU/Hadamard variance
  ~ lambda^2/shots per H element). A different H-measurement strategy (e.g. MSD itself) narrows
  the gap; the honest comparison set is recorded in `BACKLOG.md`.

## 9. Deliverables

- `odmd.py` — the method (problem builder, estimator, noisy sampler + `__main__` demo).
- `tests/test_odmd_spec.py` — gates G1–G4.
- `BACKLOG.md` entry with the measured numbers and findings.
