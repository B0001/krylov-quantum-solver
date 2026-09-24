# SPEC: Nb3(Br1-xIx)8 Alloy Magnetic Phase Diagram via Spin-ODMD

**Status:** v2.0 — IMPLEMENTED. G1, G2, G4 pass as written; **G3 FALSIFIED** and replaced by a
gate pinning the measured boundary. v1.0's claims are preserved verbatim in §9 alongside the
measurement that resolved each one.
**Depends on:** `odmd_spin.py`, `nb3x8_gaps.py` (`NB3X8_LT_BULK`, `dimer_cluster_integrals`)
**Gates:** `tests/test_nb3x8_alloy_spec.py` (G1, G2, G3, G3b, G4, G5) · `nb3x8_alloy.py` ·
`data/nb3x8_alloy_diagram.csv`

---

## 1. Goal

Locate the alloy fraction $x_c$ at which localized-moment (Heisenberg) superexchange stops being a
valid description of the interlayer magnetism of $\text{Nb}_3(\text{Br}_{1-x}\text{I}_x)_8$, by
running spin-ODMD on a **virtual-crystal** interpolation of the two endpoint cRPA parameter sets.

**Answer:** $x_c = 0.296953$ — 41% iodine. v1.0 predicted $(0.45, 0.65)$ and that is wrong by
$0.153$ in $x$. The reason is arithmetic that was available before any code was written: Heisenberg
is *already* 14.1% off at the bromide endpoint, so a 20% threshold has only ~6 points of headroom
and is crossed less than a third of the way across the series.

## 2. What this is and — importantly — what it is not

**This is the virtual crystal approximation (VCA).** Linearly interpolating ab-initio cRPA
integrals between two compounds models a fictitious averaged halide sitting on every anion site. It
is a *model prediction under a stated approximation*, not an ab-initio result for the alloy.
Specifically it contains:

- **no site disorder** — every Nb₃ trimer sees the same averaged environment;
- **no local halide configuration** — a real $x = 0.5$ crystal has Br-rich and I-rich local
  coordinations with different $t$ and $U$, and the exchange is nonlinear in those;
- **no lattice relaxation and no bowing** — VCA produces linear parameter trends *by construction*,
  so it cannot generate the bowing that real semiconductor and magnetic alloys almost always show;
- **no re-screening** — cRPA screening is not recomputed for the mixed compound.

v1.0 claimed this curve is *"wholly unmapped in the published literature"* and framed the output as
an ab-initio prediction. The first half is plausible; the second is an overclaim, and it is
withdrawn. A VCA interpolation being unpublished is not the same as an alloy cRPA calculation being
unpublished, and this repo would rather publish the boundary than the inflation. What is defensible:
**within this model, $x_c = 0.296953$, and the number is falsifiable by a supercell cRPA
calculation.** That is the claim on offer.

Everything `odmd_spin.py` disclaims is inherited: the **isolated** interlayer dimer only (no
in-plane kagome exchange, no band broadening — cf. the coordination correction in `nb3x8_gaps.py`),
density-density interactions only, exact statevector.

## 3. Approach (as built)

1. **Endpoints.** $\text{Nb}_3\text{Br}_8$ and $\text{Nb}_3\text{I}_8$ LT-bulk cRPA parameters
   straight from `NB3X8_LT_BULK` (arXiv:2501.10320 Table I). No new parameter tables.
2. **VCA interpolation.** $p(x) = (1-x)\,p_{\text{Br}} + x\,p_{\text{I}}$ for
   $p \in \{U_0, t, U_s\}$. v1.0's §3 interpolated only $t$ and $U$; $U_s$ must be interpolated too
   — it is the $J_{\text{H}}$ of v1.0's Heisenberg formula and it moves by 24% across the series.
3. **Spin-ODMD.** Rebuild the dimer through `dimer_cluster_integrals`, select the **half-filled**
   ground state by overlap with |HF⟩ (full-Fock diagonalization returns the one-electron bonding
   state for this model — `model_hamiltonians.py` warns about exactly this), kick with
   $S_{1z} - S_{2z}$, read the single line: $J(x)$.
4. **Heisenberg deviation.** $\Delta(x) = |J_{\text{Heis}}(x) - J(x)| / J(x)$ with
   $J_{\text{Heis}} = 4t^2/(U_0 - U_s)$; $x_c$ is where $\Delta$ first reaches $0.20$.

**The closed form that makes the kill airtight.** $\Delta$ depends on the three parameters only
through the single dimensionless group $a = \left(2t/(U_0 - U_s)\right)^2$:

$$\Delta(a) = \frac{2a}{\sqrt{1+4a}-1} - 1 \qquad\Longrightarrow\qquad \Delta = D \iff \left|\frac{2t}{U_0-U_s}\right| = \sqrt{D + D^2}$$

Both sides are linear in $x$, so $x_c$ is the root of a linear equation — an exact reference, not a
numerical estimate. Bisection and the closed form agree to $2.4\times10^{-15}$.

## 4. Public interface (`nb3x8_alloy.py`)

```python
alloy_parameters(x)                       -> dict[str, float]     # VCA (U0, t, Us)
alloy_integrals(x)                        -> ModelIntegrals
alloy_ground_state(x)                     -> (MolecularHamiltonian, np.ndarray)
alloy_exchange_constant(x, n=24)          -> (J_meV, local_moment_fraction)
alloy_heisenberg_exchange(x)              -> float
alloy_heisenberg_deviation(x)             -> float
critical_alloy_fraction(threshold=0.20)   -> float                # bisection
critical_alloy_fraction_closed_form(...)  -> float                # the exact reference
map_magnetic_phase_diagram(steps=21)      -> list[dict]
total_spin_operator()                     -> sparse S^2
```

**Three deliberate deviations from v1.0's signatures:**

- `alloy_integrals` returns a **`ModelIntegrals`**, not the bare
  `tuple[ndarray, ndarray, float, tuple[int,int]]`. That dataclass carries exactly those fields plus
  `norb`, it is what `dimer_cluster_integrals` already returns, and it is what every consumer in
  this repo takes. Re-deriving the cluster to satisfy a type annotation would be new code with no
  new check behind it.
- `alloy_exchange_constant` **drops `shots`**. `odmd_spin` is an exact-statevector path with no shot
  model; a silently-ignored `shots=` parameter is a lie in a signature. `krylov_dim` is renamed `n`
  to match `absorption_lines`, the primitive it forwards to.
- It returns `(J, local_moment_fraction)` rather than a bare `J` — the spectral weight
  $\lVert S_z|\psi_0\rangle\rVert^2$ comes out of the same call for free and v1.0 wanted it in the
  diagram anyway.

## 5. Acceptance gates (as gated, with the measured numbers)

- **G1 — Endpoint consistency (unchanged).** $J(0) = 119.1081$ meV, $J(1) = 245.9196$ meV. Against
  `odmd_spin`'s own gated analytic exchange: relative error $< 10^{-14}$ (the ODMD line *is* the
  exact singlet–triplet splitting, so this is machine precision, not a tolerance). Against the
  published 4-figure values 119.1 / 245.9 meV: $6.8\times10^{-5}$ and $8.0\times10^{-5}$ — the
  specced 0.1% holds with three orders of magnitude to spare, so the rounding in those figures never
  mattered. *Passes.*
- **G2 — Smooth, strictly monotonic, non-singular interpolation (unchanged).** $U_0$, $t$ and $U_s$
  are each **strictly decreasing** in $x$ ($t: -169.40 \to -218.20$, $U_0: 1186.60 \to 787.00$,
  $U_s: 342.00 \to 258.50$), so "strictly monotonic" is true for all three as signed quantities.
  The nuance v1.0 did not state, and the mechanism behind G3: $|t|$ *grows* while the
  charge-transfer scale $U_0 - U_s$ *shrinks*, so the deviation's single group $|2t/(U_0-U_s)|$
  rises $0.401 \to 0.826$ and both parameters push $\Delta$ the same way. $J(x)$ is smooth
  (max second difference $< 1\%$ of the curve's range) and strictly increasing
  $119.11 \to 245.92$ meV, with exactly one spectral line at every $x$ — no degeneracy collapse, no
  singular matrix. Local-moment fraction falls monotonically $0.890 \to 0.759$. *Passes.*
- **G3 — THE KILL (replaces v1.0's G3).** $x_c = \mathbf{0.296953}$ at the spec's own 20%
  threshold, pinned to $\pm 10^{-6}$ against the closed form of §3. The gate asserts explicitly that
  $x_c \notin (0.45, 0.65)$, so the falsification cannot silently decay back into a pass.
  $\Delta(x)$ is verified strictly increasing on a 2001-point grid ($0.1410 \to 0.4653$), so the
  crossing is unique — no second root hides inside the specced window. $J(x_c) = 150.15$ meV.
  *Passes as the kill.* **The tolerance was not widened; the number was measured and pinned.**
- **G3b — Why the window was wrong.** $x_c$ vs the threshold convention:
  $15\% \to 0.0524$, $18\% \to 0.2074$, $20\% \to 0.2970$, $25\% \to 0.4861$, $30\% \to 0.6395$ —
  about $0.039$ in $x$ per point of threshold. The spec's $(0.45, 0.65)$ is exactly the 25–30%
  band. v1.0's window was self-consistent with a *different, weaker* criterion than the 20% v1.0
  itself defined. Gated so that "just move the threshold" stays visible as the choice it is rather
  than an available tuning knob. *Passes.*
- **G4 — Singlet ground state everywhere (unchanged).** $S^2 < 10^{-6}$ at all 21 sweep points;
  measured max $\sim 10^{-30}$ (machine zero). Gated **non-vacuously**: the state the kick actually
  populates is verified to be a pure triplet, $S^2 = 2.000$, at energy exactly $J$ above the
  singlet — which doubles as an *independent* confirmation of $J$ from raw Hamiltonian expectation
  values rather than from ODMD. *Passes.*
- **G5 — Rider (added).** $J_{\text{Heis}} > J$ at every $x$: the perturbative form is a one-sided
  overestimate, never a wobble, because $4t^2/W$ is the leading term of
  $\sqrt{W^2/4 + 4t^2} - W/2$ whose remainder is strictly negative. This is what lets G3 treat
  $\Delta$ as a monotone validity measure rather than an absolute error that could cross zero.
  *Passes.*

## 6. Findings

1. **$x_c = 0.296953$, not $(0.45, 0.65)$.** The Heisenberg description of the interlayer exchange
   fails at 41% iodine substitution — much earlier in the series than the spec assumed. The error
   in v1.0 was not subtle and needed no simulation to catch: with $\Delta(0) = 14.1\%$ already at
   the pure bromide, a 20% threshold was always going to be crossed early.
2. **The whole deviation is one dimensionless group.** $\Delta$ is a function of
   $a = (2t/(U_0-U_s))^2$ alone, so a three-parameter alloy sweep collapses to a one-parameter
   curve and $x_c$ has a closed form. This is *why* the kill is clean: the gate checks a bisection
   against exact algebra, not against a finer grid.
3. **The boundary is soft, and that is the honest headline.** $x_c$ moves $0.05 \to 0.64$ as the
   threshold moves $15\% \to 30\%$. $\Delta(x)$ is smooth and strictly monotone — **nothing is
   non-analytic at $x_c$**. "Phase diagram" here means a validity boundary for an approximation, not
   a thermodynamic phase boundary, and $x_c$ is a property of the threshold convention at least as
   much as of the material. Any use of this number must quote the threshold with it.
4. **Moment collapse and Heisenberg failure are the same physics.** The local-moment fraction
   $\lVert S_z|\psi_0\rangle\rVert^2$ falls $0.890 \to 0.759$ monotonically while $\Delta$ rises
   $14.1\% \to 46.5\%$: charge fluctuations that eat the moment are what the perturbative expansion
   in $t/(U_0-U_s)$ is failing to capture. They are not independent diagnostics.
5. **$U_s$ cannot be held fixed.** v1.0's §3 interpolated only $t$ and $U$, but $U_s$ moves 24%
   across the series and enters $\Delta$ through $U_0 - U_s$. Freezing it at the bromide value
   shifts the whole curve; the gates interpolate all three.

## 7. Honest scope and caveats

- **VCA, per §2** — no disorder, no local configuration, no bowing, no relaxation, no re-screening.
  This is a model prediction, not an ab-initio alloy result, and the linearity of the parameter
  trends is an assumption rather than a finding.
- **Isolated dimer** — no in-plane kagome exchange, no band broadening. `nb3x8_gaps.py` shows that
  restoring coordination moves the *charge* gap of this same cluster by ~25%; no equivalent
  coordination correction has been computed for $J$, so the absolute exchange constants inherit an
  unquantified isolated-cluster error. $x_c$, depending only on the ratio $|2t/(U_0-U_s)|$, is
  likely more robust than $J$ itself — but that is an expectation, not a gated result.
- **Density-density interactions only**; the source paper reports non-density-density terms of a
  few meV.
- **Exact statevector.** No shot noise, no device noise, no Trotter error. The "ODMD" here is
  exercising the spectroscopy machinery, not a hardware claim.
- **The 20% threshold is a convention** (G3b), and $x_c$ is soft in it.
- **Two endpoints, one path, one family.** Nb₃Cl₈ and Nb₃F₈ alloys are not covered, and nothing
  here says the Br–I line is the relevant experimental one.

## 8. Reproducing

```bash
uv run python nb3x8_alloy.py            # writes data/nb3x8_alloy_diagram.csv, prints the table
make gates-nb3x8_alloy                  # G1, G2, G3, G3b, G4, G5  (~1 s)
```

`data/` is git-ignored; the CSV is regenerated by the command above.

## 9. v1.0, and the measurements that resolved it

> **G1 (v1.0):** *"`alloy_exchange_constant(0.0)` reproduces your verified Nb₃Br₈ exchange
> (J = 119.1 meV) and `alloy_exchange_constant(1.0)` reproduces Nb₃I₈ (J = 245.9 meV) to within
> 0.1%."*
> **Held.** Measured 119.1081 and 245.9196 meV — $6.8\times10^{-5}$ and $8.0\times10^{-5}$
> relative. The concern that the published 4-figure values were rounded too coarsely for a 0.1%
> tolerance was unfounded by three orders of magnitude.

> **G2 (v1.0):** *"Plotting hopping t(x) and Coulomb U(x) sweeps shows smooth, strictly monotonic
> linear behavior with zero discontinuities or singular matrix points."*
> **Held, with a missing statement.** All three parameters are strictly decreasing. But $|t|$ rises
> while $U_0 - U_s$ falls, which is the mechanism v1.0 never stated and the direct cause of the G3
> kill. And v1.0's §3 forgot $U_s$ entirely (Finding 5).

> **G3 (v1.0):** *"The critical concentration x_c (where Heisenberg error exceeds 20%) is resolved
> and located in the interval **x_c ∈ (0.45, 0.65)**, establishing a clear, falsifiable
> boundaries."*
> **FALSIFIED.** $x_c = 0.296953$, confirmed by closed form to $2.4\times10^{-15}$. The window was
> low by $0.153$ in $x$ and corresponds to a 25–30% threshold rather than the 20% v1.0 defined
> (G3b). The gate was rewritten to pin the measured value and to assert the old window is missed —
> not widened to accommodate it (`specs/README.md` step 5).

> **G4 (v1.0):** *"The ground state of H(x) across all x is verified to be a pure spin singlet
> (S² ≤ 10⁻⁶), ensuring the spin-response kick is structurally valid everywhere."*
> **Held**, at $\sim 10^{-30}$. Strengthened: v1.0's version would have passed against a
> broken all-zero $S^2$ operator, so the gate now also verifies $S^2 = 2$ on the kicked triplet.

> **§2 (v1.0):** *"the continuous ab-initio curve J(x) and the localized moment fraction across
> alloy concentrations are **wholly unmapped in the published literature**"* and *"By linearly
> interpolating the ab-initio cRPA hopping and Coulomb integrals … this tool maps the continuous
> phase diagram."*
> **Overclaim, withdrawn.** Interpolating ab-initio endpoints does not make the interpolant
> ab-initio; this is VCA (§2). The curve is a model prediction with a named approximation, and
> §2/§7 now say so up front rather than in a caveat appended after the numbers.
