# SPEC: Certified HF subspace overlap on molecules (SPEC-21b integration)

**Status:** Specced 2026-07-26. The molecular reachable-cluster demonstration
`SPEC_certified_overlap_subspace.md` named as its out-of-scope follow-up.
**Depends on:** `certified_overlap` (SPEC-21b block form: `certify_subspace_overlap`,
`ClusterGapCertificate`), `certified_gaps`/`QuantumKrylovSolver` (the self-mode Krylov floor),
`MolecularHamiltonian.hf_state`. Sibling of `SPEC_hf_overlap_certificate.md` (the d = 1 case).

## Claim

On a strongly-multireference molecule the Hartree–Fock state is a **poor proxy for the reachable
ground *state* but a strong proxy for the reachable ground *eigenspace***. SPEC-21b's block
certificate, fed the repo's own premise-gated self-mode Krylov E_d floor (`theta_d − sigma_d`,
the Weinstein floor generalized from d = 1), certifies `γ_min ≤ ‖P_S u‖` for the lowest-d
reachable eigenspace — **non-vacuously where the SPEC-21 single-vector certificate is vacuous**,
and with **no oracle**.

The worked system is **square H₄** (four H on a square of side a; STO-3G): as the square tightens
the HF state spreads across the two lowest reachable levels, and its overlap with the individual
reachable ground state falls below the residual-to-gap ratio — collapsing the d = 1 certificate —
while its overlap with the two-level eigenspace stays large.

## Sector honesty

Reachable-sector restricted, exactly like QKSD / temple_bounds / certified_gaps / the d = 1
HF-overlap path. P_S spans the lowest d **reachable** eigenstates (nonzero HF amplitude); β floors
the (d+1)-th **reachable** level. The block sin-θ bound carries over with "spectrum" read as
"reachable spectrum" because HF's amplitude on every unreachable level is zero — those levels
contribute neither to r² nor to sin²θ, so β may ignore them even if they sit below it.

## Gates (square H₄, sides a ∈ {1.0, 1.2, 1.4} Å unless noted; d = 2)

- G1 (validity, killable, zero-tol): `γ_min ≤ ‖P_S u‖` (dense reachable reference) for self mode
  at M ∈ {6, 8} and oracle mode, across the sweep. One escape kills the block bound on molecules.
- G2 (**the finding** — molecular vacuous-vs-useful): across the sweep the **d = 1** certificate
  (`certify_hf_overlap`) is **VACUOUS**, while the **d = 2** block certificate is **non-vacuous**
  with γ_min ≥ 0.45. The subspace view rescues a certificate the single-vector view throws away,
  on a real molecule. (Observed self-mode d=2 γ: a=1.4 → 0.50, a=1.2 → 0.80, a=1.0 → 0.92; d=1 ≡ 0.)
- G3 (no oracle needed): the self-mode γ_min matches the oracle-mode γ_min to within 0.05 at
  M ≥ 8 — the M ≥ 6 Krylov space resolves the size-2 reachable cluster, so the certificate needs
  no exact E_d input. (Self converges up in M toward oracle: M=6 0.776 → M≥12 0.809 at a=1.2.)
- G4 (premise boundary): self mode with M < 6 raises — the temple_bracket M ≥ 6 boundary inherited
  loudly, not silently ignored; oracle mode has no such premise and does not raise.
- G5 (usefulness trend): the d = 2 self-mode floor is non-decreasing as the square tightens
  (a: 1.4 → 1.2 → 1.0 ⇒ γ_min 0.50 → 0.80 → 0.92) — stronger multireference character puts *more*
  of the HF weight into the two-level cluster, and the certificate measures it.

## Out of scope

Automatic cluster-size detection (d is a caller input; picking d from certified spectral
structure is its own hypothesis). d > 2 molecular clusters (the machinery is general — G1 covers
general d in the synthetic SPEC-21b gate — but the demonstration fixes d = 2). Shot noise on
r/β. Systems beyond dense-`eigh` validation reach. Non-square multireference families (stretched
linear H₄, Nb₃X₈ singlet–triplet) — clean extensions once this lands.

## Honest caveats

Exact-statevector only (⟨H²⟩ shot cost unmodeled, inherited from temple_bounds). The self-mode
E_d floor requires the Krylov space to **resolve** the (d+1)-th reachable level; a level with
vanishing HF amplitude near the cluster boundary is the failure mode (the certificate would floor
the wrong level). The M ≥ 6 inheritance plus the observed self-vs-oracle agreement (G3) is the
guard here; a general certified resolvability test for d > 1 is deferred. Picking d is the
caller's physics judgment, not certified.
