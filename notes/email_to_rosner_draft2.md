# Draft email to Prof. Malte Rösner — v2 (DRAFT, review before sending)

**Status: DRAFT.** Supersedes the withdrawn v1 (`email_to_rosner_draft.md`), whose premise — a
~29 % Hubbard-I discrepancy — dissolved once band broadening was restored. This version has the
opposite premise: the exact-cluster check *supports* the paper's method, and along the way fixes
several constants of their own downfolded models that the paper did not report (optical gaps,
exciton bindings, inter-layer exchange couplings — two of them in closed form). That is a much
more natural unsolicited note: no challenge, just reference numbers + two small closed formulas.

Attach or link: `notes/nb3x8_cluster_gaps.md` (the note now contains all tables).

---

**To:** Malte Rösner
**Subject:** Exact reference constants for the Nb₃X₈ dimer clusters of arXiv:2501.10320

Dear Prof. Rösner,

I enjoyed your group's paper on the breathing-mode kagome Nb₃X₈ family (arXiv:2501.10320). Since
the per-bilayer impurity problem you downfold to — two inter-layer-dimerized trimer orbitals with
U₀, t_s⊥ and U_s⊥ — is only four spin-orbitals, I solved it exactly for all ten parameter sets in
your Tables I and IV, as an independent cross-check. Two results may be of interest; a short note
with all tables and the reproduction code is attached.

1. **Your Hubbard-I / cluster-DMFT gaps look robust.** On the *isolated* cluster, exact
diagonalization sits up to ~29 % above Hubbard-I for the iodides — but that offset is an artifact
of the missing band broadening: restoring inter-dimer coordination (a coordination scan plus a
DMRG chain toward the 1-D thermodynamic limit) collapses the exact Nb₃I₈ gap from 842 to
~600–650 meV, right at your cluster-DMFT value (~599). So the exact check ends up *confirming*
your method for the solid gaps; the isolated-cluster gaps themselves (842–3984 meV across the
family) may still be useful as reference values for the downfolded models.

2. **A few constants of your models the paper did not report, two in closed form.** The dimer's
only dipole-active singlet is the ionic-odd state at exactly U₀, so the cluster optical gap is
Δ_opt = U₀ − E₀ with E₀ = (U₀+U_s⊥)/2 − √(((U₀−U_s⊥)/2)² + 4t_s⊥²), and with the charge gaps this
gives eV-scale exciton bindings for the bilayers (1.1–2.0 eV) that collapse to ~68 meV for bulk
Nb₃I₈. Likewise the inter-layer exchange J = √(((U₀−U_s⊥)/2)² + 4t_s⊥²) − (U₀−U_s⊥)/2 comes out
at 246 meV (bulk I) down to ≈0 (F) — with the perturbative 4t²/(U₀−U_s⊥) overestimating the
iodides by 40–47 %, i.e. the iodide dimer is well beyond the Heisenberg regime, and its local
moment is reduced by ~24 % by charge fluctuations. All of this is isolated-cluster physics
(density–density terms only, no in-plane exchange), intended as anchors for the downfolded
Hamiltonians rather than solid-state predictions.

Everything regenerates from a small test-gated repository (exact diagonalization cross-checked
against DMRG and, where applicable, analytic limits); I would of course be happy to share it, and
to be corrected if any of this duplicates work I have missed.

With best regards,
Benjamin Hess

---

*Checklist before sending: (i) re-read arXiv:2501.10320 v-latest to confirm none of the Section-3
constants appear in it or its SM; (ii) decide whether to attach the note as PDF or link a repo;
(iii) confirm affiliation/signature line.*
