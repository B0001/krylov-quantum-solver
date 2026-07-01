# Draft email to Prof. Malte Rösner (m.roesner@science.ru.nl)

> Review notes before sending:
> - Personalize the opening and signature — say who you are and why you're interested (a cold email
>   from a real, identifiable person is taken far more seriously than an anonymous one). I've left
>   `[…]` placeholders.
> - It's deliberately framed as a *check/question*, not a correction. The finding likely confirms
>   your-own-methodology reasoning, and the honest scope (isolated cluster, not the solid) is stated
>   up front. Keep it that way.
> - Attach `notes/nb3x8_cluster_gaps.md` (or paste the table).

---

**To:** m.roesner@science.ru.nl
**Subject:** Exact-diagonalization check of the Nb₃X₈ bilayer-cluster gaps (arXiv:2501.10320)

Dear Professor Rösner,

I recently read your paper on the Nb₃X₈ breathing-kagome family (arXiv:2501.10320) — a really nice
systematic downfolding across the halide series. [*One sentence on who you are / why you're
interested — e.g. background, what drew you to the paper.*]

As a small exercise I diagonalized the per-bilayer impurity cluster exactly — the two dimerized
trimer orbitals with on-site U₀, inter-layer hopping t_s⊥, and inter-site density-density U_s⊥ — from
your Table I/IV parameters, and compared the exact charge gap to Hubbard-I. Before reading anything
into it, I wanted to check the numbers with you. Two observations, both entirely within the isolated-
cluster picture (no self-consistent bath):

1. Hubbard-I reproduces the exact cluster gap to <2% for the strongly-correlated members (Nb₃F₈,
   HT-phase Cl/Br) but underestimates the weakly-correlated iodides by ~29% (Nb₃I₈ bulk) and ~12%
   (bilayer) — I assume the same reason you use CTHYB rather than Hubbard-I for the doped/harder cases.

2. Enlarging the correlated region to two dimers plus the weak inter-bilayer link (t_w⊥, U_w⊥) shifts
   the *exact* Nb₃I₈ gap by only ~5%, so the departure doesn't appear to be a two-site artifact; on
   that enlarged cluster the Hubbard-I error actually grows slightly rather than washing out.

I'm well aware this is the isolated cluster, not your cluster-DMFT — so it's a statement about the
impurity solver on the cluster, not about the material gap. My questions are simply: does this look
right to you, and is it already understood/expected? I've attached a one-page note with the full
table (all ten LT/HT bulk and bilayer parameter sets), the method, and the caveats.

Thank you for the paper, and for any thoughts.

Best regards,
[Your name]
[Affiliation / one-line context, if any]
