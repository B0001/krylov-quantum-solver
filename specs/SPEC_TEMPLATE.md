# SPEC: <one-line title>

**Status:** DRAFT — awaiting approval. No implementation until the acceptance gates are agreed.

> Copy this file to `specs/SPEC_<slug>.md` and fill it in. A spec is a *falsifiable hypothesis*,
> not a contract: if implementation shows a gate is wrong, **change the gate and record why**
> (that mismatch is the finding). Keep it short — the gates carry the weight, not the prose.

---

## 1. Goal

One paragraph: the single claim this spec makes, stated so it could be proven false.

## 2. Background and honest framing

- Why this is worth doing; prior art / references.
- **What you can claim** if the gates pass.
- **What you cannot claim** (state up front, not after): scope limits, "reproduction not novelty",
  "no quantum advantage at this scale", finite-cluster/minimal-basis, etc. Be specific.

## 3. Approach

The method, in enough detail that the gates below are well-defined. Name the ground-truth
**reference** every result will be checked against (FCI, DMRG, experiment, an analytic limit).
If there is no checkable reference, this is not yet a spec — go find one.

## 4. Public interface

Exact functions/modules/CLIs to add or change. Prefer **composing validated primitives** over new
code. List signatures and return types.

```
<module>.<function>(args) -> <result type>
<driver/CLI>             -> <artifact>
```

## 5. Acceptance criteria (validation gates)

Each gate is an automated check in `tests/test_<slug>_spec.py` (test-first), and each must be
**falsifiable and cheap to run** (a small/fast case that would fail loudly if the claim is wrong).

- **G1 — <name>.** <concrete, numeric pass condition, e.g. `|E − E_ref| < 1e-4 Ha`>.
- **G2 — <name>.** <…>
- **G3 — <name>.** <…>
- (Stretch goals that are *not* pass/fail gates go here, labelled as such.)

> State which gate is the definition of "done". If a gate proves unsatisfiable during
> implementation, revise it here with a short note explaining what reality showed (see
> `specs/SPEC_hchain_tdl.md` G1 for a worked example).

## 6. Implementation plan (test-first)

1. Write `tests/test_<slug>_spec.py` encoding the gates (initially failing).
2. Implement the minimum code to pass them (reuse > write).
3. Iterate to green. Run block2/DMRG gates in their own process (`make gates` handles this).

## 7. Out of scope

Bullet list of what this explicitly does **not** do (a follow-up spec can).

## 8. Caveats and risks

- **R1:** <main technical risk + mitigation/fallback>.
- Honest limitations of the result (so they are never presented as more than they are).

## 9. Deliverables

- `<module>` — <new/changed code>.
- `<driver>` — <artifact it produces>.
- `tests/test_<slug>_spec.py` — gates.
- Results summary (with the §2/§7 caveats) in the PR description.
