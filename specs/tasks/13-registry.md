# Task Breakdown 13 — #18 Certificate Registry (v1 = your own ledger)
Goal: append-only ledger of certified results + `verify` command that recomputes and diffs. Depends: CertChem-M1 determinism.

1. **Entry format** — one JSON per entry: request-tuple hash, bracket, certificate, solver_version, timestamp, entry hash = SHA-256(canonical entry). Ledger = git repo of entries (append-only enforced by CI rule: PRs may only add files).
   ✓ Format doc; CI rejects a PR that modifies an existing entry. (S)
2. **`publish` command** — takes a `CertifiedResult` + original request → writes entry, computes hash, opens the commit.
   ✓ Golden-suite results published as entries 1–4. (S)
3. **`verify` command** — given an entry: re-run the calculation from the request tuple at the pinned solver_version, byte-diff result.
   ✓ All golden entries verify in CI on every solver release (catches accidental nondeterminism forever — the registry doubles as a determinism canary). (M)
4. **Version-pin strategy** — verify requires the entry's solver_version; document the policy: old versions installable from tags; entries never re-verified against newer physics.
   ✓ Verify against a deliberately wrong version fails with a clear message. (S)
5. **Citation format** — `certreg:ENTRYHASH` short-form + human page per entry (static site generated from the ledger).
   ✓ CertChem docs cite their own golden results by registry ID — dogfooding the citation loop. (S)
