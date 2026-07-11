# Task Breakdown 1 — CertChem M1: Solver-as-Library
Goal: `certified_energy(mol, basis, cas, mode) -> CertifiedResult` per the library contract. Foundation for #1, #5, #10, #12, #15, #17, #18, #19, SenseForge.

1. **Freeze the contract types** — copy `Bracket`, `Certificate`, `CertifiedResult`, `Mode`, and the three exceptions from `interfaces/solver-library-contract.md` into `src/certchem/contract.py`, frozen dataclasses, full type hints.
   ✓ mypy clean; importable with zero solver deps. (S)
2. **Limits module** — `limits.py`: `MAX_SPIN_ORBITALS=16`, allowed bases, `check_caps(mol, basis, cas)` raising `CapExceededError` naming the violated cap.
   ✓ Unit tests: one pass case, one fail per cap, error message names the cap. (S)
3. **Wrap the existing pipeline** — `core.py`: route molecule → PySCF Hamiltonian → ODMD estimate → floor guard → Temple bounds → assemble `CertifiedResult`. No new physics; only plumbing existing repo functions behind the contract.
   ✓ H₂/STO-3G returns a bracket containing the known FCI value. (M)
4. **Floor guard as chokepoint** — ensure the ONLY path from estimator to bounds passes `floor_guard()`; grep-able assertion + a test that monkeypatches a bad estimate and asserts `FloorViolationError`, never a return.
   ✓ Injected −999 Ha estimate raises; no `CertifiedResult` constructed. (S)
5. **Determinism audit** — same inputs twice → byte-identical serialized results. Pin every RNG seed; expose `solver_version()` from package metadata.
   ✓ Round-trip hash equality test in CI. (M)
6. **Golden regression gate** — pytest module running H₂, H₄, LiH, N₂-CAS(6,6); asserts (a) bracket contains FCI reference, (b) |estimate−FCI| ≤ documented tol, (c) floor pass.
   ✓ All four green; wired as required CI check. One containment failure = red build. (M)
7. **`Mode.FAST` path** — returns bare float, no bracket ever constructible from fast path (type-level).
   ✓ Test asserts return type is `float`, not `CertifiedResult`. (S)
8. **README + one worked example** — 30-line quickstart computing N₂ with the half-Hartree HF-failure comparison shown.
   ✓ Example runs top-to-bottom in CI. (S)

Definition of done: `pip install -e . && pytest` green; tag `v0.1.0`. Everything downstream imports this tag.
