# CertChem — certified chemistry, solver-as-library

CertChem-M1: a small, dependency-light contract over this repo's validated quantum-Krylov
solver. It answers one question honestly — *what is the ground-state energy, and how sure are
we?* — by returning a **certified two-sided bracket** that provably encloses the exact
active-space ground energy, not just a point estimate you have to take on faith.

Foundation for the portfolio (SenseForge, CertLabel, the registry, screening loops, …). See
`specs/tasks/01-certchem-m1.md` and `architecture/interfaces/solver-library-contract.md`.

## The contract

```python
from certchem import certified_energy, Mode

# CERTIFIED (default): floor-checked bracket + certificate
result = certified_energy("H 0 0 0; H 0 0 0.735", "sto-3g", cas=(2, 2))
print(result.bracket.lower_hartree, result.bracket.upper_hartree)  # encloses FCI
print(result.certificate.method)                                   # temple_bound + variational_floor

# FAST: a bare float, no guarantee — the return TYPE makes the two impossible to confuse
e = certified_energy("H 0 0 0; H 0 0 0.735", "sto-3g", cas=(2, 2), mode=Mode.FAST)
```

Guarantees a caller may rely on (architecture ADR-0001/0004):

1. Holding a `CertifiedResult` means the variational **floor check passed** — a sub-floor
   estimate raises `FloorViolationError`, it is never returned as a number.
2. `bracket.lower_hartree <= best_estimate_hartree <= upper_hartree`, always (enforced at
   construction).
3. Same inputs → byte-identical serialized result (energies quantized to 1e-12 Ha; basis of
   content-hash caching).
4. `Mode.FAST` never yields a `Bracket`.

Outside the validated envelope (`MAX_SPIN_ORBITALS = 16`, validated bases) `check_caps` raises
`CapExceededError` naming the violated cap — before any solve runs.

## Quickstart — watch Hartree–Fock fail

```
python -m certchem.examples.n2_quickstart
```

For N₂ CAS(6,6), Hartree–Fock misses ~0.13 Ha of correlation energy (>80× chemical accuracy).
CertChem returns a certified bracket ~0.04 mHa wide that encloses the exact answer — and HF
sits provably *above* it. Source: `certchem/examples/n2_quickstart.py`.

## How it works (no new physics)

`certified_energy` is plumbing over validated primitives: PySCF integrals →
`build_molecular_hamiltonian` → real-time quantum Krylov → the Temple/Weinstein certified
bracket (`temple_bounds.krylov_bracket`, Invention #20). The **only** path from an energy
estimate to a result passes through `floor_guard()` (`certchem/core.py`).

## Honest scope

Certification is **sector-restricted** to the lowest reachable levels — the same scope as QKSD
itself (see `specs/SPEC_temple_bracket.md`). The self-mode Temple premise is gated at Krylov
dim ≥ 6. A vacuous (infinite) lower bound raises `ConvergenceError` rather than quoting an
`inf`-wide bracket. Statevector path only; hardware shot cost of ⟨H²⟩ is not modeled here.
