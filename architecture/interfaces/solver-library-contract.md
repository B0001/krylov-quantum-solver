# Interface: Solver Core public Python contract

The in-process boundary consumed by CertChem's workers and SenseForge (ADR-0005).
Semantic-versioned; breaking changes bump major and invalidate the result cache
(cache key includes solver_version, ADR-0008).

```python
from dataclasses import dataclass
from enum import Enum

class Mode(Enum):
    FAST = "fast"            # point estimate, no guarantee
    CERTIFIED = "certified"  # bracket + certificate, floor-checked

@dataclass(frozen=True)
class Bracket:
    lower_hartree: float
    upper_hartree: float
    best_estimate_hartree: float
    @property
    def width(self) -> float: ...

@dataclass(frozen=True)
class Certificate:
    method: str                # e.g. "temple_bound + variational_floor"
    floor_check: str           # always "pass" if you are holding this object
    krylov_dim: int
    convergence: str           # "converged" | "converged_marginal"
    solver_version: str
    manifest: dict | None      # ShallowForge provenance, if compiled path used

@dataclass(frozen=True)
class CertifiedResult:
    bracket: Bracket
    certificate: Certificate

def certified_energy(molecule, basis: str, cas: tuple[int, int],
                     mode: Mode = Mode.CERTIFIED) -> CertifiedResult | float: ...
def certified_gap(molecule, basis: str, cas: tuple[int, int],
                  gap_type: str, mode: Mode = Mode.CERTIFIED) -> CertifiedResult | float: ...
def certified_reaction(species: list[tuple[float, "Molecule", tuple[int, int]]],
                       basis: str) -> CertifiedResult: ...
```

## Exceptions (the contract's failure vocabulary)
| Exception | Raised when | Guarantee |
|---|---|---|
| `FloorViolationError` | estimate falls below the variational floor | never returns a number; carries diagnostics |
| `CapExceededError` | system outside validated envelope | names the violated cap |
| `ConvergenceError` | ODMD signal insufficient at configured Krylov dim | partial data attached, no estimate |

## Invariants callers may rely on
1. A returned `CertifiedResult` implies the floor check passed (ADR-0001).
2. `bracket.lower <= best_estimate <= bracket.upper`, always.
3. Same inputs + same `solver_version` → byte-identical result (determinism; basis
   of content-hash caching).
4. `Mode.FAST` never returns a `Bracket` — the types make blurring impossible (ADR-0004).
