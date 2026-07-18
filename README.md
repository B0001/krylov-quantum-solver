# Krylov Quantum Solver

**A hybrid quantum–classical pipeline for molecular ground-state energies via real-time quantum Krylov subspace diagonalization.**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Qiskit](https://img.shields.io/badge/Qiskit-1.0+-6133BD.svg)](https://qiskit.org/)
[![Qiskit Nature](https://img.shields.io/badge/Qiskit_Nature-0.7+-6133BD.svg)](https://qiskit-community.github.io/qiskit-nature/)
[![PySCF](https://img.shields.io/badge/PySCF-2.0+-red.svg)](https://pyscf.org/)

## Overview

This project extracts molecular electronic-structure integrals with **PySCF**, maps them to a
qubit Hamiltonian with a **vetted Jordan–Wigner transform** (Qiskit Nature), and estimates the
ground-state energy with a **real-time quantum Krylov subspace** method. It is validated to
reproduce Full Configuration Interaction (FCI) to sub-milli-Hartree accuracy on small molecules.

> **History / honesty note.** An earlier version of this repository made much larger claims
> (bypassing the FCI wall for transition-metal superconductors, a "QCIVET" cryptographic noise
> guard, etc.). Those results were artifacts of a broken physics core: an incomplete Jordan–Wigner
> mapping, a reference state fixed to the empty vacuum, a "Krylov basis" of near-identical states,
> and asymmetric noise added directly to the Hamiltonian. The full audit and the rebuild are
> documented in [`REFACTOR_PLAN.md`](REFACTOR_PLAN.md). This README describes only what is now
> validated by the test suite.

On top of this pipeline sits the spec-gated **ODMD suite** — ground/excited energies from the
survival amplitude alone, certified two-sided brackets, Trotter/device-noise circuit paths,
single-signal error bars, and photoemission/optical/spin spectroscopy, applied end-to-end to the
Nb₃X₈ dimer clusters. The ladder, the material scorecard, and the recorded boundaries:
[`docs/ODMD_SUITE.md`](docs/ODMD_SUITE.md).

## Why certified bounds

This solver refuses to return an unaccompanied ground-state energy. Instead, it always emits two-sided brackets: a Rayleigh–Ritz upper bound (unconditional) and a Temple/Lehmann-type lower bound (explicit about its conditional gap input). For the full theoretical justification — why the problem itself is QMA-complete in our regime, why physical instance hardness remains open, and why this posture matters — see [`docs/theory/WHY_CERTIFIED_BOUNDS.md`](docs/theory/WHY_CERTIFIED_BOUNDS.md) and the arXiv-verified complexity proof graph in [`docs/theory/complexity-proof-graph-v2.cypher`](docs/theory/complexity-proof-graph-v2.cypher).

## Method

1. **Classical pre-processing** (`hybrid_quantum_solver/chemistry_gateway.py`) — PySCF HF/CASCI
   gives active-space one- and two-body integrals.
2. **Qubit Hamiltonian** (`hybrid_quantum_solver/molecular_hamiltonian.py`) — Qiskit Nature's
   `JordanWignerMapper` builds a `SparsePauliOp`; the Hartree–Fock state is the reference; the
   nuclear/core constant is tracked as an energy offset.
3. **Real-time quantum Krylov** (`hybrid_quantum_solver/quantum_krylov_solver.py`) — builds the
   Krylov space |φₖ⟩ = e^(−ikΔtH)|φ_HF⟩, forms the Hermitian subspace matrices H and S, and solves
   the generalized eigenproblem `H c = E S c` by thresholded canonical orthogonalization. Because
   every basis vector is a genuine state, the estimate is **variationally bounded below by the true
   ground state**.
4. **Noise** (`hybrid_quantum_solver/noise.py`) — optional finite-shot sampling noise (added
   Hermitian-symmetrically) and a real Qiskit Aer `NoiseModel` builder, with a noise-aware overlap
   cutoff for stability.

The end-to-end pipeline is `hybrid_quantum_solver/pipeline.py` (`run_geometry`, `run_from_integrals`).

## Validation

Reproduced against an **independent PySCF FCI** solve (`benchmark_krylov.py`, Krylov dimension M = 10):

| System | Qubits | Hartree–Fock error | Quantum Krylov error | ≥ FCI (variational) |
| :--- | :--- | :--- | :--- | :--- |
| H₂ (0.74 Å, equilibrium) | 4 | 20.5 mHa | < 1×10⁻⁶ mHa | ✅ |
| H₂ (2.0 Å, stretched / multireference) | 4 | 164.8 mHa | < 1×10⁻⁶ mHa | ✅ |
| H₄ chain (1.0 Å) | 8 | 67.8 mHa | 0.0002 mHa | ✅ |
| LiH (1.6 Å) | 12 | 20.5 mHa | 0.14 mHa | ✅ |

*(1 kcal/mol "chemical accuracy" ≈ 1.6 mHa.)* The stretched-H₂ row is the informative one: Hartree–Fock
is wrong by 165 mHa, and the Krylov solver still recovers FCI. Energies decrease toward FCI as M grows
and never drop below it. Under shot noise the error grows smoothly as ~1/√shots and stays bounded —
in contrast to the previous code, which returned values hundreds of Hartree below the true minimum.

### N₂ dissociation (the multireference rung)

`benchmark_n2.py` breaks the N₂ triple bond in a CAS(6,6) active space (12 qubits). Hartree–Fock
collapses as the bond stretches; quantum Krylov tracks the exact CASCI curve throughout (exact
reference cross-checked against an independent PySCF CASCI solve, agreeing to 1×10⁻⁶ Ha):

| Bond length | Hartree–Fock error | Quantum Krylov error (M = 12) |
| :--- | :--- | :--- |
| 1.0 Å | 100.9 mHa | 0.0006 mHa |
| 1.1 Å (≈ equilibrium) | 126.6 mHa | 0.008 mHa |
| 1.3 Å | 192.9 mHa | 0.042 mHa |
| 1.6 Å | 328.5 mHa | 0.649 mHa |
| 2.1 Å (stretched) | 625.1 mHa | 0.357 mHa |

As the molecule becomes strongly multireference the Krylov space needs a larger dimension M to
converge (honest behaviour — Hartree–Fock has poor overlap with the true ground state there), but it
still reaches sub-milli-Hartree accuracy where single-reference HF is off by half a Hartree.

For active spaces beyond exact FCI's reach, `hybrid_quantum_solver/dmrg_reference.py` provides a
`reference_energy(method="auto")` that uses **DMRG (block2)** when installed and falls back to exact
FCI otherwise (the shared integral convention is validated by the test suite).

### Resource accounting

`benchmark_resources.py` reports the honest circuit cost of the on-hardware solver — replacing the
previous code's fictitious wall-clock figures (e.g. a 19-hour 16-qubit job). The numbers are sobering
and transparent:

| System | Trotter step (CX) | Deepest M=6 circuit (CX) | Circuits @ M=6 | Shots @ 8192 |
| :--- | :--- | :--- | :--- | :--- |
| H₂ (5 qubits) | 70 | 2,850 | 21 | 688k |
| N₂ CAS(6,6) (13 qubits) | 6,534 | 238,210 | 21 | 688k |

N₂'s ~6.5k two-qubit gates **per Trotter step** is already far beyond what current hardware can run
with useful fidelity — which is exactly why the next research lever is shallower time evolution
(higher-order/longer-step Trotter, qubitization, circuit optimization) rather than more qubits.

## Installation

```bash
conda create -n chem python=3.11 && conda activate chem
pip install -r requirements.txt
```

## Quickstart

From a molecular geometry:

```python
from hybrid_quantum_solver.pipeline import run_geometry

result = run_geometry(atom="Li 0 0 0; H 0 0 1.6", basis="sto3g", krylov_dim=10)
print(result.summary())   # E=-7.882... Ha (HF=-7.861..., ref=-7.882..., ...)
```

From a CIF via CASCI active space (materials path):

```bash
python run_hybrid_solver.py --input_file data/nb_structures/NbN_mp-2634.cif \
    --active_space 8,8 --krylov_dim 8           # add --shots 8192 to model sampling noise
```

> **Scientific caveat (materials).** The CIF path builds a *finite molecular cluster* from the
> unit-cell atoms with no periodic boundary conditions — it is **not** a calculation of the periodic
> solid. Transition-metal systems are a research target, to be benchmarked against DMRG/AFQMC where
> tractable, not a validated result. See `REFACTOR_PLAN.md`, Phase 4.

## Tests

```bash
pytest tests/ -v
```

Three gates (18 tests): `test_reference_energies.py` (qubit Hamiltonian == FCI; HF == RHF; active
space == CASCI), `test_krylov_convergence.py` (converges to FCI, respects the variational floor),
and `test_noise_resilience.py` (shot noise is bounded and improves with shots).

## Repository layout

```
hybrid_quantum_solver/
  chemistry_gateway.py        # PySCF integral extraction (CIF / geometry -> CASCI)
  molecular_hamiltonian.py    # vetted Jordan-Wigner qubit Hamiltonian + HF reference
  quantum_krylov_solver.py    # real-time quantum Krylov subspace diagonalization (exact evolution)
  trotter_krylov.py           # Trotter-circuit Krylov + qiskit-aer device-noise expectation path
  hardware_krylov.py          # on-hardware Krylov: Hij/Sij measured via ancilla Hadamard tests
  noise.py                    # shot-noise model, Hermitisation, Aer NoiseModel, ZNE (fold + extrapolate)
  pipeline.py                 # end-to-end: run_geometry / run_from_integrals
  dmrg_reference.py           # classical reference: exact FCI, or DMRG (block2) for larger spaces
  orchestrate_hybrid_pipeline.py, quantum_sampler.py   # LEGACY/BROKEN, regression fixtures only
benchmark_krylov.py           # honest FCI benchmark table (H2 / H4 / LiH)
benchmark_n2.py               # N2 dissociation curve vs exact CASCI (multireference rung)
benchmark_resources.py        # circuit-level resource accounting (depth / CX / shot budget)
REFACTOR_PLAN.md              # full scientific audit + phased rebuild
tests/                        # validation gates
```

## References

- Parrish & McMahon, *Quantum filter diagonalization* (2019).
- Stair, Huang & Evangelista, *A multireference quantum Krylov algorithm*, JCTC 16, 2236 (2020).
- Klymko et al., *Real-time evolution for ultracompact Hamiltonian eigenstates*, PRX Quantum 3, 020323 (2022).
- Epperly, Lin & Nakatsukasa, *A theory of quantum subspace diagonalization*, SIAM J. Matrix Anal. Appl. 43, 1263 (2022).
