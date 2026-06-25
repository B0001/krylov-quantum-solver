#!/usr/bin/env python3
"""
Phase 2 validation gate (see REFACTOR_PLAN.md).

The real-time quantum Krylov solver (hybrid_quantum_solver.quantum_krylov_solver) must,
on H2 and LiH:

  1. respect the VARIATIONAL FLOOR -- no estimate ever dips below FCI
     (the original pipeline returned energies hundreds of Ha *below* the true minimum);
  2. CONVERGE to FCI (best estimate within 1 mHa) as the Krylov dimension M grows;
  3. capture real CORRELATION beyond Hartree-Fock;
  4. decrease MONOTONICALLY on every step where the effective Krylov rank actually grows
     (exact step-monotonicity is not asserted: once the Krylov space saturates, canonical
     orthogonalisation introduces ~1e-5 Ha conditioning ripples -- these are above the
     floor and do not affect convergence).

FCI is taken from an INDEPENDENT PySCF solve.

Run:  pytest tests/test_krylov_convergence.py -v
"""
import pytest

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

# name : (geometry, basis, max Krylov dimension)
CASES = {
    "H2":  ("H 0 0 0; H 0 0 0.74", "sto3g", 6),
    "LiH": ("Li 0 0 0; H 0 0 1.6", "sto3g", 12),
}


def _pyscf_fci_total(atom, basis):
    from pyscf import gto, scf, fci
    mol = gto.M(atom=atom, basis=basis, verbose=0)
    mf = scf.RHF(mol)
    mf.kernel()
    e_fci, _ = fci.FCI(mf).kernel()
    return float(e_fci)


@pytest.mark.parametrize("name", list(CASES))
def test_krylov_converges_to_fci(name):
    atom, basis, max_dim = CASES[name]
    mh = build_molecular_hamiltonian(atom=atom, basis=basis)
    fci = _pyscf_fci_total(atom, basis)

    steps = QuantumKrylovSolver(mh).convergence(max_dim)
    energies = [s.energy for s in steps]

    # 1. variational floor -- never below the true ground state
    assert min(energies) >= fci - 1e-7, (
        f"{name}: variational violation, min={min(energies):.8f} < FCI={fci:.8f}")

    # 2. converges to FCI within 1 mHa
    best_err = min(abs(e - fci) for e in energies)
    assert best_err < 1e-3, f"{name}: best error {best_err*1e3:.3f} mHa exceeds 1 mHa"

    # 3. captures real correlation beyond Hartree-Fock (M=1 is the HF reference)
    hf_err = abs(energies[0] - fci)
    assert hf_err - best_err > 5e-3, (
        f"{name}: no meaningful improvement over Hartree-Fock "
        f"(HF err {hf_err*1e3:.3f} mHa, best err {best_err*1e3:.3f} mHa)")

    # 4. monotone decrease on every step where the Krylov rank grows
    for a, b in zip(steps, steps[1:]):
        if b.rank > a.rank:
            assert b.energy <= a.energy + 1e-9, (
                f"{name}: energy rose on rank increase "
                f"M{a.dim}(r{a.rank})={a.energy:.8f} -> M{b.dim}(r{b.rank})={b.energy:.8f}")


def test_first_step_is_hartree_fock():
    """M=1 (the bare reference, no evolution) must equal the Hartree-Fock energy."""
    mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")
    step1 = QuantumKrylovSolver(mh).solve(1)
    assert abs(step1.energy - mh.hf_energy) < 1e-9


def test_energy_offset_frame():
    """The returned energy is a physical total (electronic eigenvalue + offset)."""
    mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")
    step = QuantumKrylovSolver(mh).solve(4)
    # H2/STO-3G total energy is ~ -1.137 Ha, nowhere near the offset-free electronic frame.
    assert -1.20 < step.energy < -1.10
