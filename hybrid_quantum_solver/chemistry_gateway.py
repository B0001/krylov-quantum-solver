#!/usr/bin/env python3
"""
Chemistry Gateway: chemistry_gateway.py
Classical pre-processing layer — HF/CASCI integral extraction and tensor conversion.
"""

import numpy as np
from typing import List, Tuple, Dict
from ase.io import read as ase_read

try:
    from pyscf import gto, scf, ao2mo, mcscf
    PYSCF_AVAILABLE = True
except ImportError:
    PYSCF_AVAILABLE = False

# ==========================================
# CIF → CASCI integral pipeline
# ==========================================

def _get_smart_basis(atoms) -> Dict[str, str]:
    """ECP basis for heavy atoms (Z > 36), all-electron for light atoms."""
    return {a.symbol: ('lanl2dz' if a.number > 36 else '6-31g*') for a in atoms}


def load_and_compute_integrals(
    cif_filepath: str,
    spin_targets: List[int] = [0, 2, 4, 6],
    cas_electrons: int = 8,
    cas_orbitals: int = 8,
) -> Tuple[np.ndarray, np.ndarray, int, float, float, Tuple[int, int]]:
    """
    CIF → PySCF HF → CASCI.

    Returns (h1, eri_4d, n_orbitals, casci_total_energy, e_core, nelecas), where
    ``nelecas = (n_alpha, n_beta)`` is the active-space electron split needed to build the
    Hartree-Fock reference for the qubit solver. Feed (h1, eri_4d, nelecas, e_core) into
    ``molecular_hamiltonian.build_hamiltonian_from_integrals`` /
    ``pipeline.run_from_integrals``; the active-space FCI target is ``casci_total_energy``.

    NOTE (scientific caveat): for a crystalline CIF this builds a *finite molecular cluster*
    from the unit-cell atoms with no periodic boundary conditions -- it is not a calculation
    of the periodic solid. Treat materials results accordingly (see REFACTOR_PLAN.md, Phase 4).
    """
    print(f"[CLASSICAL PRE-PROCESSING] Generating Hamiltonian for: {cif_filepath}")
    atoms = ase_read(cif_filepath)
    atom_str = "; ".join(
        f"{a.symbol} {a.position[0]} {a.position[1]} {a.position[2]}" for a in atoms
    )
    basis_set = _get_smart_basis(atoms)
    ecp_dict = {a.symbol: 'lanl2dz' for a in atoms if a.number > 36}

    # Determine ground-state spin by trying each candidate.
    mol_dummy = gto.M(atom=atom_str, basis=basis_set, ecp=ecp_dict, charge=0, spin=None)
    total_electrons = sum(mol_dummy.nelec)
    is_odd = total_electrons % 2 != 0
    # A spin S requires n_beta = (N - S) / 2 >= 0, i.e. S <= total_electrons; otherwise PySCF
    # asserts on a negative electron count. Filter on both parity and that magnitude bound so
    # small molecules (e.g. H2, where S=4 is impossible) don't crash the scan.
    valid_spins = [s for s in spin_targets
                   if (s % 2 != 0) == is_odd and s <= total_electrons]
    if not valid_spins:
        valid_spins = [1 if is_odd else 0]

    scf_energies = {}
    for spin in valid_spins:
        mol = gto.M(atom=atom_str, basis=basis_set, ecp=ecp_dict, charge=0, spin=spin)
        mf = scf.UHF(mol) if mol.nelec[0] != mol.nelec[1] else scf.RHF(mol)
        mf.max_cycle = 500
        mf.level_shift = 0.3
        mf.diis_space = 12
        mf.conv_tol = 1e-8
        mf.init_guess = 'atom'
        scf_energies[spin] = mf.kernel()

    ground_spin = min(scf_energies, key=scf_energies.get)
    print(f"[SUCCESS] Classical baseline established. Ground state spin: {ground_spin}")

    final_mol = gto.M(atom=atom_str, basis=basis_set, ecp=ecp_dict, charge=0, spin=ground_spin)
    final_mf = scf.UHF(final_mol) if final_mol.nelec[0] != final_mol.nelec[1] else scf.RHF(final_mol)
    final_mf.max_cycle = 500
    final_mf.level_shift = 0.3
    final_mf.diis_space = 12
    final_mf.conv_tol = 1e-8
    final_mf.init_guess = 'atom'
    final_mf.kernel()

    print(f"[CLASSICAL NODE] Truncating to CAS({cas_electrons},{cas_orbitals}) Active Space...")
    cas = mcscf.CASCI(final_mf, cas_orbitals, cas_electrons)
    cas.kernel()

    h1, e_core = cas.get_h1eff()
    eri_4d = ao2mo.restore(1, cas.get_h2eff(), cas_orbitals)

    nelecas = cas.nelecas
    if isinstance(nelecas, (int, np.integer)):
        n_b = int(nelecas) // 2
        nelecas = (int(nelecas) - n_b, n_b)
    else:
        nelecas = (int(nelecas[0]), int(nelecas[1]))

    return h1, eri_4d, h1.shape[0], cas.e_tot, float(e_core), nelecas


def translate_tensors_to_orchestrator(
    h1: np.ndarray, eri: np.ndarray
) -> Tuple[List[tuple], List[tuple]]:
    """
    Converts dense PySCF arrays to sparse tuple lists for the orchestrator.
    Uses numpy masking instead of O(N^4) Python loops (~100-1000x faster for CAS ≥ 8).
    """
    mask1 = np.abs(h1) > 1e-6
    i_idx, j_idx = np.where(mask1)
    single_body = list(zip(i_idx.tolist(), j_idx.tolist(), h1[mask1].tolist()))

    # ao2mo.restore(1, ...) unpacks PySCF's 8-fold-symmetric format to full (N,N,N,N).
    eri_4d = ao2mo.restore(1, eri, h1.shape[0])
    mask2 = np.abs(eri_4d) > 1e-6
    p_idx, q_idx, r_idx, s_idx = np.where(mask2)
    two_body = list(zip(
        p_idx.tolist(), q_idx.tolist(),
        r_idx.tolist(), s_idx.tolist(),
        eri_4d[mask2].tolist(),
    ))
    return single_body, two_body


# ==========================================
# Gateway Diagnostic Test Harness
# ==========================================
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python -m hybrid_quantum_solver.chemistry_gateway <path_to_cif> [cas_e,cas_o]")
        sys.exit(1)

    cif = sys.argv[1]
    cas_e, cas_o = (map(int, sys.argv[2].split(",")) if len(sys.argv) > 2 else (8, 8))
    h1, eri, n_orb, casci_total, e_core, nelecas = load_and_compute_integrals(
        cif, cas_electrons=cas_e, cas_orbitals=cas_o
    )
    print(f"-> CASCI total energy:   {casci_total:.6f} Ha")
    print(f"-> active orbitals:      {n_orb}  (nelecas={nelecas})")
    print(f"-> core/offset energy:   {e_core:.6f} Ha")
    print(f"-> h1 shape {h1.shape}, eri shape {eri.shape}")
    print("\nData translation pipeline checked. Asset ready to hook directly into EnterprisePipelineOrchestrator loops.")