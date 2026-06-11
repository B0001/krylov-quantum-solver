#!/usr/bin/env python3
"""
Master Integration Pipeline: ASE -> PySCF -> Hybrid Quantum Orchestrator
"""

import sys
import numpy as np
from ase.io import read
from pyscf import gto, scf, ao2mo
from pyscf.data import elements as pyscf_elements
from hybrid_quantum_solver.orchestrate_hybrid_pipeline import EnterprisePipelineOrchestrator


def count_electrons_from_atom_str(atom_str: str, charge: int = 0) -> int:
    """
    Counts total electrons by summing atomic numbers from a PySCF-format
    atom string (e.g. "Nb 0.0 0.0 0.0; Nb 1.5 1.5 1.5") before any
    gto.M call, avoiding the spin/electron inconsistency RuntimeError.
    """
    total = sum(
        pyscf_elements.charge(segment.strip().split()[0])
        for segment in atom_str.split(';')
        if segment.strip()
    )
    return total - charge


def load_and_compute_integrals(cif_filepath, spin_targets=[0, 1, 2], basis_set: str = 'lanl2dz', ecp: str = 'lanl2dz'):
    """
    Dynamically determines ground state spin and runs SCF.
    """
    print(f"================================================================================")
    print(f"[CLASSICAL NODE] Initializing simulation for: {cif_filepath}")
    print(f"================================================================================\n")

    # 1. Geometry Parser
    atom_str = parse_cif_to_pyscf_format(cif_filepath)

    results = {}

    # Pre-count electrons from the atom string BEFORE calling gto.M.
    # This avoids the RuntimeError "Electron number N and spin 0 are not consistent"
    # that fires for odd-electron systems (e.g. a single Nb atom, Z=41) when spin
    # is blindly set to 0. The minimum valid spin is 1 for odd-electron systems.
    total_electrons = count_electrons_from_atom_str(atom_str, charge=0)
    is_odd = total_electrons % 2 != 0
    default_spin = 1 if is_odd else 0

    mol_dummy = gto.M(atom=atom_str, basis=basis_set, ecp=ecp, charge=0, spin=default_spin)
    current_spin = mol_dummy.spin
    print(f"Detected {total_electrons} electrons → minimum valid spin={default_spin}")

    # 2. Filter spin targets based on parity
    # For odd systems, spin must be 1, 3, 5... (1, 3, 5 unpaired e-)
    # For even systems, spin must be 0, 2, 4... (0, 2, 4 unpaired e-)
    valid_spin_targets = [s for s in spin_targets if (s % 2 != 0) == is_odd]

    print(f"-> Detected {total_electrons} electrons (Odd: {is_odd}). Scanning spins: {valid_spin_targets}")

    for spin in valid_spin_targets:
        print(f"-> Testing spin_target: {spin} (Multiplicity {spin+1})")

        mol = gto.M(atom=atom_str, basis=basis_set, ecp=ecp, charge=0, spin=spin)
        n_alpha, n_beta = mol.nelec

        # 3. Dynamic Solver Selection
        if n_alpha != n_beta:
            mf = scf.UHF(mol)
        else:
            mf = scf.RHF(mol)

        energy = mf.kernel()
        results[spin] = energy
        print(f"   Energy: {energy:.6f} Ha")

    # 4. Identify Ground State
    ground_spin = min(results, key=results.get)
    print(f"\n[SUCCESS] Ground state found at spin_target={ground_spin}")

    # Final production run
    final_mol = gto.M(atom=atom_str, basis=basis_set, ecp=ecp, charge=0, spin=ground_spin)
    final_mf = scf.UHF(final_mol) if ground_spin != 0 else scf.RHF(final_mol)
    final_mf.kernel()

    # get_eri() is RHF-only. Use mol.intor('int2e') for AO-basis ERIs,
    # which works for both RHF and UHF and respects the ECP.
    h1 = final_mf.get_hcore()
    eri = final_mol.intor('int2e')
    n_orbitals = h1.shape[0]

    return h1, eri, n_orbitals


def get_basis_set(atoms):
    """
    Dynamically assigns basis sets based on atomic number.
    Uses LANL2DZ (with ECP) for heavy elements (Z > 36),
    and 6-31g* (all-electron) for lighter elements.
    """
    basis_map = {}
    for atom in atoms:
        symbol = atom.symbol
        # atomic_number is a standard ASE property
        z = atom.number 
        
        if z > 36:  # Kr and beyond (Transition metals/Heavy elements)
            basis_map[symbol] = 'lanl2dz'
        else:
            basis_map[symbol] = '6-31g*'
            
    return basis_map


def parse_cif_to_pyscf_format(cif_filepath):
    # Pymatgen CIF-to-string logic
    # 1. Parse coordinates via ASE
    print("================================================================================")
    print(f"[CLASSICAL NODE] Parsing crystal structure from {cif_filepath}...")
    atoms = read(cif_filepath)

    # 2. Get the smart basis map
    basis_set = get_basis_set(atoms)
    print(f"[CLASSICAL NODE] Basis Set Strategy deployed: {basis_set}")

    # 3. Build the molecule
    atom_str = "; ".join([
        f"{atom.symbol} {atom.position[0]} {atom.position[1]} {atom.position[2]}"
        for atom in atoms
    ])
    return atom_str


def translate_tensors_to_orchestrator(h1: np.ndarray, eri: np.ndarray):
    """Converts dense PySCF arrays into sparse tuple lists for the orchestrator."""
    single_body = []
    for i in range(h1.shape[0]):
        for j in range(h1.shape[1]):
            if abs(h1[i, j]) > 1e-6:
                single_body.append((i, j, float(h1[i, j])))

    two_body = []
    # eri array from ao2mo is typically returned in a packed 2D format.
    # For a full 4D tensor unpacking:
    eri_4d = ao2mo.restore(1, eri, h1.shape[0])
    for p in range(eri_4d.shape[0]):
        for q in range(eri_4d.shape[1]):
            for r in range(eri_4d.shape[2]):
                for s in range(eri_4d.shape[3]):
                    if abs(eri_4d[p, q, r, s]) > 1e-6:
                        two_body.append((p, q, r, s, float(eri_4d[p, q, r, s])))

    return single_body, two_body


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python run_real_material_simulation.py <path_to_cif_file>")
        sys.exit(1)

    target_file = sys.argv[1]

    # Generate Physics Tensors
    h1, eri, n_orbitals = load_and_compute_integrals(target_file)
    single_body, two_body = translate_tensors_to_orchestrator(h1, eri)

    # Initialize Quantum Orchestrator
    orchestrator = EnterprisePipelineOrchestrator(
        enterprise_id=f"REAL_MAT_{target_file}",
        n_spin_orbitals=n_orbitals * 2,  # Spin orbitals = 2 * spatial orbitals
        subspace_dim=4                    # Keep small for initial classical test
    )

    # Execute Pipeline
    result = orchestrator.execute_molecular_query(single_body, two_body)
    print(f"\n[PIPELINE COMPLETE] Ground State Energy: {result['computed_energy']} Ha")
