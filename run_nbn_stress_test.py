#!/usr/bin/env python3
"""
Sprint 8: Algorithmic Stress Testing & Subspace Scaling
Target: NbN (Niobium Nitride)
"""

import sys
import time
import numpy as np
import pandas as pd
from ase.io import read
from pyscf import (
    gto,
    scf,
    ao2mo,
    mcscf,
)

from hybrid_quantum_solver.orchestrate_hybrid_pipeline import EnterprisePipelineOrchestrator


def get_smart_basis(atoms):
    """Assigns ECP to heavy atoms, all-electron to light atoms."""
    basis_map = {}
    for atom in atoms:
        if atom.number > 36:
            basis_map[atom.symbol] = 'lanl2dz'
        else:
            basis_map[atom.symbol] = '6-31g*'
    return basis_map


def load_and_compute_integrals(cif_filepath, spin_targets=[0, 1, 2], cas_electrons=8, cas_orbitals=8):

    print(f"================================================================================")
    print(f"[CLASSICAL PRE-PROCESSING] Generating Hamiltonian for: {cif_filepath}")
    
    atoms = read(cif_filepath)
    atom_str = "; ".join([f"{a.symbol} {a.position[0]} {a.position[1]} {a.position[2]}" for a in atoms])
    basis_set = get_smart_basis(atoms)
    ecp_dict = {a.symbol: 'lanl2dz' for a in atoms if a.number > 36}
    
    # 1. Determine ground state spin
    mol_dummy = gto.M(atom=atom_str, basis=basis_set, ecp=ecp_dict, charge=0, spin=None)
    total_electrons = sum(mol_dummy.nelec)
    is_odd = total_electrons % 2 != 0
    valid_spin_targets = [s for s in spin_targets if (s % 2 != 0) == is_odd]
    
    results = {}
    for spin in valid_spin_targets:
        mol = gto.M(atom=atom_str, basis=basis_set, ecp=ecp_dict, charge=0, spin=spin)
        mf = scf.UHF(mol) if mol.nelec[0] != mol.nelec[1] else scf.RHF(mol)
        mf.max_cycle = 200 # Stabilize early checks
        results[spin] = mf.kernel()

    ground_spin = min(results, key=results.get)
    print(f"[SUCCESS] Classical baseline established. Ground state spin: {ground_spin}")
    
    # 2. Final production run
    final_mol = gto.M(atom=atom_str, basis=basis_set, ecp=ecp_dict, charge=0, spin=ground_spin)
    final_mf = scf.UHF(final_mol) if final_mol.nelec[0] != final_mol.nelec[1] else scf.RHF(final_mol)
    final_mf.max_cycle = 200
    final_mf.kernel()
    
    # [NEW] Active Space Truncation (CASSCF)
    # We select 8 spatial orbitals and 8 active electrons (CAS(8,8))
    # This translates to exactly 16 spin-orbitals (16 qubits) - perfect for simulation!
    # [NEW] Active Space Truncation (CASCI)
    print("[CLASSICAL NODE] Truncating to CAS(8,8) Active Space...")
    cas = mcscf.CASCI(final_mf, cas_orbitals, cas_electrons)
    
    # Execute the classical active space solver to establish your exact baseline
    cas.kernel()
    cas_energy = cas.e_tot

    # Extract the effective 1-body Hamiltonian and the frozen core energy
    h1, e_core = cas.get_h1eff()
    
    # Extract the effective 2-body integrals for the active space
    eri_cas = cas.get_h2eff()
    
    # PySCF returns the ERI in a condensed 1D array. Restore it to 4D (pq|rs).
    eri = ao2mo.restore(1, eri_cas, 8) 
    
    n_orbitals = h1.shape[0] # Will be exactly 8 (16 spin-orbitals)
    
    return h1, eri, n_orbitals, cas_energy

def translate_tensors_to_orchestrator(h1: np.ndarray, eri: np.ndarray):
    # Standard translation logic
    single_body = [(i, j, float(h1[i, j])) for i in range(h1.shape[0]) for j in range(h1.shape[1]) if abs(h1[i, j]) > 1e-6]
    
    two_body = []
    eri_4d = ao2mo.restore(1, eri, h1.shape[0])
    for p in range(eri_4d.shape[0]):
        for q in range(eri_4d.shape[1]):
            for r in range(eri_4d.shape[2]):
                for s in range(eri_4d.shape[3]):
                    if abs(eri_4d[p, q, r, s]) > 1e-6:
                        two_body.append((p, q, r, s, float(eri_4d[p, q, r, s])))
    return single_body, two_body

if __name__ == "__main__":
    target_file = "data/nb_structures/NbN_mp-2634.cif"
    
    # 1. Classical Setup (Done once)
    h1, eri, n_orbitals, exact_classical_energy = load_and_compute_integrals(target_file)
    single_body, two_body = translate_tensors_to_orchestrator(h1, eri)
    
    # 2. Experimental Parameters
    subspace_dimensions = [4, 8, 16, 32]
    noise_variances = [0.0, 0.01, 0.05, 0.1, 0.2]
    
    telemetry_data = []

    print(f"\n================================================================================")
    print(f"[SPRINT 8] Initiating 2D Sweep: Subspace vs. Hardware Noise")
    print(f"================================================================================")

    # 1. NEW: Initialize the orchestrator outside the loop
    orchestrator = EnterprisePipelineOrchestrator(
        enterprise_id="NbN_STRESS_TEST",
        n_spin_orbitals=n_orbitals * 2,
        subspace_dim=4 # Temporary, overridden in loop
    )
    
    # 2. NEW: Compile the target system ONCE (The heavy O(N^4) operation)
    orchestrator.compile_target_system(single_body, two_body)

    for dim in subspace_dimensions:
        for noise in noise_variances:
            print(f"\n-> Running: Subspace Dim [{dim}] | Noise Variance [{noise}]")
            
            start_time = time.time()
            
            # 3. NEW: Call the fast subspace sweep inside the loop
            result = orchestrator.execute_subspace_sweep(target_dim=dim, noise_variance=noise)
            
            execution_time = time.time() - start_time
            
            computed_energy = result['computed_energy']
            energy_delta = abs(computed_energy - exact_classical_energy)
            
            telemetry_data.append({
                'subspace_dim': dim,
                'noise_variance': noise,
                'computed_energy': computed_energy,
                'energy_delta': energy_delta,
                'execution_time_s': execution_time,
                'state': result.get('status', 'STABILIZED')
            })

    # 3. Export to CSV for Zettelkasten/Grafana
    df = pd.DataFrame(telemetry_data)
    df.to_csv("sprint_8_nbn_stress_test.csv", index=False)
    
    print(f"\n================================================================================")
    print(f"[COMPLETE] Sweep finished. Telemetry saved to 'sprint_8_nbn_stress_test.csv'")
