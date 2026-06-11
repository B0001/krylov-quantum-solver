#!/usr/bin/env python3
"""
Bridge Script: PySCF to Hybrid Quantum Solver
Extracts integrals for Nb compounds and formats them for the orchestrator.
"""

from pyscf import gto, scf, ao2mo
import numpy as np

def get_nb_integrals(atom_string: str = 'Nb 0 0 0'):
    # Define the molecule (Niobium cluster/atom)
    mol = gto.M(atom=atom_string, basis='sto-3g', spin=1, charge=0)
    mf = scf.RHF(mol).run()
    
    # Extract 1-body integrals
    h1 = mf.get_hcore()
    
    # Extract 2-body integrals (ERIs)
    eri = ao2mo.kernel(mol, mol.intor('int2e'))
    
    return h1, eri

if __name__ == "__main__":
    h1, eri = get_nb_integrals()
    print(f"Generated 1-body integral matrix shape: {h1.shape}")
    print(f"Generated 2-body integral tensor shape: {eri.shape}")
    
    # You can now feed these into your orchestrator
    np.save("data/nbn_1body.npy", h1)
    np.save("data/nbn_2body.npy", eri)
