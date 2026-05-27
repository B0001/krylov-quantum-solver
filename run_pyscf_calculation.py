#!/usr/bin/env python3
"""
Production Execution Script: run_pyscf_calculation.py
Bridges classical quantum chemistry calculations (PySCF) with the 
upgraded, SVD-stabilized hybrid quantum-classical orchestrator pipeline.
"""

import sys
import json
from hybrid_quantum_solver.chemistry_gateway import PySCFDataGateway
from hybrid_quantum_solver.orchestrate_hybrid_pipeline import EnterprisePipelineOrchestrator

def main():
    print("================================================================================")
    print("[SYSTEM INITIALIZATION] Commencing Live PySCF-Orchestrator Gateway Pipeline Run")
    print("================================================================================")

    # 1. Define Molecular Target Specification (Hydrogen dimer stretch at equilibrium)
    # Target Geometry format: Element X Y Z (in Angstroms)
    h2_equilibrium_geometry = "H 0.0 0.0 0.0; H 0.0 0.0 0.7414"
    basis_set = "sto-3g"
    
    print(f"\n[PHASE 1] Configuring Classical Chemistry Baseline...")
    print(f"  -> Target Molecule: H2 Dimer")
    print(f"  -> Coordinates:\n{h2_equilibrium_geometry}")
    print(f"  -> Atomic Basis Set: {basis_set.upper()}")

    # Initialize data gateway connector
    gateway = PySCFDataGateway(molecule_geometry=h2_equilibrium_geometry, basis_set=basis_set)

    # 2. Execute Classical Restricted Hartree-Fock (RHF) Baseline
    print("\n[PHASE 2] Launching Self-Consistent Field (SCF) Mean-Field Kernels...")
    try:
        scf_reference_energy = gateway.execute_baseline_scf()
        print(f"  -> Classical Mean-Field Baseline Resolved: {scf_reference_energy:.6f} Hartrees")
    except Exception as e:
        print(f"[FATAL GATEWAY ERROR] Classical baseline computation aborted: {e}")
        sys.exit(1)

    # 3. Extract Core and Electron Repulsion Tensors
    print("\n[PHASE 3] Extracting Spatial Molecular Orbital (MO) Integrals...")
    one_body_integrals, two_body_integrals = gateway.extract_and_parse_integrals()
    
    print(f"  -> Parsed {len(one_body_integrals)} Active One-Body (h_pq) Integrals.")
    print(f"  -> Parsed {len(two_body_integrals)} Active Two-Body (h_pqrs) ERI Terms.")

    # 4. Initialize Enterprise Hybrid Pipeline Orchestrator
    # Hydrogen dimer in STO-3G yields 4 spin-orbitals (alpha/beta pairs for 2 spatial orbitals)
    # We allocate a dimension 4 Krylov subspace to process the stochastic samples
    print("\n[PHASE 4] Bootstrapping Upgraded Enterprise Hybrid Pipeline...")
    orchestrator = EnterprisePipelineOrchestrator(
        enterprise_id="LIVE_GATEWAY_H2_BENCHMARK",
        n_spin_orbitals=4,
        subspace_dim=4,
        accuracy=1e-3  # Seeking 1 milli-Hartree chemical accuracy threshold
    )

    # 5. Execute End-to-End Molecular Query
    print("\n[PHASE 5] Executing Hybrid Core Processing Loop...")
    job_receipt = orchestrator.execute_molecular_query(
        single_body=one_body_integrals,
        two_body=two_body_integrals
    )

    # 6. Outbound Telemetry Logging
    print("\n[PHASE 6] Intercepting Production Pipeline Outbound Stream Receipts:")
    if job_receipt["status"] == "SUCCESS":
        print(f"  -> Status: Job Successfully Executed.")
        print(f"  -> Classical Reference Energy:    {scf_reference_energy:.6f} Ha")
        print(f"  -> Hybrid Ground State Resolved:  {job_receipt['computed_energy']:.6f} Ha")
        print(f"  -> Global Spectral Norm Bounds (λ): {job_receipt['system_norm_lambda']:.5f}")
        
        # Verify electron correlation energy delta subtraction
        correlation_delta = job_receipt['computed_energy'] - scf_reference_energy
        print(f"  -> Extracted Electron Correlation Energy: {correlation_delta:+.6f} Hartrees")
    else:
        print(f"  -> Status: Job Terminated or Suspended by Internal Security Pipelines.")
        print(f"  -> Diagnostic Data Log: {json.dumps(job_receipt, indent=2)}")
        sys.exit(1)

    print("\n================================================================================")
    print("[GATEWAY STATUS] End-to-End Live Integration Sweep Verified: Success.")
    print("================================================================================")

if __name__ == "__main__":
    main()
