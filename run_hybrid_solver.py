#!/usr/bin/env python3
"""
SqDRIFT Hybrid Orchestrator - Execution Entry Point
Provides a generalized CLI for executing hybrid quantum-classical subspace projections.
"""

import click
import time
import csv
import os
import numpy as np

# Import your pipeline modules
from hybrid_quantum_solver.orchestrate_hybrid_pipeline import EnterprisePipelineOrchestrator
from hybrid_quantum_solver.chemistry_gateway import load_and_compute_integrals

@click.command()
@click.option('--input_file', '-i', required=True, type=click.Path(exists=True), 
              help='Path to the molecular .cif structure file (e.g., data/NbN.cif).')
@click.option('--active_space', '-a', default='8,8', show_default=True,
              help='CASSCF active space formatted as "electrons,spatial_orbitals".')
@click.option('--subspace_dim', '-d', default=16, type=int, show_default=True,
              help='Krylov subspace projection dimension (M).')
@click.option('--noise_variance', '-n', default=0.0, type=float, show_default=True,
              help='Variance for Gaussian hardware noise injection.')
@click.option('--output', '-o', default='sqdrift_telemetry.csv', show_default=True,
              help='Output CSV filename for telemetry data.')
def main(input_file, active_space, subspace_dim, noise_variance, output):
    """
    Executes a high-performance Hybrid Quantum-Classical Simulation using SqDRIFT.
    """
    click.secho("================================================================================", fg="blue", bold=True)
    click.secho(f"[INIT] SqDRIFT Hybrid Orchestrator", fg="blue", bold=True)
    click.secho("================================================================================", fg="blue", bold=True)
    
    # 1. Parse Active Space
    try:
        cas_elec, cas_orb = map(int, active_space.split(','))
    except ValueError:
        raise click.BadParameter('active_space must be in the format "electrons,orbitals", e.g., "8,8"')

    click.echo(f"-> Target Structure:  {os.path.basename(input_file)}")
    click.echo(f"-> Active Space:      CAS({cas_elec}, {cas_orb}) -> {cas_orb * 2} Qubits")
    click.echo(f"-> Krylov Dimension:  {subspace_dim}")
    click.echo(f"-> Hardware Noise:    {noise_variance}")
    
    # 2. Classical Pre-Processing (PySCF -> CASCI)
    # *Note: Ensure your load_and_compute_integrals function is updated to accept cas_elec and cas_orb*
    click.echo("\n[PHASE 1] Initializing Classical Node (PySCF Tensor Extraction)...")
    h1, eri, n_orbitals, exact_classical_energy = load_and_compute_integrals(
        input_file, 
        cas_electrons=cas_elec, 
        cas_orbitals=cas_orb
    )
    
    enterprise_name = os.path.splitext(os.path.basename(input_file))[0]
    
    # 3. Hybrid Pipeline Orchestration
    click.echo("\n[PHASE 2] Compiling Quantum Subspace...")
    orchestrator = EnterprisePipelineOrchestrator(
        enterprise_id=enterprise_name,
        n_spin_orbitals=n_orbitals * 2, # Total qubits
        subspace_dim=subspace_dim
    )
    
    # Run O(N^4) Compilation exactly once
    orchestrator.compile_target_system(h1, eri)
    
    # 4. Quantum Execution & SVD Shift
    click.echo("\n[PHASE 3] Executing QPU Oracle & Stabilized Subspace Shift...")
    start_time = time.time()
    result = orchestrator.execute_subspace_sweep(target_dim=subspace_dim, noise_variance=noise_variance)
    execution_time = time.time() - start_time
    
    computed_energy = result['computed_energy']
    energy_delta = abs(computed_energy - exact_classical_energy)
    
    # 5. Telemetry Output
    click.secho("\n================================================================================", fg="green", bold=True)
    click.secho(f"[COMPLETE] Ground State Resolved: {computed_energy:.6f} Ha", fg="green", bold=True)
    click.secho(f"           Classical Baseline:    {exact_classical_energy:.6f} Ha", fg="white")
    click.secho(f"           Delta (Error):         {energy_delta:.6e} Ha", fg="white")
    click.secho(f"           QCIVET Audit State:    {result.get('status')}", fg="yellow")
    click.secho("================================================================================", fg="green", bold=True)

    # Write Telemetry to CSV
    file_exists = os.path.isfile(output)
    with open(output, mode='a', newline='') as csv_file:
        fieldnames = ['molecule', 'qubits', 'subspace_dim', 'noise_variance', 'computed_energy', 'delta', 'time_s', 'audit_state']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        
        if not file_exists:
            writer.writeheader()
            
        writer.writerow({
            'molecule': enterprise_name,
            'qubits': cas_orb * 2,
            'subspace_dim': subspace_dim,
            'noise_variance': noise_variance,
            'computed_energy': computed_energy,
            'delta': energy_delta,
            'time_s': round(execution_time, 4),
            'audit_state': result.get('status')
        })
    click.echo(f"Telemetry appended to {output}")

if __name__ == '__main__':
    main()