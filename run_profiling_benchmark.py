#!/usr/bin/env python3
"""
Sprint 9: Dynamic Profiling Benchmark
Scales the electronic Hamiltonian complexity with orbital count 
to measure real-world performance scaling.
"""

import time
import csv
import numpy as np
from hybrid_quantum_solver.orchestrate_hybrid_pipeline import EnterprisePipelineOrchestrator

def generate_dynamic_integrals(n_orbitals: int):
    """Generates a Hamiltonian that grows in complexity with N."""
    # Kinetic energy terms grow with orbital index
    single = [(i, i, -0.5 * (1.0 + 0.05 * i)) for i in range(n_orbitals)]
    
    # 2-body Coulomb interactions scale with orbital count
    two = []
    for i in range(n_orbitals - 1):
        two.append((i, i+1, i, i+1, 0.1))
    return single, two

def profile_system(n_orbitals: int, subspace_dim: int):
    # Initialize the orchestrator with dynamic dimensions
    orchestrator = EnterprisePipelineOrchestrator(
        enterprise_id=f"SCALING_RUN_{n_orbitals}",
        n_spin_orbitals=n_orbitals,
        subspace_dim=subspace_dim
    )
    
    single, two = generate_dynamic_integrals(n_orbitals)
    
    # Run loop to average out OS noise
    iterations = 20
    times = []
    for _ in range(iterations):
        start = time.time()
        result = orchestrator.execute_molecular_query(single, two)
        times.append(time.time() - start)
        
    return np.mean(times), result['computed_energy']

if __name__ == "__main__":
    results_file = "benchmark_results.csv"
    with open(results_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["n_orbitals", "time_seconds", "energy_hartrees"])
        
        # Test across increasing orbital counts
        for n in [4, 8, 16, 24, 32]:
            print(f"Scaling benchmark: N={n} orbitals...")
            t, e = profile_system(n, subspace_dim=4)
            writer.writerow([n, t, e])
            print(f"  -> Mean execution time: {t:.4f}s")
            
    print(f"\n[SPRINT 9 COMPLETE] Dynamic benchmarks serialized to {results_file}")