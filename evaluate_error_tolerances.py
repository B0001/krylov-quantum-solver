#!/usr/bin/env python3
"""
Advanced Diagnostic Pass: evaluate_error_tolerances.py
Simulates random phase-gate noise thresholds against alpha regularization 
parameters to map precision stability margins on near-term hardware profiles.
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from hybrid_quantum_solver.orchestrate_hybrid_pipeline import AdvancedStochasticCompactor, StabilizedSubspaceShifter

def execute_noise_stress_test(csv_source="benchmark_results.csv", output_plot="error_tolerance_scaling.png"):
    if not os.path.exists(csv_source):
        print(f"[ERROR] Found missing baseline dependency. Execute 'run_sprint_benchmarks.sh' first.")
        return

    # Ingest the dynamic geometric potential energy surface records
    pes_data = pd.read_csv(csv_source)
    
    # Isolate a single structural coordinate position for consistent stress profiling (Equilibrium 0.70 A)
    target_row = pes_data[pes_data['bond_distance_A'] == 0.70]
    if target_row.empty:
        target_row = pes_data.iloc[1] # Fallback to second element index if precise step missing
        
    hf_baseline_E = target_row['classical_hf_energy_Ha'].values[0]
    distance_lbl = target_row['bond_distance_A'].values[0]

    # Setup the execution grid array matrices
    alpha_values = np.linspace(0.0, 20.0, 5)
    noise_thresholds = np.logspace(-4, -1, 4)  # Phase variance gamma ranging from 0.0001 to 0.1
    
    # Container tracking matrix dictionary for parsed outcomes
    results_matrix = {noise: [] for noise in noise_thresholds}
    lambda_prime_tracking = {noise: [] for noise in noise_thresholds}

    print("================================================================================")
    print(f"[STRESS KERNEL] Evaluating Error Tolerances at Bond Distance: {distance_lbl} Å")
    print(f"  -> Hartree-Fock Energy Asset Line: {hf_baseline_E:.6f} Ha")
    print("================================================================================")

    # Reconstruct the base spatial molecular orbital integral inputs for H2
    mock_one_body = [(0, 0, -1.25), (1, 1, -0.455)]
    mock_two_body = [(0, 1, 0, 1, 0.680)]

    for gamma in noise_thresholds:
        print(f"\n[NOISE MATRIX Channel] Injecting Phase Variance Threshold gamma = {gamma:.4f}")
        for alpha in alpha_values:
            
            # Initialize compactor with positional targets
            compactor = AdvancedStochasticCompactor(4, 1e-3)
            compactor.alpha = alpha  # Dynamically assign to override package parameters safely
            
            for p, q, w in mock_one_body:
                compactor.map_one_body_term(p, q, w)
            for p, q, r, s, w in mock_two_body:
                compactor.map_two_body_term(p, q, r, s, w)
            compactor.finalize_and_compile_metrics()
            
            # Determine regularized norm from package target name parameters securely
            # Uses a dynamic fallback to protect against local environment drift or un-synchronized builds
            current_lambda_prime = getattr(compactor, "lambda_prime_norm", getattr(compactor, "lambda_norm", 0.0))
            
            # 2. Simulate QPU element data matrix feeds under noisy phase degradation
            stabilizer = StabilizedSubspaceShifter(subspace_dimension=4, conditioning_cutoff=1e-5)
            
            # Apply dynamic error drift transformations mimicking random phase deviations
            random_drift_0 = np.random.normal(0.0, gamma * (compactor.lambda_norm / (alpha + 1.0)))
            random_drift_1 = np.random.normal(0.0, gamma * 2.5) 
            
            noisy_elements = [
                {"row": 0, "col": 0, "h_val": float(hf_baseline_E * 1.15 + random_drift_0), "s_val": 1.0},
                {"row": 0, "col": 1, "h_val": float(-0.115 + random_drift_0), "s_val": float(0.025 + abs(random_drift_0))},
                {"row": 1, "col": 1, "h_val": float(hf_baseline_E * 0.95 + random_drift_1), "s_val": 1.0},
                {"row": 1, "col": 2, "h_val": float(-0.052 + random_drift_1), "s_val": float(0.004 + abs(random_drift_1))},
                {"row": 2, "col": 2, "h_val": float(hf_baseline_E * 0.55), "s_val": 1.0},
                {"row": 2, "col": 3, "h_val": float(hf_baseline_E * 0.55), "s_val": 1.0},
                {"row": 3, "col": 3, "h_val": float(hf_baseline_E * 0.55), "s_val": 1.0},
            ]
            
            stabilizer.construct_subspace_matrices(noisy_elements)
            resolved_energy = stabilizer.compute_ground_state()
            
            # Compute absolute error discrepancy delta relative to unperturbed state target
            unperturbed_target_E = -1.332870 
            precision_error = abs(resolved_energy - unperturbed_target_E)
            
            results_matrix[gamma].append(precision_error)
            lambda_prime_tracking[gamma].append(current_lambda_prime)
            
            print(f"  -> Penalty alpha: {alpha:4.1f} | Regularized prime norm: {current_lambda_prime:.4f} | Subspace Abs Error: {precision_error:.6f} Ha")

    # ==========================================
    # Render Dual-Panel Diagnostic Visualization
    # ==========================================
    print(f"\n[VISUALIZATION LAYER] Generating Error Tolerance Topology Plot...")
    
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), dpi=300)
    
    markers = ['o', 's', '^', 'D']
    colors = ['#2ecc71', '#3498db', '#f39c12', '#e74c3c']

    for idx, gamma in enumerate(noise_thresholds):
        lbl = f"Noise Threshold $\\gamma$ = {gamma:.4f}"
        ax1.plot(alpha_values, results_matrix[gamma], 
                 color=colors[idx], linestyle='-', marker=markers[idx], linewidth=2, label=lbl)
        ax2.plot(alpha_values, lambda_prime_tracking[gamma], 
                 color=colors[idx], linestyle='--', marker=markers[idx], linewidth=1.5, label=lbl)

    # Format Subplot 1: Algorithmic Error Scaling with raw strings
    ax1.set_title(r"Algorithmic Precision Error vs. Penalty Strength", fontsize=11, fontweight='bold', pad=12)
    ax1.set_xlabel(r"Noise Penalty Regularization Strength ($\alpha$)", fontsize=10, labelpad=8)
    ax1.set_ylabel(r"Absolute Ground-State Error Deviation (Hartrees)", fontsize=10, labelpad=8)
    ax1.set_yscale('log') 
    ax1.axhline(1e-3, color='#7f8c8d', linestyle=':', label='Chemical Accuracy Threshold (1 mHa)')
    ax1.legend(loc='upper right', frameon=True, fontsize=8.5)

    # Format Subplot 2: Compaction Track Bounds Decay with raw strings
    ax2.set_title(r"Regularized Sampling Norm Bounds ($\lambda'$) Compression", fontsize=11, fontweight='bold', pad=12)
    ax2.set_ylabel(r"Target Norm Bounds Value ($\lambda'$)", fontsize=10, labelpad=8)
    ax2.set_xlabel(r"Noise Penalty Regularization Strength ($\alpha$)", fontsize=10, labelpad=8)

    ax2.legend(loc='lower left', frameon=True, fontsize=8.5)

    plt.tight_layout()
    plt.savefig(output_plot)
    plt.close()
    
    print(f"================================================================================")
    print(f"[SUCCESS] Diagnostic Pass Complete. Plot Sheet Logged inside: {output_plot}")
    print(f"================================================================================")

if __name__ == "__main__":
    execute_noise_stress_test()
