#!/usr/bin/env python3
"""
Production Visualization Harness: plot_pes_curve.py
Parses logged CSV telemetry data to render publication-grade Potential 
Energy Surface (PES) curves comparing classical and hybrid solver variants.
"""

import os
import polars as pl
import matplotlib.pyplot as plt

def generate_pes_chart(csv_path="benchmark_results.csv", output_png="pes_curve_comparison.png"):
    if not os.path.exists(csv_path):
        print(f"[ERROR] Telemetry file '{csv_path}' not found. Execute your benchmark sweep first.")
        return

    # Ingest the telemetry database
    df = pl.read_csv(csv_path)

    # Configure publication styling options
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=300)

    # Plot Classical Hartree-Fock Baseline (RHF)
    ax.plot(df['bond_distance_A'].to_numpy(), df['classical_hf_energy_Ha'].to_numpy(),
            color='#e74c3c', linestyle='--', marker='o', linewidth=2, markersize=5,
            label='Classical Hartree-Fock Baseline (RHF)')

    # Plot Upgraded SVD-Stabilized Hybrid Solver Curve
    ax.plot(df['bond_distance_A'].to_numpy(), df['hybrid_ground_energy_Ha'].to_numpy(),
            color='#2c3e50', linestyle='-', marker='s', linewidth=2, markersize=5,
            label='Sample-Based Krylov Subspace Shifter')

    # Labeling and Typography configuration
    ax.set_title("Potential Energy Surface (PES) Comparison: $H_2$ Dimer Stretch", 
                 fontsize=13, fontweight='bold', pad=15)
    ax.set_xlabel("Interatomic Nuclear Distance ($R$, Angstroms $\AA$)", fontsize=11, labelpad=10)
    ax.set_ylabel("Total Ground-State Energy ($E$, Hartrees)", fontsize=11, labelpad=10)
    
    # Position the legend layout cleanly away from the data curves
    ax.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='#bdc3c7', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_png)
    plt.close()
    
    print(f"[SUCCESS] Publication-grade chart rendered and exported to: {output_png}")

if __name__ == "__main__":
    generate_pes_chart()