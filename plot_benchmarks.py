#!/usr/bin/env python3
import polars as pl
import matplotlib.pyplot as plt

def plot_performance(csv_file):
    df = pl.read_csv(csv_file)
    
    fig, ax1 = plt.subplots(figsize=(8, 5))
    
    # Plot Timing (Scaling)
    color = 'tab:blue'
    ax1.set_xlabel('System Size (Orbitals)')
    ax1.set_ylabel('Execution Time (s)', color=color)
    ax1.plot(df['n_orbitals'].to_numpy(), df['time_seconds'].to_numpy(), marker='o', color=color, label='Time')
    ax1.tick_params(axis='y', labelcolor=color)
    
    # Plot Energy (Verification of Stability)
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('Energy (Hartrees)', color=color)
    ax2.plot(df['n_orbitals'].to_numpy(), df['energy_hartrees'].to_numpy(), marker='x', linestyle='--', color=color, label='Energy')
    ax2.tick_params(axis='y', labelcolor=color)
    
    plt.title('Orchestrator Performance & Stability Scaling')
    plt.grid(True, alpha=0.3)
    plt.savefig('performance_scaling.svg')
    print("[SUCCESS] Performance plot saved as performance_scaling.svg")

if __name__ == "__main__":
    plot_performance("benchmark_results.csv")