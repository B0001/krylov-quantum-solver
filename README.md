# SqDRIFT Hybrid Orchestrator
**A High-Performance Hybrid Quantum-Classical Pipeline for Strongly Correlated Electron Systems**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Qiskit](https://img.shields.io/badge/Qiskit-1.0+-6133BD.svg)](https://qiskit.org/)
[![PySCF](https://img.shields.io/badge/PySCF-2.4+-red.svg)](https://pyscf.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Abstract
The **SqDRIFT Hybrid Orchestrator** is an enterprise-grade computational chemistry pipeline designed to bypass the classical Full Configuration Interaction (FCI) wall. It bridges high-performance classical tensor processing (e.g., PySCF on HPC clusters) with quantum algorithmic execution (Qiskit/QPU). 

Originally developed to simulate complex topological anomalies (e.g., 90° twisted $\pi$-orbital half-Möbius systems) and highly correlated transition-metal superconductors (e.g., Niobium Nitride), this architecture natively handles multi-body Jordan-Wigner mappings, compiles Stochastic qDRIFT sequences, and resolves exact ground states via quantum Krylov subspace projection.

## Core Architecture

This framework decouples classical pre-processing from quantum execution to ensure maximum resource efficiency across heterogeneous compute nodes (such as SLURM-managed HPCs and cloud-based QPUs). 

The pipeline consists of three interconnected modules:

### 1. The Advanced Stochastic Compactor
Executes $O(N^4)$ classical pre-processing.
* Ingests molecular `.cif` structures and extracts exact classical integrands via PySCF.
* Supports active space truncation (CASCI/CASSCF) to map valence electrons into computable qubit spaces.
* Performs rigorous multi-body Jordan-Wigner transformations to generate Pauli strings and bounds the spectral norm ($\lambda$).
* Compiles the Hamiltonian into a weight-proportional stochastic qDRIFT instruction set.

### 2. The Qiskit Krylov Oracle
Executes the quantum simulation and measurement grid.
* Ingests the qDRIFT instruction set and synthesizes abstract Pauli evolutions into native quantum operations (`LieTrotter` synthesis).
* Generates the $M$-dimensional Krylov basis states via parameterized Trotter steps ($e^{-i \tau \hat{P}}$).
* Measures the $O(M^2)$ tensor grid to extract exact Hamiltonian ($H_{ij}$) and Overlap ($S_{ij}$) expectation values.

### 3. QCIVET Guard & Subspace Shifter
Ensures network transit integrity and resolves the final energy.
* **QCIVET** (Quantum Contract-based Integrity Verification and Error Transit): A cryptographic auditing layer that hashes outbound compiled slices and semantically verifies the symmetry of inbound quantum matrices to detect hardware noise drift.
* **Stabilized Subspace Shifter**: Ingests the verified $H_{ij}$ and $S_{ij}$ arrays and performs a regularized Singular Value Decomposition (SVD) to drop linearly dependent/noise-corrupted basis dimensions, finally solving the generalized eigenvalue problem $H\mathbf{c} = ES\mathbf{c}$.

## Installation

We recommend using a dedicated Conda environment to manage dependencies:

```bash
# Create and activate environment
conda create -n sqdrift-env python=3.11
conda activate sqdrift-env

# Clone the repository
git clone [https://github.com/yourusername/sqdrift-hybrid-orchestrator.git](https://github.com/yourusername/sqdrift-hybrid-orchestrator.git)
cd sqdrift-hybrid-orchestrator

# Install core quantum and chemistry libraries
pip install -r requirements.txt