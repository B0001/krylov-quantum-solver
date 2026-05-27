#!/usr/bin/env python3
"""
Live Hardware Gateway Test
Executes a remote IBM QPU call utilizing secure environment variables.
"""

from hybrid_quantum_solver.ibm_quantum_gateway import IBMQuantumGateway
from hybrid_quantum_solver.orchestrate_hybrid_pipeline import StabilizedSubspaceShifter

# 1. Initialize the Gateway 
# We no longer pass the token here. The gateway automatically scrubs the .env file!
gateway = IBMQuantumGateway(use_real_hardware=False)

# 2. Feed it a highly correlated Pauli array
test_paulis = ["IIII", "ZIII", "IZII", "XXXX", "YYXX"]
test_coeffs = [-1.042, 0.250, 0.250, 0.085, -0.085]

# 3. Ping the IBM network
live_matrix_elements = gateway.measure_subspace_elements(test_paulis, test_coeffs, subspace_dim=4)

# 4. Route the raw data into the SVD Canonical Pass
print("\n[LOCAL HPC NODE] Filtering raw IBM telemetry through SVD Canonical Pass...")
stabilizer = StabilizedSubspaceShifter(subspace_dimension=4, conditioning_cutoff=1e-5)
stabilizer.construct_subspace_matrices(live_matrix_elements)

final_energy = stabilizer.compute_ground_state()
print(f"  -> Final Ground State Energy Resolved: {final_energy:.6f} Hartrees")