#!/usr/bin/env python3
"""Scout probe #2 for chem-70c: how sensitive is centered-frame PDS to a badly estimated mu?

device_odmd.centered_frame computes mu (band center) from dense diagonalization -- a
validation-scale-only technique (specs/BACKLOG.md caveat: "gate PDS's sensitivity to a badly
estimated mu, not assume mu is free to compute in general"). This probe manually shifts by
mu_ideal + delta for a range of delta and checks how much PDS(7)/PDS(8) degrades, on N2 CAS(6,6).
"""
import dataclasses

import numpy as np
from qiskit.quantum_info import SparsePauliOp

from device_odmd import centered_frame
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from moment_expansion import hamiltonian_moments, pds_energy

mh = build_molecular_hamiltonian(
    atom="N 0 0 0; N 0 0 1.10", basis="sto3g",
    active_electrons=6, active_orbitals=6,
)
fci = mh.ground_state_energy()
_, _, mu_ideal = centered_frame(mh)
print(f"N2 CAS(6,6): FCI={fci:.6f}  mu_ideal={mu_ideal:.6f} Ha")

# Reachable band half-width (for context: how big is delta relative to the band itself).
w_eig, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
psi0 = np.asarray(mh.hf_state().data, dtype=complex)
pops = np.abs(V.conj().T @ psi0) ** 2
reach = w_eig[pops > 1e-8].real
half_width = 0.5 * (reach.max() - reach.min())
print(f"reachable band half-width W/2 = {half_width:.4f} Ha "
      f"(mu error as a fraction of this is the honest scale)")
print()

deltas = [0.0, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0]
header = (f"{'delta(Ha)':>10} {'delta/W':>9} {'cond_K7':>12} {'cond_K8':>12} "
          f"{'err_K7(mHa)':>12} {'err_K8(mHa)':>12}")
print(header)
print("-" * len(header))
for delta in deltas:
    mu_shift = mu_ideal + delta
    shifted = (mh.qubit_hamiltonian
               - SparsePauliOp("I" * mh.num_qubits, coeffs=[mu_shift])).simplify()
    mh_shift = dataclasses.replace(mh, qubit_hamiltonian=shifted,
                                    energy_offset=mh.energy_offset + mu_shift)
    mu, off = hamiltonian_moments(mh_shift, 16)

    def cond_at(K):
        M = np.array([[mu[2 * K - i - j] for j in range(1, K + 1)] for i in range(1, K + 1)])
        return np.linalg.cond(M)

    try:
        e7 = pds_energy(mu, 7, off)
        err7 = (e7 - fci) * 1e3
    except (ValueError, np.linalg.LinAlgError):
        err7 = float("nan")
    try:
        e8 = pds_energy(mu, 8, off)
        err8 = (e8 - fci) * 1e3
    except (ValueError, np.linalg.LinAlgError):
        err8 = float("nan")

    print(f"{delta:10.3f} {delta / half_width:9.4f} {cond_at(7):12.3e} {cond_at(8):12.3e} "
          f"{err7:12.4f} {err8:12.4f}")
