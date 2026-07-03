#!/usr/bin/env python3
"""
Spin spectroscopy via ODMD -- the interlayer magnetic exchange J of the Nb3X8 dimers.

The third channel of the response trilogy (charge = odmd_spectral, optical = odmd_optical):
kick the ground state with the staggered magnetization S1z - S2z. That operator cannot move
charge, so on the inversion-symmetric dimer it exposes exactly the state the polarization P
leaves dark -- the m=0 triplet -- as a SINGLE line at omega = J, the singlet-triplet splitting:

    J = sqrt(((U0-Us)/2)^2 + 4 t^2) - (U0-Us)/2      (exact; -> 4t^2/(U0-Us) as t -> 0).

From the paper's own cRPA parameters this gives the interlayer exchange constants of the family
(J = 0.051 / 66.2 / 119.1 / 245.9 meV for F/Cl/Br/I -- numbers arXiv:2501.10320 did not report;
its focus was charge gaps), and a falsifiable physics rider: the textbook Heisenberg
superexchange errs by 0% (F) -> 6% (Cl) -> 14% (Br) -> 46.5% (I) -- the iodide's interlayer
dimer is beyond the Heisenberg regime. The Sz spectral weight ||Sz|psi0>||^2 is the local-moment
fraction: 1.000 (F, pure spins) -> 0.759 (I) -- charge fluctuations eat a quarter of the
iodide's moment.

Machinery is 100% reuse: absorption_lines with a spin kick; Sz|psi0> is an exact eigenstate
(Sz annihilates the ionic components), exercising the SPEC_odmd_optical degenerate-reference
short-circuit on a second operator.

HONEST SCOPE (specs/SPEC_odmd_spin.md): the ISOLATED dimer's interlayer J only (no in-plane
kagome exchange); density-density interactions only -- Nb3F8's J (0.051 meV) is below the
model's own neglected few-meV terms, so quote it as "~0", not a prediction; exact statevector.
"""
from __future__ import annotations

import math

from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from odmd_optical import absorption_lines

_MAPPER = JordanWignerMapper()


def dimer_staggered_moment():
    """S1z - S2z of the two-orbital dimer (JW block ordering: alpha_1 alpha_2 beta_1 beta_2)."""
    op = FermionicOp({"+_0 -_0": 0.5, "+_2 -_2": -0.5, "+_1 -_1": -0.5, "+_3 -_3": 0.5},
                     num_spin_orbitals=4)
    return _MAPPER.map(op).to_matrix(sparse=True).tocsc()


def spin_excitation_lines(mh: MolecularHamiltonian, reference=None, n: int = 24):
    """(omegas, weights) of the staggered-magnetization kick -- the spin excitation spectrum."""
    return absorption_lines(mh, dimer_staggered_moment(), reference=reference, n=n)


def dimer_exchange_analytic(U0: float, t: float, Us: float) -> float:
    """Exact singlet-triplet splitting of the extended-Hubbard dimer (the interlayer J)."""
    return math.sqrt(0.25 * (U0 - Us) ** 2 + 4.0 * t * t) - 0.5 * (U0 - Us)


def dimer_exchange_heisenberg(U0: float, t: float, Us: float) -> float:
    """The perturbative (Heisenberg superexchange) estimate 4 t^2 / (U0 - Us)."""
    return 4.0 * t * t / (U0 - Us)


if __name__ == "__main__":
    import numpy as np

    from hybrid_quantum_solver.model_hamiltonians import ModelIntegrals
    from nb3x8_gaps import NB3X8_LT_BULK, dimer_cluster_integrals

    print("Nb3X8 interlayer exchange from the staggered-magnetization kick [meV]:")
    print(f"{'material':>8} | {'J (ODMD=exact)':>14} | {'Heisenberg':>10} | {'Heis err':>8} | "
          f"{'local moment':>12}")
    Sz = dimer_staggered_moment()
    for name, p in NB3X8_LT_BULK.items():
        base = dimer_cluster_integrals(**p)
        mh = ModelIntegrals(base.h1, base.eri, 0.0, (1, 1), 2).to_hamiltonian()
        w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
        psi_hf = np.asarray(mh.hf_state().data, dtype=complex)
        pops = np.abs(V.conj().T @ psi_hf) ** 2
        psi0 = V[:, int(np.flatnonzero(pops > 1e-8)[0])]
        om, wt = spin_excitation_lines(mh, reference=psi0)
        J, Jh = dimer_exchange_analytic(**p), dimer_exchange_heisenberg(**p)
        print(f"{name:>8} | {om[0]:14.3f} | {Jh:10.3f} | {100 * (Jh / J - 1):7.1f}% | "
              f"{wt[0]:12.3f}")
