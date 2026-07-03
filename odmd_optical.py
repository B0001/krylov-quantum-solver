#!/usr/bin/env python3
"""
Optical absorption and exciton binding via ODMD -- the two-particle side of odmd_spectral.

Kick the ground state with a same-sector operator O (the dipole mu for a molecule, the
polarization P = n_1 - n_2 for the dimer clusters): the survival amplitude of O|psi0> has DMD
poles at the BRIGHT excited states and weights |<E_n|O|psi0>|^2 -- the absorption spectrum, with
selection rules emerging as missing lines (dark states simply never appear in the signal). The
elastic line at omega = 0 carries <O>^2 when nonzero (Rayleigh -- physical).

On the inversion-symmetric Nb3X8 dimers, P is odd and the singlet sector has exactly ONE odd
state (the ionic-odd combination, at energy exactly U0), so the whole optical spectrum is a
single line and the optical gap has a closed form:

    omega_opt = U0 - E0,    E0 = (U0+Us)/2 - sqrt(((U0-Us)/2)^2 + 4 t^2).

Exciton binding = charge gap - optical gap (the capstone's sector-FCI charge gap minus this):
it collapses from ~Us in the atomic limit (Nb3F8, 0.986 Us) to 0.26 Us for Nb3I8 -- the exciton
unbinds with hopping. Numbers the source paper (arXiv:2501.10320) did not report.

HONEST SCOPE (specs/SPEC_odmd_optical.md): exact statevector signals; isolated-dimer optics (no
band broadening / non-density-density terms -- the capstone caveats); P is the standard stand-in
for a dipole in a geometry-free model (trends physical, absolute intensities model-defined);
"binding" is the cluster's charge-vs-neutral gap difference, not a solid-state exciton.
"""
from __future__ import annotations

import math

import numpy as np
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from nb3x8_device_gap import exact_gap
from odmd_spectral import lines_from_signal, reference_signal

_MAPPER = JordanWignerMapper()


def dimer_polarization():
    """P = n_1 - n_2 of the two-orbital dimer (JW block ordering: alpha_1 alpha_2 beta_1 beta_2)."""
    op = FermionicOp({"+_0 -_0": 1.0, "+_2 -_2": 1.0, "+_1 -_1": -1.0, "+_3 -_3": -1.0},
                     num_spin_orbitals=4)
    return _MAPPER.map(op).to_matrix(sparse=True).tocsc()


def absorption_lines(mh: MolecularHamiltonian, kick_op, reference=None, n: int = 24,
                     amp_floor: float = 1e-6):
    """(omegas, weights) of the ``kick_op``-kicked reference: omega relative to the
    reference-sector ground, weights |<E_n|O|ref>|^2. ``reference`` defaults to |HF> (pass the
    exact ground state for true absorption intensities at validation scale)."""
    psi = (np.asarray(mh.hf_state().data, dtype=complex) if reference is None
           else np.asarray(reference, dtype=complex))
    # ground energy of the reference's sector (for the omega convention)
    w_eig, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    pops = np.abs(V.conj().T @ (psi / np.linalg.norm(psi))) ** 2
    e0 = float(w_eig[pops > 1e-10].min())
    raw = kick_op @ psi
    s, tau, mu, nrm2 = reference_signal(mh, raw, n)
    poles, wts = lines_from_signal(s, tau, mu, nrm2, amp_floor=amp_floor)
    return np.asarray(poles) - e0, np.asarray(wts)


def dimer_optical_gap(U0: float, t: float, Us: float) -> float:
    """Analytic optical gap of the extended-Hubbard dimer: the only bright (odd-singlet) state
    sits at exactly U0; the even-block ground is E0 = (U0+Us)/2 - sqrt(((U0-Us)/2)^2 + 4 t^2)."""
    e0 = 0.5 * (U0 + Us) - math.sqrt(0.25 * (U0 - Us) ** 2 + 4.0 * t * t)
    return U0 - e0


def dimer_exciton_binding(U0: float, t: float, Us: float) -> float:
    """Charge gap (sector FCI, the capstone reference) minus the analytic optical gap."""
    return exact_gap(U0, t, Us) - dimer_optical_gap(U0, t, Us)


if __name__ == "__main__":
    from nb3x8_gaps import NB3X8_LT_BULK

    print("Nb3X8 dimer optics [meV] -- one bright line each (the odd singlet at U0):")
    print(f"{'material':>8} | {'optical gap':>11} | {'charge gap':>10} | {'exciton binding':>15} |"
          f" {'binding/Us':>10} | {'oscillator wt':>13}")
    from hybrid_quantum_solver.model_hamiltonians import ModelIntegrals
    from nb3x8_gaps import dimer_cluster_integrals

    for name, p in NB3X8_LT_BULK.items():
        base = dimer_cluster_integrals(**p)
        mh = ModelIntegrals(base.h1, base.eri, 0.0, (1, 1), 2).to_hamiltonian()
        w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
        psi_hf = np.asarray(mh.hf_state().data, dtype=complex)
        pops = np.abs(V.conj().T @ psi_hf) ** 2
        psi0 = V[:, int(np.flatnonzero(pops > 1e-8)[0])]
        P = dimer_polarization()
        om, wt = absorption_lines(mh, P, reference=psi0)
        print(f"{name:>8} | {dimer_optical_gap(**p):11.2f} | {exact_gap(**p):10.2f} | "
              f"{dimer_exciton_binding(**p):15.2f} | {dimer_exciton_binding(**p) / p['Us']:10.3f}"
              f" | {wt[0]:13.3e}")
