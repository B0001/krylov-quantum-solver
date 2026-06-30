#!/usr/bin/env python3
"""
The rodeo algorithm -- ground/eigenstate energy estimation by a stochastic spectral filter
(Choi, Lee, Bonitati, Qian, Watkins & Lee, Phys. Rev. Lett. 127, 040505, 2021; arXiv:2009.04092).

A distinct family from subspace / variational / moment methods: each "rodeo cycle" applies a
controlled time evolution e^{-iHt} with a random time t and an ancilla phase e^{iE t}, then measures
the ancilla. Conditioning on the ancilla returning 0 multiplies the state by cos((H-E)t/2): a band-
pass filter centered on the target energy E. After K cycles with independent random times the
expected survival probability (probability all K ancillas read 0), for a reference
|psi> = sum_i c_i |E_i>, is

    P(E) = sum_i |c_i|^2  prod_{k=1}^K cos^2((E_i - E) t_k / 2).

Averaging the random times t_k ~ N(0, sigma^2) gives the closed form this module evaluates,

    P_bar(E) = sum_i |c_i|^2 [ (1 + exp(-(E_i - E)^2 sigma^2 / 2)) / 2 ]^K ,

which peaks at each eigenvalue with height |c_i|^2 (the reference overlap), the peaks narrowing as K
grows and off-resonance energies suppressed as (<1)^K. Scanning E and taking the lowest peak yields
the ground-state energy.

Reuses the validated qubit Hamiltonian and Hartree-Fock reference (``MolecularHamiltonian``); the
reference is FCI (dense diagonalization). Exact-statevector / expected-value simulation; small
systems. See specs/SPEC_rodeo.md.

Honest scope: the ground peak's height is the reference overlap |<HF|E_0>|^2, so a poor reference
gives a weak, hard-to-resolve ground peak; and the cycle count K for a target resolution grows with
the spectral range -- the cost the algorithm trades circuit repetitions for.
"""
from __future__ import annotations

import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian


def rodeo_survival(eigvals: np.ndarray, overlaps: np.ndarray, e_target: float,
                   sigma: float, n_cycles: int) -> float:
    """Expected rodeo survival probability ``P_bar(E)`` at target energy ``e_target``.

    ``eigvals``/``overlaps`` are the Hamiltonian eigenvalues and reference overlaps |<E_i|psi>|^2
    (electronic frame); ``sigma`` is the random-time standard deviation, ``n_cycles`` = K.
    """
    factor = (1.0 + np.exp(-0.5 * (np.asarray(eigvals) - e_target) ** 2 * sigma ** 2)) / 2.0
    return float(np.sum(np.asarray(overlaps) * factor ** n_cycles))


def reference_spectrum(mh: MolecularHamiltonian, overlap_tol: float = 1e-8):
    """Eigenvalues (electronic frame), reference overlaps |<E_i|HF>|^2, and the energy offset.

    Only the HF-reachable eigenvalues (overlap above ``overlap_tol``) carry rodeo signal; the full
    arrays are returned (overlaps are ~0 for unreachable states, so they do not contribute).
    """
    w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    hf = np.asarray(mh.hf_state().data, dtype=complex)
    overlaps = np.abs(V.conj().T @ hf) ** 2
    return w.real, overlaps, float(mh.energy_offset)


def rodeo_ground_energy(mh: MolecularHamiltonian, sigma: float = 2.0, n_cycles: int = 12,
                        n_grid: int = 4000, window=(-0.5, 1.5)) -> float:
    """Total ground-state energy from the lowest (here: dominant) rodeo peak.

    Scans ``e_target`` over a grid around the lowest reachable eigenvalue and returns the peak
    location plus the offset. For a Hartree-Fock reference near equilibrium the ground state carries
    the largest overlap, so the global maximum of ``P_bar`` over the low-energy window is the ground
    peak.
    """
    w, overlaps, offset = reference_spectrum(mh)
    e0 = w[overlaps > 1e-8].min()
    grid = np.linspace(e0 + window[0], e0 + window[1], n_grid)
    p = np.array([rodeo_survival(w, overlaps, e, sigma, n_cycles) for e in grid])
    return float(grid[p.argmax()]) + offset


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    for name, kw in {
        "H2": dict(atom="H 0 0 0; H 0 0 0.74"),
        "H4": dict(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"),
    }.items():
        mh = build_molecular_hamiltonian(**kw)
        fci = mh.ground_state_energy()
        w, ov, off = reference_spectrum(mh)
        e0 = w[ov > 1e-8].min()
        print(f"{name}: FCI={fci:.6f}  ground overlap |<HF|E0>|^2={ov[w.argmin()]:.3f}")
        for K in (3, 6, 12):
            e = rodeo_ground_energy(mh, n_cycles=K)
            bg = rodeo_survival(w, ov, e0 + 0.7, 2.0, K)
            print(f"   K={K:2d}: E_rodeo={e:.6f}  err={(e - fci) * 1e3:+.2f} mHa  "
                  f"off-res background={bg:.4f}")
