#!/usr/bin/env python3
"""
ODMD spectroscopy -- one-particle Green's-function poles and weights from survival amplitudes.

The Lehmann representation is ODMD-shaped: for a particle-removed reference a_i|ref>,

    s_k = <ref|a_i^dag e^{-i k tau H} a_i|ref> / ||a_i|ref>||^2
        = sum_n  p_n e^{-i E_n^{N-1} k tau},   p_n = |<E_n^{N-1}| a_i |ref>|^2 / ||...||^2,

so the DMD eigenphases are the (N-1)-sector eigenvalues (ionization lines), the Vandermonde
amplitudes are the spectral weights, and a_i^dag gives the electron-addition side -- the
photoemission / inverse-photoemission spectrum A(omega) from the same 1-D signals the ODMD stack
already measures. Real-time propagation of a_i|psi> for Green's functions is established (e.g.
Kosugi & Matsushita, PRA 101, 012330 (2020)); ODMD supplies the pole/weight extraction -- a
composition of pinned primitives.

FINDINGS (gated in tests/test_odmd_spectral_spec.py): poles and degeneracy-aggregated weights
match exact Lehmann to machine precision; the HF-referenced weights differ from TRUE Lehmann
(exact-ground-state reference) by ~2% on H4 -- a real, bounded approximation; the Nb3I8 dimer's
band gap min(omega+) - max(omega-) reproduces the capstone 842.44 meV exactly; and uniform
damping s -> f^k s enters |lambda| ONLY, so INTENSITIES are damping-immune along with energies
(extending SPEC_device_odmd's phase immunity to the whole spectrum).

HONEST SCOPE (specs/SPEC_odmd_spectral.md): exact statevector signals (no shot noise/circuits
here); exactly degenerate lines merge with summed weight (correct for A(omega), individual
components unresolved); lines below the amplitude floor are invisible (the visibility law);
Nb3X8 numbers are the isolated dimer's (no band broadening -- the capstone caveat).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np
import scipy.sparse as sp
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp
from scipy.sparse.linalg import expm_multiply

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from odmd import odmd_spectrum

_MAPPER = JordanWignerMapper()


@dataclass
class SpectralLine:
    """One line of A(omega). ``omega`` is relative to the reference sector's ground energy:
    removal ('-'): omega = E0^N - pole; addition ('+'): omega = pole - E0^N."""
    omega: float
    weight: float
    orbital: int
    kind: str


def ladder_operator(kind: str, i: int, n_so: int):
    """Sparse qubit matrix of the JW-mapped ladder operator: kind '-' = a_i, '+' = a_i^dag."""
    if kind not in ("-", "+"):
        raise ValueError("kind must be '-' (removal) or '+' (addition)")
    op = FermionicOp({f"{kind}_{i}": 1.0}, num_spin_orbitals=n_so)
    return _MAPPER.map(op).to_matrix(sparse=True).tocsc()


def reference_signal(mh: MolecularHamiltonian, psi_raw, n: int = 24):
    """Survival amplitude of the normalized reference; returns (s, tau, mu, nrm2).

    Each reference gets its own centered frame (mu = center, tau = pi/W of ITS reachable
    spectrum -- evolution conserves particle number, so the sector is pinned by the reference).
    Dense diagonalization for the frame: validation scale, as everywhere in the ODMD stack.
    """
    psi_raw = np.asarray(psi_raw, dtype=complex)
    nrm2 = float(np.real(np.vdot(psi_raw, psi_raw)))
    if nrm2 < 1e-12:
        raise ValueError("reference state has (near-)zero norm -- empty/full orbital")
    psi = psi_raw / np.sqrt(nrm2)
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    w_eig, V = np.linalg.eigh(H.toarray())
    pops = np.abs(V.conj().T @ psi) ** 2
    reach = w_eig[pops > 1e-10].real
    mu = float(0.5 * (reach.max() + reach.min()))
    width = float(reach.max() - reach.min())
    full_width = float(w_eig.max() - w_eig.min())
    if width < 1e-9 * max(full_width, 1.0):
        # The reference is (numerically) an EIGENSTATE -- e.g. an operator kick that lands on
        # the single symmetry-allowed state (P|psi0> on the inversion-symmetric dimer). The
        # signal is exactly constant; evolving it with tau = pi/width would demand ~1/width
        # expm_multiply substeps (the pre-fix hang -- see specs/SPEC_odmd_optical.md G2).
        return np.ones(n, dtype=complex), float(np.pi / max(full_width, 1e-12)), mu, nrm2
    tau = float(np.pi / width)
    Hs = (H - mu * sp.identity(H.shape[0], format="csc")).tocsc()
    s = np.array([psi.conj() @ expm_multiply(-1j * (k * tau) * Hs, psi) for k in range(n)])
    return s, tau, mu, nrm2


def lines_from_signal(s, tau: float, mu: float, nrm2: float, amp_floor: float = 1e-4,
                      mod_window: float = 0.2):
    """(poles, weights) from one reference signal: ODMD eigenphases (electronic frame, mu added
    back) and Vandermonde amplitudes rescaled by the reference norm ||a_i|ref>||^2.

    ``amp_floor`` filters normalized amplitudes; ``mod_window`` widens for damped signals
    (uniform damping moves |lambda| off the unit circle but not the phases or amplitudes)."""
    E, a, _ = odmd_spectrum(s, tau, cutoff=0.0, mod_window=mod_window, amp_floor=amp_floor)
    return E + mu, a * nrm2


def greens_function_lines(mh: MolecularHamiltonian, kind: str,
                          orbitals: Optional[Sequence[int]] = None, reference=None,
                          n: int = 24, amp_floor: float = 1e-4) -> List[SpectralLine]:
    """All visible A(omega) lines of one side ('-' removal / '+' addition).

    ``reference`` defaults to |HF> (HF-referenced weights -- see the honesty gate G2); pass the
    exact ground state for true Lehmann weights at validation scale."""
    psi = (np.asarray(mh.hf_state().data, dtype=complex) if reference is None
           else np.asarray(reference, dtype=complex))
    n_so = 2 * mh.num_spatial_orbitals
    # ground energy of the reference's own sector (for the omega convention)
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    w_eig, V = np.linalg.eigh(H.toarray())
    pops = np.abs(V.conj().T @ (psi / np.linalg.norm(psi))) ** 2
    e0 = float(w_eig[pops > 1e-10].min())
    lines: List[SpectralLine] = []
    for i in (orbitals if orbitals is not None else range(n_so)):
        raw = ladder_operator(kind, i, n_so) @ psi
        if float(np.real(np.vdot(raw, raw))) < 1e-10:
            continue
        s, tau, mu, nrm2 = reference_signal(mh, raw, n)
        poles, wts = lines_from_signal(s, tau, mu, nrm2, amp_floor=amp_floor)
        for p, wt in zip(poles, wts):
            omega = (e0 - p) if kind == "-" else (p - e0)
            lines.append(SpectralLine(omega=float(omega), weight=float(wt), orbital=i, kind=kind))
    return lines


def photoemission_gap(mh: MolecularHamiltonian, reference=None, n: int = 24) -> float:
    """min(omega+) - max(omega-): the particle gap read off A(omega) -- equals the charge gap
    E(N+1) + E(N-1) - 2 E(N) when the reference sector's ground is resolved."""
    rem = greens_function_lines(mh, "-", reference=reference, n=n)
    add = greens_function_lines(mh, "+", reference=reference, n=n)
    if not rem or not add:
        raise ValueError("no visible lines on one side -- raise n or lower amp_floor")
    return min(line.omega for line in add) - max(line.omega for line in rem)


if __name__ == "__main__":
    from hybrid_quantum_solver.model_hamiltonians import ModelIntegrals
    from nb3x8_gaps import NB3X8_LT_BULK, dimer_cluster_integrals

    base = dimer_cluster_integrals(**NB3X8_LT_BULK["Nb3I8"])
    mh = ModelIntegrals(base.h1, base.eri, 0.0, (1, 1), 2).to_hamiltonian()
    w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    psi_hf = np.asarray(mh.hf_state().data, dtype=complex)
    pops = np.abs(V.conj().T @ psi_hf) ** 2
    psi0 = V[:, int(np.flatnonzero(pops > 1e-8)[0])]
    print("Nb3I8 dimer photoemission spectrum A(omega) [meV] (exact N=2 reference):")
    for kind, label in (("-", "removal  (lower Hubbard band)"),
                        ("+", "addition (upper Hubbard band)")):
        agg = {}
        for line in greens_function_lines(mh, kind, reference=psi0):
            agg[round(line.omega, 2)] = agg.get(round(line.omega, 2), 0.0) + line.weight
        print(f"  {label}:")
        for om, wt in sorted(agg.items()):
            print(f"    omega = {om:9.2f}   total weight = {wt:.3f}")
    print(f"  gap = min(omega+) - max(omega-) = {photoemission_gap(mh, reference=psi0):.2f} "
          "(capstone exact charge gap: 842.44)")
