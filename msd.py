#!/usr/bin/env python3
"""
Mirror subspace diagonalization (MSD) -- a quantum Krylov algorithm that estimates the projected
Hamiltonian matrix from *overlap* measurements at symmetrically shifted timesteps instead of
measuring each Pauli term of H, lowering the sampling cost (Shibata/Yoshioka-style MSD,
arXiv:2511.20998).

Standard KQD evaluates the Toeplitz Hamiltonian element H_{0k} = <phi_0|H e^{-iHk tau}|phi_0> by an
LCU/Hadamard-test sum over the Pauli decomposition of H, whose sampling variance scales with the
1-norm lambda = sum_l |c_l|. MSD instead writes H = i d/dt e^{-iHt}|_{t=0}, so

    H_{0k} = i S'(k tau),   S(t) = <phi_0| e^{-iHt} |phi_0>,

and estimates the derivative by a central finite-difference stencil of the overlap function at the
*shifted* times {k tau + j delta}. The resulting estimator's variance scales with the stencil
1-norm  ||w||_1 / delta  ("fd1") rather than lambda. Two ingredients make ||w||_1/delta << lambda:
  * an **energy-level shift** H -> H - mu I that centers the reachable spectrum, so the finite
    difference only has to resolve the spectral *width* W (not the absolute energy), allowing a
    larger delta;
  * a **high-order central stencil**, which holds the finite-difference bias ~ delta^order down at a
    larger delta (delta ~ bias^{1/order}), shrinking fd1 = ||w||_1/delta.

HONEST SCOPE (see specs/SPEC_msd_sampling.md): the sampling advantage is a (lambda / W) effect. On
minimal-basis validatable systems lambda/W is O(1-3), so the realized advantage is modest (~2x RMS
=> ~4-5x shots on N2 CAS(6,6) with an order-8 stencil) and *absent* on H2/LiH (lambda/W ~ 1-2,
fd1 >= lambda). The paper's 10-10^4x regime needs much larger lambda/W than dense-diagonalizable
systems provide. What is delivered is the exact MSD construction, a faithful per-element shot-noise
model, and a quantified, falsifiable demonstration + scaling boundary.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import toeplitz
from scipy.sparse import identity
from scipy.sparse.linalg import expm_multiply

from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import solve_generalized_eig

# Central finite-difference weights for the first derivative: f'(0) ~ sum_j w_j f(j*delta) / delta.
_STENCILS = {
    2: (np.array([-1, 1]), np.array([-1 / 2, 1 / 2])),
    4: (np.array([-2, -1, 1, 2]), np.array([1 / 12, -2 / 3, 2 / 3, -1 / 12])),
    6: (np.array([-3, -2, -1, 1, 2, 3]),
        np.array([-1 / 60, 3 / 20, -3 / 4, 3 / 4, -3 / 20, 1 / 60])),
    8: (np.array([-4, -3, -2, -1, 1, 2, 3, 4]),
        np.array([1 / 280, -4 / 105, 1 / 5, -4 / 5, 4 / 5, -1 / 5, 4 / 105, -1 / 280])),
}


def central_difference_weights(order: int):
    """Offsets (in units of delta) and weights of the central first-derivative stencil of ``order``.

    ``f'(0) ~ sum_j weights[j] * f(offsets[j]*delta) / delta``; the stencil 1-norm is
    ``sum_j |weights[j]| / delta``.
    """
    if order not in _STENCILS:
        raise ValueError(f"order must be one of {sorted(_STENCILS)}")
    return _STENCILS[order]


def _overlaps(H_sparse, psi0, times):
    return np.array([psi0.conj() @ expm_multiply(-1j * t * H_sparse, psi0) for t in times])


@dataclass
class MSDProblem:
    """Noiseless KQD/MSD Toeplitz matrix elements + the per-element shot-noise scales.

    Energies are in the shifted electronic frame (H - mu); add ``offset`` to recover the total
    energy. ``ref`` is the exact (noiseless) KQD ground eigenvalue in that frame.
    """
    n: int
    tau: float
    mu: float
    delta: float
    order: int
    s: np.ndarray          # overlap elements S_{0k}, k=0..n-1 (s[0] = 1)
    h_kqd: np.ndarray      # exact H_{0k}
    h_msd: np.ndarray      # finite-difference H_{0k}
    lam: float             # Pauli 1-norm of the shifted Hamiltonian (KQD noise scale)
    fd1: float             # stencil 1-norm ||w||_1 / delta (MSD noise scale)
    width: float           # reachable spectral width W
    dim: int               # Hilbert-space dimension d
    offset: float          # mh.energy_offset + mu  -> total energy
    ref: float             # noiseless KQD ground eigenvalue (shifted frame)

    @property
    def msd_bias(self) -> float:
        """|noiseless-MSD ground - noiseless-KQD ground| -- the finite-difference bias."""
        return abs(_solve(self.h_msd, self.s, 0.0) - self.ref)


def _solve(hvec, svec, floor):
    """Ground Ritz value of the Hermitian-Toeplitz (H, S) built from their first columns."""
    H = toeplitz(np.conj(hvec), hvec)
    S = toeplitz(np.conj(svec), svec)
    return solve_generalized_eig(H, S, 1e-12, floor)[0]


def build_msd_problem(mh: MolecularHamiltonian, n: int = 8, order: int = 8,
                      delta: float | None = None, bias_target: float = 1.6e-3) -> MSDProblem:
    """Construct the exact KQD and MSD matrix elements with an energy-level shift.

    ``mu`` centers the HF-reachable spectrum; ``tau = pi / W``. If ``delta`` is None, the largest
    stencil spacing keeping the finite-difference bias below ``bias_target`` is chosen.
    """
    H_full = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    psi0 = np.asarray(mh.hf_state().data, dtype=complex)
    w_eig, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    reach = w_eig[(np.abs(V.conj().T @ psi0) ** 2) > 1e-8].real
    width = float(reach.max() - reach.min())
    mu = float(0.5 * (reach.max() + reach.min()))
    H_s = (H_full - mu * identity(H_full.shape[0], format="csc")).tocsc()

    is_id = np.array([str(p) == "I" * mh.num_qubits for p in mh.qubit_hamiltonian.paulis])
    lam = float(np.sum(np.abs(mh.qubit_hamiltonian.coeffs - mu * is_id)))
    tau = float(np.pi / width)

    s = _overlaps(H_s, psi0, [k * tau for k in range(n)])
    h_kqd = np.array([psi0.conj() @ (H_s @ expm_multiply(-1j * (k * tau) * H_s, psi0))
                      for k in range(n)])
    ref = _solve(h_kqd, s, 0.0)

    offs, wts = central_difference_weights(order)
    w1 = float(np.sum(np.abs(wts)))

    def msd_elements(d):
        return np.array([1j * np.sum(wts * _overlaps(H_s, psi0, [k * tau + o * d for o in offs])) / d
                         for k in range(n)])

    if delta is None:
        delta = 0.02
        for cand in np.linspace(0.02, 0.6, 30):
            if abs(_solve(msd_elements(cand), s, 0.0) - ref) < bias_target:
                delta = float(cand)
    h_msd = msd_elements(delta)

    return MSDProblem(n=n, tau=tau, mu=mu, delta=float(delta), order=order, s=s,
                      h_kqd=h_kqd, h_msd=h_msd, lam=lam, fd1=w1 / delta, width=width,
                      dim=int(H_full.shape[0]), offset=mh.energy_offset + mu, ref=ref)


def sample_ground_energy(prob: MSDProblem, shots: int, method: str, seed: int) -> float:
    """One shot-noisy ground-energy estimate (shifted electronic frame) for ``method`` in
    {'kqd', 'msd'}.

    Per the Hadamard-test sampling analysis (arXiv:2511.20998 Eq. 19 and the LCU variance): each
    overlap element has variance ~ 2(2 - 1/d)/shots, and the Hamiltonian element has variance
    ~ scale^2 (2 - 1/d)/shots with ``scale = lambda`` (KQD) or ``fd1`` (MSD). Noise is added to the
    ``n`` independent Toeplitz elements, then the Hermitian-Toeplitz GEVP is solved with a
    noise-aware overlap floor.
    """
    rng = np.random.default_rng(seed)
    c = 2.0 - 1.0 / prob.dim
    sig_s = np.sqrt(2.0 * c / shots)
    scale = prob.lam if method == "kqd" else prob.fd1
    sig_h = scale * np.sqrt(c / shots)
    hvec = prob.h_kqd if method == "kqd" else prob.h_msd

    def cgauss(sig, size):
        return rng.normal(0, sig / np.sqrt(2), size) + 1j * rng.normal(0, sig / np.sqrt(2), size)

    s_noisy = prob.s + np.concatenate([[0.0], cgauss(sig_s, prob.n - 1)])   # S_00 = 1 exactly
    h_noisy = hvec + cgauss(sig_h, prob.n)
    return _solve(h_noisy, s_noisy, 5.0 * sig_s)


def rms_error(prob: MSDProblem, shots: int, method: str, seeds: int = 200) -> float:
    """RMS ground-energy error (vs the noiseless KQD value) over ``seeds`` noise realizations."""
    errs = [sample_ground_energy(prob, shots, method, s) - prob.ref for s in range(seeds)]
    return float(np.sqrt(np.mean(np.square(errs))))


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    cases = {
        "H2": build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74"),
        "N2 CAS(6,6)": build_molecular_hamiltonian(
            atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
    }
    for name, mh in cases.items():
        prob = build_msd_problem(mh, n=8, order=8)
        print("=" * 72)
        print(f"{name}: lambda={prob.lam:.2f}  W={prob.width:.2f}  lambda/W={prob.lam/prob.width:.1f}"
              f"  delta={prob.delta:.3f}  fd1={prob.fd1:.2f}  MSD bias={prob.msd_bias*1e3:.2f} mHa")
        for shots in (1e5, 1e6):
            rk = rms_error(prob, int(shots), "kqd")
            rm = rms_error(prob, int(shots), "msd")
            print(f"  shots={int(shots):>8}: RMS KQD={rk:.4f}  MSD={rm:.4f}  "
                  f"advantage(KQD/MSD)={rk/rm:.2f}x")
    print("=" * 72)
