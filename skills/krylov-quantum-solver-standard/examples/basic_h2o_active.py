"""qDRIFT real-time Krylov for H2O/STO-3G in an active space, checked against PySCF CASCI.

Standalone (pyscf + qiskit-nature + numpy + scipy). Run:  python basic_h2o_active.py [ne no L]
Exits non-zero if any self-check fails -- a result that cannot be checked is not reported.

Method: one random walk of single-Pauli rotations exp(-i tau sgn(c_j) P_j), j ~ |c_j|/lambda,
tau = lambda*dt/L, supplies a WHOLE Krylov basis |phi_k> (read off every L rotations from the
Hartree-Fock state). H_ij, S_ij are averaged over R walks and the generalized eigenproblem is
solved by thresholded canonical orthogonalization. The averaged pair is a Gram pair, so the
estimate is variationally bounded below; postselecting each state onto the reference's
(N_alpha, N_beta) sector makes that bound the N-electron CASCI energy.
"""
import sys

import numpy as np
from pyscf import gto, mcscf, scf
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit_nature.second_q.circuit.library import HartreeFock
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer

GEO = "O 0 0 0; H 0 0.757 0.586; H 0 -0.757 0.586"
CHEM_ACC = 1.6e-3  # Ha


def active_space_hamiltonian(geo, basis, n_el, n_orb):
    """JW qubit Hamiltonian of the CAS(n_el, n_orb) problem, its constant offset, and the HF state."""
    problem = PySCFDriver(atom=geo, basis=basis).run()
    problem = ActiveSpaceTransformer(n_el, n_orb).transform(problem)
    mapper = JordanWignerMapper()
    op = mapper.map(problem.hamiltonian.second_q_op())
    offset = float(sum(problem.hamiltonian.constants.values()))  # nuclear + frozen core
    hf = Statevector(HartreeFock(problem.num_spatial_orbitals, problem.num_particles, mapper)).data
    return op, offset, hf, problem.num_spatial_orbitals, problem.num_particles


def casci_reference(geo, basis, n_el, n_orb):
    """Independent reference: PySCF CASCI on the same RHF orbitals (same active space)."""
    mf = scf.RHF(gto.M(atom=geo, basis=basis, verbose=0)).run()
    return float(mcscf.CASCI(mf, n_orb, n_el).kernel()[0])


def sector_mask(n_orb, n_particles):
    idx = np.arange(2 ** (2 * n_orb))
    count = np.vectorize(lambda x: int(x).bit_count())
    return (count(idx & ((1 << n_orb) - 1)) == n_particles[0]) & (count(idx >> n_orb) == n_particles[1])


def thresholded_ground(H, S, rel_threshold=1e-6):
    """Lowest root of H c = E S c after dropping overlap directions below rel_threshold*max."""
    s, U = np.linalg.eigh(S)
    keep = s > rel_threshold * s.max()
    X = U[:, keep] / np.sqrt(s[keep])
    Hp = X.conj().T @ H @ X
    return float(np.linalg.eigvalsh(0.5 * (Hp + Hp.conj().T))[0]), int(keep.sum())


def qdrift_krylov(op, psi0, mask, dt, M=6, L=2000, R=20, seed=0, postselect=True):
    coeffs = np.asarray(op.coeffs)
    assert np.max(np.abs(coeffs.imag)) < 1e-10, "qDRIFT here needs real Pauli coefficients"
    w = np.abs(coeffs.real)
    lam = w.sum()
    p, sgn = w / lam, np.sign(coeffs.real)
    P = [SparsePauliOp(op.paulis[j]).to_matrix(sparse=True) for j in range(len(w))]
    Hs = op.to_matrix(sparse=True)
    tau = lam * dt / L
    c, s = np.cos(tau), np.sin(tau)
    rng = np.random.default_rng(seed)
    H = np.zeros((M, M), complex)
    S = np.zeros((M, M), complex)
    for _ in range(R):
        x, cols = psi0, [psi0]
        for _k in range(1, M):
            for j in rng.choice(len(w), size=L, p=p):
                x = c * x - 1j * sgn[j] * s * (P[j] @ x)          # exp(-i tau sgn_j P_j) x
            cols.append(x)
        V = np.array(cols).T
        if postselect:
            V = V * mask[:, None]
        S += V.conj().T @ V / R
        H += V.conj().T @ (Hs @ V) / R
    return thresholded_ground(0.5 * (H + H.conj().T), 0.5 * (S + S.conj().T)), lam


def main(n_el=4, n_orb=4, L=2000):
    op, offset, hf, norb, npart = active_space_hamiltonian(GEO, "sto3g", n_el, n_orb)
    mask = sector_mask(norb, npart)
    assert np.sum(np.abs(hf[~mask]) ** 2) < 1e-12, "reference must be in the N-electron sector"
    Hd = op.to_matrix()
    e_sector = float(np.linalg.eigvalsh(Hd[np.ix_(mask, mask)])[0]) + offset
    e_casci = casci_reference(GEO, "sto3g", n_el, n_orb)
    w = np.linalg.eigvalsh(Hd)                         # dense: fine for <= ~12 qubits
    dt = np.pi / (w[-1] - w[0])
    (e, rank), lam = qdrift_krylov(op, hf, mask, dt, L=L)
    E = e + offset
    print(f"CAS({n_el},{n_orb}) {op.num_qubits} qubits, lambda = {lam:.3f}")
    print(f"PySCF CASCI           : {e_casci:.6f} Ha")
    print(f"exact qubit sector min: {e_sector:.6f} Ha   (|diff| {abs(e_sector - e_casci):.1e})")
    print(f"qDRIFT Krylov (L={L}) : {E:.6f} Ha   error {E - e_sector:+.2e}, rank {rank}")
    checks = {
        "qubit Hamiltonian reproduces CASCI (1e-6 Ha)": abs(e_sector - e_casci) < 1e-6,
        "variational: E >= E_CASCI": E >= e_sector - 1e-9,
        "within chemical accuracy (1.6 mHa)": E - e_sector < CHEM_ACC,
    }
    for name, ok in checks.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    return all(checks.values())


if __name__ == "__main__":
    args = [int(a) for a in sys.argv[1:4]]
    sys.exit(0 if main(*args) else 1)
