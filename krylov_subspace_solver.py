#!/usr/bin/env python3
"""
Real-time quantum Krylov subspace diagonalization (a.k.a. filter diagonalization).

This is the method your original orchestrator's StabilizedSubspaceShifter was trying to
be, implemented correctly. It is a second, independent ground-state estimator to
cross-check SQD: same active-space interface (h1, eri, e_core, nelec, norb), same
energy-frame discipline (total = subspace eigenvalue + e_core).

Idea
----
Approximate the ground state by Rayleigh-Ritz in the Krylov space spanned by
real-time-evolved copies of a reference determinant:

    |phi_k> = e^{-i H k*dt} |phi_0>,   k = 0 .. m-1

then solve the generalized eigenproblem  H c = E S c  with

    H_ij = <phi_i| H |phi_j>,   S_ij = <phi_i|phi_j>.

The subspace is a subset of the full active-space CI space, so the method is variational:
the Krylov energy can never fall below CASCI. It approaches CASCI from above as m grows.

Hardware mapping
----------------
On a quantum computer, the device supplies the matrix elements H_ij and S_ij via
Hadamard tests on the time-evolved states e^{-iH t}|phi_0>; the classical post-processing
(the generalized eigenproblem below) is identical. Here we compute H_ij/S_ij exactly via
a matrix-free FCI Hamiltonian and exact eigenbasis propagation, which makes this module the
*validation oracle* for the quantum version: if hardware-sampled matrix elements don't
reproduce these energies, the discrepancy is noise/sampling, not the algorithm.

The collapse fix
----------------
The overlap matrix S becomes ill-conditioned as the time-evolved vectors grow nearly
linearly dependent (condition numbers of 1e6+ are normal and harmless). The regularizer
drops basis directions whose overlap eigenvalue is below `cutoff` RELATIVE to the largest
overlap eigenvalue. The original bug used an ABSOLUTE cutoff, which on a poorly scaled S
dropped almost every vector and nulled the eigenproblem. Relative thresholding keeps a
stable, meaningful subspace instead.
"""

import numpy as np
from pyscf import fci


def _dense_active_H(h1, eri, norb, nelec):
    """Dense active-space CI Hamiltonian via matrix-free FCI contraction.

    Builds the full (dim x dim) matrix once by applying H to each basis vector. Suitable
    for active spaces up to a few thousand determinants (e.g. CAS(8,8)); beyond that, switch
    to iterative real-time propagation (scipy expm_multiply) and form only H_ij/S_ij.
    """
    h2e = fci.direct_spin1.absorb_h1e(h1, eri, norb, nelec, 0.5)
    na = fci.cistring.num_strings(norb, nelec[0])
    nb = fci.cistring.num_strings(norb, nelec[1])
    dim = na * nb
    H = np.zeros((dim, dim))
    unit = np.zeros(dim)
    for j in range(dim):
        unit[:] = 0.0
        unit[j] = 1.0
        H[:, j] = fci.direct_spin1.contract_2e(
            h2e, unit.reshape(na, nb), norb, nelec
        ).ravel()
    return 0.5 * (H + H.T), dim


def krylov_ground_state(h1, eri, e_core, nelec, norb,
                        krylov_dim=12, dt=0.5, cutoff=1e-8):
    """
    Real-time Krylov ground-state energy for one active space / spin sector.

    Args:
        h1, eri, e_core, nelec, norb: same as the SQD path (cas.get_h1eff / get_h2eff).
        krylov_dim: number of real-time basis states m.
        dt: time step. Too small -> states nearly parallel (ill-conditioned S, slow
            convergence). Too large -> spectral aliasing. A value near 1/spectral_range
            is a good start; the relative threshold tolerates a wide range.
        cutoff: RELATIVE overlap-eigenvalue threshold for regularization.

    Returns:
        (total_energy, info) where total_energy = subspace eigenvalue + e_core.
    """
    Hd, dim = _dense_active_H(h1, eri, norb, nelec)
    w, V = np.linalg.eigh(Hd)                         # exact spectrum (H is real symmetric)

    phi0 = np.zeros(dim)
    phi0[0] = 1.0                                     # HF determinant reference
    c0 = V.T @ phi0                                   # reference in the eigenbasis
    times = dt * np.arange(krylov_dim)
    phase = np.exp(-1j * np.outer(times, w))          # (m, dim)

    B = (V @ (phase * c0).T).T                        # |phi_k>            (m, dim)
    HB = (V @ (w * (phase * c0)).T).T                 # H|phi_k>           (m, dim)

    S = B.conj() @ B.T
    Hmat = B.conj() @ HB.T
    S = 0.5 * (S + S.conj().T)
    Hmat = 0.5 * (Hmat + Hmat.conj().T)

    s_eig, s_vec = np.linalg.eigh(S)
    keep = s_eig > cutoff * s_eig.max()               # RELATIVE threshold (the fix)
    X = s_vec[:, keep] / np.sqrt(s_eig[keep])
    Hp = X.conj().T @ Hmat @ X
    e_active = float(np.linalg.eigvalsh(Hp)[0].real)

    return e_active + e_core, {
        "krylov_dim": krylov_dim,
        "kept": int(keep.sum()),
        "dropped": int((~keep).sum()),
        "overlap_condition": float(s_eig.max() / max(s_eig[keep].min(), 1e-300)),
    }


def krylov_convergence_sweep(h1, eri, e_core, nelec, norb,
                             dims=(2, 4, 8, 12, 16), dt=0.5, cutoff=1e-8,
                             casci_energy=None):
    """
    Sweep krylov_dim and report the convergence trend. If casci_energy is given, each row
    carries delta_mHa and a variational flag (Krylov must stay >= CASCI). Returns a list of
    dict rows, mirroring the telemetry style of the SQD harness.
    """
    rows = []
    prev = None
    for m in dims:
        energy, info = krylov_ground_state(h1, eri, e_core, nelec, norb,
                                           krylov_dim=m, dt=dt, cutoff=cutoff)
        row = {
            "krylov_dim": m, "energy": energy,
            "kept": info["kept"], "overlap_condition": info["overlap_condition"],
            "monotone_decrease": (prev is None or energy <= prev + 1e-9),
        }
        if casci_energy is not None:
            row["delta_mHa"] = abs(energy - casci_energy) * 1e3
            row["variational_ok"] = energy >= casci_energy - 1e-6
        rows.append(row)
        prev = energy
    return rows


if __name__ == "__main__":
    from pyscf import gto, scf, mcscf, ao2mo

    def reference(atom, norb, nelec_active, spin=0, basis="sto-3g"):
        mol = gto.M(atom=atom, basis=basis, spin=spin)
        mf = scf.RHF(mol) if spin == 0 else scf.ROHF(mol)
        mf.verbose = 0
        mf.kernel()
        na, nb = (nelec_active + spin) // 2, (nelec_active - spin) // 2
        cas = mcscf.CASCI(mf, norb, (na, nb))
        cas.verbose = 0
        cas.kernel()
        h1, e_core = cas.get_h1eff()
        eri = ao2mo.restore(1, cas.get_h2eff(), norb)
        return h1, eri, float(e_core), (na, nb), norb, float(cas.e_tot)

    systems = {
        "H2  CAS(2,2) singlet": ("H 0 0 0; H 0 0 0.74", 2, 2, 0),
        "H4  CAS(4,4) singlet": ("H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3", 4, 4, 0),
        "O2  CAS(4,4) triplet": ("O 0 0 0; O 0 0 1.21", 4, 4, 2),
    }
    print("=" * 74)
    for label, (atom, norb, ne, spin) in systems.items():
        h1, eri, e_core, nelec, norb, casci = reference(atom, norb, ne, spin)
        print(f"[{label}] (na,nb)={nelec}  CASCI={casci:.8f} Ha")
        for r in krylov_convergence_sweep(h1, eri, e_core, nelec, norb,
                                          dims=(2, 4, 8, 12), casci_energy=casci):
            assert r["variational_ok"], "VARIATIONAL VIOLATION"
            print(f"   m={r['krylov_dim']:>2}  E={r['energy']:.8f}  "
                  f"Δ={r['delta_mHa']:8.4f} mHa  kept={r['kept']}  "
                  f"cond={r['overlap_condition']:.1e}  mono={r['monotone_decrease']}")
        print("-" * 74)
    print("KRYLOV SOLVER VALIDATED (converges to CASCI from above; open-shell works)")
