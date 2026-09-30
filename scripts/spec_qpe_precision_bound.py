#!/usr/bin/env python3
"""
Pre-registration script for specs/SPEC_qpe_precision_bound.md (chem-0ke).

Computes theta_0 = arccos(E_0/lambda) and the predicted precision-bound constant
pi*sin(theta_0) for each of the three target systems, from the EXACT spectrum alone --
no QPE sweep runs here. This is the committed prediction; tests/test_qpe_precision_bound_spec.py
runs the actual t=4..14 sweep afterwards and checks the measured max ratio against it.

Run: uv run python scripts/spec_qpe_precision_bound.py
"""
import numpy as np
from pyscf import ao2mo, gto, mcscf, scf

from qubitization_blueprint import build_qubit_hamiltonian, pauli_decompose

# CAS(n, m) here follows this repo's existing convention (taper_qubits.py, qpe_walk_readout.py):
# mcscf.CASCI(mf, n, m) -- n orbitals, m electrons.
SYSTEMS = {
    "H2 CAS(2,2)":  dict(atom="H 0 0 0; H 0 0 0.74", symmetry=None, norb=2, ne=2),
    "LiH CAS(2,2)": dict(atom="Li 0 0 0; H 0 0 1.6", symmetry=None, norb=2, ne=2),
    "N2 CAS(3,4)":  dict(atom="N 0 0 0; N 0 0 1.10", symmetry="D2h", norb=3, ne=4),
}


def build(spec):
    kwargs = dict(atom=spec["atom"], basis="sto-3g")
    if spec["symmetry"]:
        kwargs["symmetry"] = spec["symmetry"]
    mol = gto.M(**kwargs)
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    cas = mcscf.CASCI(mf, spec["norb"], spec["ne"])
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), spec["norb"])
    return h1, eri, e_core, spec["norb"], float(cas.e_tot)


def predict(spec):
    """Return (n_qubits, lambda, E0_electronic, theta0, bound) from the exact spectrum."""
    h1, eri, e_core, norb, casci = build(spec)
    H, n = build_qubit_hamiltonian(h1, eri, norb)
    Ek = np.linalg.eigvalsh(H)
    lam = sum(abs(c) for _, c in pauli_decompose(H, n))
    E0 = float(Ek[0])
    theta0 = float(np.arccos(np.clip(E0 / lam, -1, 1)))
    bound = float(np.pi * np.sin(theta0))
    return n, lam, E0, theta0, bound


if __name__ == "__main__":
    print(f"{'system':14s} {'qubits':>6} {'lambda':>10} {'E0(elec)':>12} "
          f"{'theta0(rad)':>12} {'pi*sin(theta0)':>15}")
    for name, spec in SYSTEMS.items():
        n, lam, E0, theta0, bound = predict(spec)
        print(f"{name:14s} {n:6d} {lam:10.6f} {E0:12.6f} {theta0:12.6f} {bound:15.6f}")
