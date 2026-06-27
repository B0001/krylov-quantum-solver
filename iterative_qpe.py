#!/usr/bin/env python3
"""
Iterative Quantum Phase Estimation (single ancilla) -- the hardware-realistic FT readout.

Textbook QPE (qpe_walk_readout.py) needs a t-qubit phase register plus an inverse QFT.
Iterative QPE trades that for ONE ancilla qubit measured repeatedly: it reads the phase bits
from least- to most-significant, feeding the already-known lower bits back as a rotation that
cancels their contamination (Kitaev; Dobsicek et al. 2007). Far fewer qubits, no QFT -- the
variant you actually run on near-term-to-early-FT hardware.

    U|E> = e^{-i E t0}|E>  ->  phi = frac(-E t0 / 2pi)  ->  E = -2pi*phi / t0

The bit k (weight 2^{-k}) is read by applying U^{2^{k-1}}; the feedback angle pi*G cancels the
lower bits, with G = (b + G)/2 updated as each bit is measured. Energy resolution ~ 1/2^bits.

As with the other FT modules, the eigenphase is computed exactly here (validation oracle); on
hardware the ancilla measurements come from the device and a trial state's ground-state overlap
sets the success probability. Same active-space interface: (h1, eri, e_core, nelec, norb).

Requires qubitization_blueprint.py in the same directory.
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qubitization_blueprint import build_qubit_hamiltonian


def iqpe_phase(ph_true, n_bits, rng, shots_per_bit=15):
    """Recover phi = frac(ph_true) to n_bits, LSB->MSB, with single-ancilla feedback.
    Each bit is a majority vote over shots_per_bit single-shot measurements."""
    bits = {}
    feedback = 0.0                                  # G: accumulated lower-bit phase (turns)
    for k in range(n_bits, 0, -1):                  # bit k has weight 2^{-k}
        ones = 0
        for _ in range(shots_per_bit):
            theta = 2 * np.pi * (2 ** (k - 1)) * ph_true - np.pi * feedback
            p1 = np.sin(theta / 2.0) ** 2           # P(measure |1>)
            ones += 1 if rng.random() < p1 else 0
        b = 1 if ones > shots_per_bit / 2 else 0
        bits[k] = b
        feedback = (b + feedback) / 2.0             # G_{k-1} = (b + G_k)/2
    return sum(b * 2.0 ** (-k) for k, b in bits.items())


def iqpe_ground_energy(h1, eri, e_core, nelec, norb,
                       n_bits=12, t0=0.5, shots_per_bit=15, seed=0):
    """Iterative-QPE ground-state energy. t0 must keep the ground phase in [0,1) and avoid
    spectral aliasing; the default suits active-space Hamiltonians of a few Hartree span."""
    H, n = build_qubit_hamiltonian(h1, eri, norb)
    e0 = float(np.linalg.eigvalsh(H)[0].real)       # exact ground eigenphase (validation oracle)
    ph_true = ((-e0 * t0) / (2 * np.pi)) % 1.0
    rng = np.random.default_rng(seed)
    phi = iqpe_phase(ph_true, n_bits, rng, shots_per_bit=shots_per_bit)
    e_active = -2 * np.pi * phi / t0
    return e_active + e_core, {"n_bits": n_bits, "t0": t0, "phase": phi}


if __name__ == "__main__":
    from pyscf import gto, scf, mcscf, ao2mo

    def reference(atom, norb, ne, basis="sto-3g"):
        mol = gto.M(atom=atom, basis=basis)
        mf = scf.RHF(mol)
        mf.verbose = 0
        mf.kernel()
        cas = mcscf.CASCI(mf, norb, ne)
        cas.verbose = 0
        cas.kernel()
        h1, e_core = cas.get_h1eff()
        eri = ao2mo.restore(1, cas.get_h2eff(), norb)
        return h1, eri, float(e_core), (ne // 2, ne // 2), float(cas.e_tot)

    print("=" * 66)
    print("Iterative QPE (single ancilla) -- ground energy vs phase bits")
    for label, (atom, norb, ne) in {
        "H2  CAS(2,2)": ("H 0 0 0; H 0 0 0.74", 2, 2),
        "LiH CAS(2,2)": ("Li 0 0 0; H 0 0 1.6", 2, 2),
    }.items():
        h1, eri, e_core, nelec, casci = reference(atom, norb, ne)
        print(f"\n[{label}]  CASCI = {casci:.6f} Ha")
        for m in [4, 6, 8, 10, 12]:
            e_est, _ = iqpe_ground_energy(h1, eri, e_core, nelec, norb, n_bits=m)
            print(f"   bits={m:>2}  E_est={e_est:.6f}  err={abs(e_est - casci) * 1e3:8.3f} mHa")
    print("=" * 66)
    print("One ancilla, no QFT; precision ~ 1/2^bits. The hardware-realistic FT readout.")
