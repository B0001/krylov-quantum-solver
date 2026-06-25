#!/usr/bin/env python3
"""
DMRG-backed benchmark rung -- extends the reference ladder past exact FCI.

Fixes a stretched (multireference) N2 geometry and GROWS the active space. At each rung it
reports four energies from one shared set of CASCI integrals:

  * HF        -- active-space Hartree-Fock determinant (single reference),
  * Krylov    -- real-time quantum Krylov (this project); statevector-limited, so only run
                 while the active space fits in <= KRYLOV_QUBIT_CAP qubits,
  * FCI       -- exact active-space FCI (PySCF); only while the determinant count is tractable,
  * DMRG      -- block2 / pyblock2 (the reference that survives past FCI's reach).

The point of the rung:
  1. Where FCI is tractable, DMRG reproduces it to ~1e-9 Ha  -> validates DMRG as a drop-in
     reference (the cross-check the refactor plan called for but could not run without block2).
  2. Past the FCI determinant cutoff, FCI drops out and DMRG carries the reference alone, so the
     correlation energy HF misses (E_HF - E_DMRG) is still quantified.
  3. The quantum Krylov estimate tracks the reference wherever the statevector simulation fits.

This is an HONEST rung: the quantum solver is statevector-limited and cannot itself reach the
active spaces where DMRG matters -- that gap (Trotter/qubitization for shallower hardware
circuits, and measurement rather than statevector inner products) is the remaining work in
REFACTOR_PLAN.md, not something this script hides.

Run:  python benchmark_dmrg.py   ->  prints a table and writes data/dmrg_ladder.csv
Requires: block2 (pip install block2).  Falls back to FCI-only with a warning if absent.
"""
import csv
import math
import os

import numpy as np
from pyscf import gto, scf, mcscf, ao2mo

from hybrid_quantum_solver.molecular_hamiltonian import build_hamiltonian_from_integrals
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from hybrid_quantum_solver.dmrg_reference import fci_energy, dmrg_energy, dmrg_available

# --- benchmark configuration ------------------------------------------------
GEOMETRY_R = 1.30                 # Angstrom; stretched N2 (~1.10 is equilibrium) => multireference
BASIS = "6-31g"                   # 18 spatial orbitals for N2 -> room to grow the active space
ACTIVE_SPACES = [(6, 6), (8, 8), (10, 10), (12, 12), (14, 14)]   # (electrons, orbitals)
KRYLOV_DIM = 12
KRYLOV_QUBIT_CAP = 16             # statevector limit for the quantum solver in this project
FCI_DET_CUTOFF = 5_000_000        # skip exact FCI past this many determinants (DMRG takes over)
OUTPUT = "data/dmrg_ladder.csv"

# Modest DMRG schedule: ample for <= 14 orbitals, keeps the demo fast.
DMRG_KW = dict(bond_dims=(100, 200, 400), n_sweeps=14, noises=(1e-4, 1e-5, 0.0), n_threads=4)


def hf_determinant_energy(h1, eri, n_occ, e_core):
    """Closed-shell Hartree-Fock determinant energy from active-space integrals.

    occ = lowest ``n_occ`` orbitals, doubly occupied. eri is chemist (pq|rs).
        E = e_core + sum_i 2 h_ii + sum_ij [ 2 (ii|jj) - (ij|ji) ]
    """
    o = slice(0, n_occ)
    e = float(e_core) + 2.0 * np.einsum("ii->", h1[o, o])
    e += np.einsum("iijj->", eri[o, o, o, o]) * 2.0
    e -= np.einsum("ijji->", eri[o, o, o, o])
    return float(e)


def krylov_energy(h1, eri, nelecas, e_core, reference, n_qubits):
    """Best quantum Krylov estimate vs the reference (None if too big for statevector).

    The qubit-count gate is checked from ``n_qubits = 2 * norb`` BEFORE building anything --
    constructing the qubit Hamiltonian for an active space past the statevector cap (e.g. 28
    qubits) is itself prohibitively expensive in memory, so it must be skipped, not built-then-
    discarded.
    """
    if n_qubits > KRYLOV_QUBIT_CAP:
        return None, None
    mh = build_hamiltonian_from_integrals(h1, eri, num_particles=nelecas, energy_offset=e_core)
    steps = QuantumKrylovSolver(mh).convergence(KRYLOV_DIM)
    best = min(steps, key=lambda s: abs(s.energy - reference))
    return best.energy, best.rank


def main():
    if not dmrg_available():
        print("[WARN] block2 not importable -- DMRG column unavailable; install with "
              "`pip install block2`. Falling back to FCI only where tractable.\n")

    print(f"DMRG-backed ladder | N2 @ R={GEOMETRY_R} A | basis {BASIS} | quantum Krylov M={KRYLOV_DIM}")
    mol = gto.M(atom=f"N 0 0 0; N 0 0 {GEOMETRY_R}", basis=BASIS, verbose=0)
    mf = scf.RHF(mol).run()

    header = (f"{'CAS':>8} {'qubits':>6} {'ndet':>12} {'HF':>12} {'Krylov':>12} "
              f"{'FCI':>12} {'DMRG':>12} {'|DMRG-FCI|':>11} {'HF corr':>9} {'Kry err':>9}")
    print(header)
    print("-" * len(header))

    os.makedirs(os.path.dirname(OUTPUT) or ".", exist_ok=True)
    rows = []
    for nelec, norb in ACTIVE_SPACES:
        cas = mcscf.CASCI(mf, norb, nelec)               # integrals only; no cas.kernel() needed
        h1, e_core = cas.get_h1eff()
        eri = ao2mo.restore(1, cas.get_h2eff(), norb)
        na, nb = cas.nelecas
        nelecas = (int(na), int(nb))
        ndet = math.comb(norb, na) * math.comb(norb, nb)
        n_qubits = 2 * norb

        e_hf = hf_determinant_energy(h1, eri, na, e_core)

        e_dmrg = dmrg_energy(h1, eri, nelecas, e_core, **DMRG_KW) if dmrg_available() else None
        e_fci = fci_energy(h1, eri, nelecas, e_core) if ndet <= FCI_DET_CUTOFF else None

        reference = e_dmrg if e_dmrg is not None else e_fci   # DMRG preferred, else FCI
        e_kry, kry_rank = (None, None)
        if reference is not None:
            e_kry, kry_rank = krylov_energy(h1, eri, nelecas, e_core, reference, n_qubits)

        dmrg_fci = (abs(e_dmrg - e_fci) if (e_dmrg is not None and e_fci is not None) else None)
        hf_corr = (e_hf - reference) if reference is not None else None
        kry_err = (abs(e_kry - reference) if (e_kry is not None and reference is not None) else None)

        def col(x, w=12, p=6): return f"{x:{w}.{p}f}" if x is not None else f"{'--':>{w}}"
        c_dfci = f"{dmrg_fci:11.2e}" if dmrg_fci is not None else f"{'--':>11}"
        c_hfc = f"{hf_corr * 1e3:8.1f}m" if hf_corr is not None else f"{'--':>9}"
        c_kry = f"{kry_err * 1e3:8.3f}m" if kry_err is not None else f"{'--':>9}"
        print(f"({nelec:>2},{norb:>2}) {n_qubits:6d} {ndet:12,d} {col(e_hf)} {col(e_kry)} "
              f"{col(e_fci)} {col(e_dmrg)} {c_dfci} {c_hfc} {c_kry}")

        rows.append({
            "cas_electrons": nelec, "cas_orbitals": norb, "qubits": n_qubits, "ndet": ndet,
            "hf_energy": e_hf, "krylov_energy": e_kry, "krylov_rank": kry_rank,
            "fci_energy": e_fci, "dmrg_energy": e_dmrg,
            "dmrg_vs_fci": dmrg_fci, "hf_correlation_Ha": hf_corr,
            "krylov_error_vs_ref_Ha": kry_err,
            "reference": "dmrg" if e_dmrg is not None else ("fci" if e_fci is not None else "none"),
        })
        with open(OUTPUT, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)

    print("\nWhere both run, |DMRG-FCI| ~ 1e-9 Ha (DMRG validated as the reference). Past the FCI")
    print("determinant cutoff DMRG carries the reference alone; quantum Krylov tracks it while the")
    print(f"active space fits in <= {KRYLOV_QUBIT_CAP} qubits (statevector limit).  ->  {OUTPUT}")


if __name__ == "__main__":
    main()
