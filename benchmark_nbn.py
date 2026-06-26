#!/usr/bin/env python3
"""
NbN transition-metal active-space study -- where the statevector quantum solver runs out of room
and DMRG takes over.

The quantum Krylov solver here is a STATEVECTOR simulation, capped at ~16 qubits = CAS(_,8). A Nb
d-manifold wants more orbitals than that. This script holds the (one, expensive) NbN Hartree-Fock
fixed and GROWS the active space, reporting at each rung:

  * Krylov   -- real-time quantum Krylov (this project); only while <= KRYLOV_QUBIT_CAP qubits,
  * FCI      -- exact active-space FCI (PySCF); only while the determinant count is tractable,
  * DMRG     -- block2 reference; runs at every rung (this is the point).

The scientific question it answers: how much active-space correlation lies BEYOND the CAS(8,8)
ceiling the statevector quantum method is limited to?  ``corr_beyond_cap`` = E_DMRG(rung) -
E_DMRG(smallest rung) quantifies it.

HONEST CAVEATS (see also REFACTOR_PLAN.md):
  * This is NbN as a FINITE 2-atom cluster (CIF unit cell, no periodic boundary conditions) --
    NOT the periodic superconductor. It says nothing about superconductivity.
  * At 16 qubits the exact answer (CASCI) is classically trivial, so there is NO quantum
    advantage here; this is a methodology study, not a discovery.
  * block2 (DMRG) required: pip install block2.

Run:  python benchmark_nbn.py   ->  prints a table and writes data/nbn_active_space.csv
"""
import csv
import math
import os

import numpy as np
from ase.io import read as ase_read
from pyscf import gto, scf, mcscf, ao2mo, lib

from hybrid_quantum_solver.chemistry_gateway import _get_smart_basis
from hybrid_quantum_solver.molecular_hamiltonian import build_hamiltonian_from_integrals
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from hybrid_quantum_solver.dmrg_reference import fci_energy, dmrg_energy, dmrg_available

CIF = "data/nb_structures/NbN_mp-2634.cif"
ACTIVE_SPACES = [(8, 8), (10, 10), (12, 12), (14, 14)]   # (electrons, orbitals)
SPIN_TARGETS = [0, 2, 4, 6]
KRYLOV_DIM = 8
KRYLOV_QUBIT_CAP = 16
FCI_DET_CUTOFF = 5_000_000
OUTPUT = "data/nbn_active_space.csv"
CHKFILE = "data/nbn_scf.chk"          # cached ground-state SCF for resume (skips the spin scan)
DMRG_KW = dict(bond_dims=(150, 300, 500), n_sweeps=16, noises=(1e-4, 1e-5, 0.0), n_threads=4)

FIELDS = ["cas_electrons", "cas_orbitals", "qubits", "ndet", "krylov_energy",
          "fci_energy", "dmrg_energy", "dmrg_vs_fci", "corr_beyond_cap_Ha", "reference"]


def load_existing(path):
    """Resume support: return (done_keys, cap_reference_energy) from a prior partial CSV.

    ``cap_reference_energy`` is the reference of the smallest active space (the CAS(8,8) baseline
    for the ``corr_beyond_cap`` column), recovered so the metric stays consistent across resumes.
    """
    done, cap = set(), None
    if not os.path.exists(path):
        return done, cap
    with open(path) as f:
        for row in csv.DictReader(f):
            key = (int(row["cas_electrons"]), int(row["cas_orbitals"]))
            done.add(key)
            if key == ACTIVE_SPACES[0]:
                ref = row.get("dmrg_energy") or row.get("fci_energy")
                cap = float(ref) if ref else None
    return done, cap


def _tight_scf(mf):
    """Damped DIIS then second-order (Newton/SOSCF) tightening, following any internal instability.

    Plain DIIS stalls on open-shell d-metal UHF ('not converged after 500 cycles'). Newton/SOSCF
    converges it quadratically from the damped solution; the stability loop guards against settling
    on a saddle point (a real risk for transition metals) by following a lower broken-symmetry
    solution if one exists.
    """
    mf.max_cycle, mf.level_shift, mf.diis_space = 500, 0.5, 12
    mf.conv_tol, mf.init_guess = 1e-8, "atom"
    mf.kernel()

    mf = mf.newton()
    mf.conv_tol = 1e-10
    mf.kernel(mf.mo_coeff, mf.mo_occ)

    for _ in range(3):
        try:
            mo = mf.stability(return_status=True)[0]
        except Exception:
            break
        if np.allclose(np.asarray(mo), np.asarray(mf.mo_coeff)):
            break  # internally stable
        dm = mf.make_rdm1(mo, mf.mo_occ)
        mf = mf.newton()
        mf.kernel(dm0=dm)
    return mf


def ground_state_mf(cif):
    """One transition-metal-aware SCF at the lowest-energy spin (reused for every active space).

    Resumable: the converged ground-state SCF is cached to CHKFILE and restored on rerun, so the
    (slow) spin scan is not repeated.
    """
    if os.path.exists(CHKFILE):
        try:
            mol = lib.chkfile.load_mol(CHKFILE)
            res = lib.chkfile.load(CHKFILE, "scf")
            mf = scf.UHF(mol) if mol.nelec[0] != mol.nelec[1] else scf.RHF(mol)
            mf.mo_coeff, mf.mo_occ = res["mo_coeff"], res["mo_occ"]
            mf.mo_energy, mf.e_tot = res["mo_energy"], res["e_tot"]
            print(f"[CLASSICAL] restored SCF from {CHKFILE} "
                  f"(spin={mol.spin}, e_tot={mf.e_tot:.6f} Ha)")
            return mf
        except Exception as exc:
            print(f"[CLASSICAL] chkfile restore failed ({exc}); recomputing SCF.")

    atoms = ase_read(cif)
    atom_str = "; ".join(f"{a.symbol} {a.position[0]} {a.position[1]} {a.position[2]}" for a in atoms)
    basis = _get_smart_basis(atoms)
    ecp = {a.symbol: "lanl2dz" for a in atoms if a.number > 36}

    dummy = gto.M(atom=atom_str, basis=basis, ecp=ecp, charge=0, spin=None)
    n_elec = sum(dummy.nelec)
    is_odd = n_elec % 2 != 0
    spins = [s for s in SPIN_TARGETS if (s % 2 != 0) == is_odd and s <= n_elec] or [1 if is_odd else 0]

    def build(spin):
        mol = gto.M(atom=atom_str, basis=basis, ecp=ecp, charge=0, spin=spin)
        mf = scf.UHF(mol) if mol.nelec[0] != mol.nelec[1] else scf.RHF(mol)
        mf.max_cycle, mf.level_shift, mf.diis_space = 500, 0.3, 12
        mf.conv_tol, mf.init_guess = 1e-8, "atom"
        return mf

    # Cheap damped SCF per spin just to rank them; tighten only the winner.
    energies = {s: build(s).kernel() for s in spins}
    ground = min(energies, key=energies.get)
    print(f"[CLASSICAL] {n_elec} electrons; ground spin = {ground} (2S+1 = {ground + 1})")

    mol = gto.M(atom=atom_str, basis=basis, ecp=ecp, charge=0, spin=ground)
    mf = _tight_scf(scf.UHF(mol) if mol.nelec[0] != mol.nelec[1] else scf.RHF(mol))
    try:
        s2 = mf.spin_square()[0]
    except Exception:
        s2 = float("nan")
    print(f"[CLASSICAL] tightened SCF: e_tot={mf.e_tot:.8f} Ha  converged={mf.converged}  "
          f"<S^2>={s2:.4f}")
    scf.chkfile.dump_scf(mol, CHKFILE, mf.e_tot, mf.mo_energy, mf.mo_coeff, mf.mo_occ)
    return mf


def main():
    if not dmrg_available():
        print("[WARN] block2 not importable -- install with `pip install block2`. "
              "DMRG column unavailable.\n")
    print(f"NbN active-space growth | {CIF} | quantum Krylov M={KRYLOV_DIM} (statevector <= {KRYLOV_QUBIT_CAP} qubits)")

    os.makedirs(os.path.dirname(OUTPUT) or ".", exist_ok=True)
    done, e_dmrg_cap = load_existing(OUTPUT)
    if all(a in done for a in ACTIVE_SPACES):
        print(f"All {len(ACTIVE_SPACES)} rungs already in {OUTPUT}; nothing to recompute.")
        return
    mf = ground_state_mf(CIF)

    header = (f"{'CAS':>8} {'qubits':>6} {'ndet':>14} {'Krylov':>12} {'FCI/CASCI':>12} "
              f"{'DMRG':>12} {'|DMRG-FCI|':>11} {'corr>cap':>10}")
    print(header)
    print("-" * len(header))

    file_exists = os.path.exists(OUTPUT)
    for nelec, norb in ACTIVE_SPACES:
        if (nelec, norb) in done:
            print(f"({nelec:>2},{norb:>2}) already in {OUTPUT} -- skipping (resumed)")
            continue
        cas = mcscf.CASCI(mf, norb, nelec)
        h1, e_core = cas.get_h1eff()
        eri = ao2mo.restore(1, cas.get_h2eff(), norb)
        na, nb = cas.nelecas
        nelecas = (int(na), int(nb))
        ndet = math.comb(norb, na) * math.comb(norb, nb)
        nq = 2 * norb

        e_dmrg = dmrg_energy(h1, eri, nelecas, e_core, **DMRG_KW) if dmrg_available() else None
        e_fci = fci_energy(h1, eri, nelecas, e_core) if ndet <= FCI_DET_CUTOFF else None
        reference = e_dmrg if e_dmrg is not None else e_fci

        e_kry = None
        if nq <= KRYLOV_QUBIT_CAP and reference is not None:
            mh = build_hamiltonian_from_integrals(h1, eri, num_particles=nelecas, energy_offset=e_core)
            steps = QuantumKrylovSolver(mh).convergence(KRYLOV_DIM)
            e_kry = min(steps, key=lambda s: abs(s.energy - reference)).energy

        if e_dmrg_cap is None and reference is not None:
            e_dmrg_cap = e_dmrg if e_dmrg is not None else e_fci
        corr_beyond = (reference - e_dmrg_cap) if (reference is not None and e_dmrg_cap is not None) else None
        dmrg_fci = abs(e_dmrg - e_fci) if (e_dmrg is not None and e_fci is not None) else None

        def col(x, w=12, p=6): return f"{x:{w}.{p}f}" if x is not None else f"{'--':>{w}}"
        c_df = f"{dmrg_fci:11.2e}" if dmrg_fci is not None else f"{'--':>11}"
        c_cb = f"{corr_beyond * 1e3:9.1f}m" if corr_beyond is not None else f"{'--':>10}"
        print(f"({nelec:>2},{norb:>2}) {nq:6d} {ndet:14,d} {col(e_kry)} {col(e_fci)} "
              f"{col(e_dmrg)} {c_df} {c_cb}")

        row = {
            "cas_electrons": nelec, "cas_orbitals": norb, "qubits": nq, "ndet": ndet,
            "krylov_energy": e_kry, "fci_energy": e_fci, "dmrg_energy": e_dmrg,
            "dmrg_vs_fci": dmrg_fci, "corr_beyond_cap_Ha": corr_beyond,
            "reference": "dmrg" if e_dmrg is not None else ("fci" if e_fci is not None else "none"),
        }
        with open(OUTPUT, "a", newline="") as fh:        # append so resumed rungs accumulate
            w = csv.DictWriter(fh, fieldnames=FIELDS)
            if not file_exists:
                w.writeheader()
                file_exists = True
            w.writerow(row)

    print("\n'corr>cap' = how much lower the active-space energy goes once the space grows past the")
    print(f"CAS(8,8)/{KRYLOV_QUBIT_CAP}-qubit ceiling the statevector quantum solver is limited to --")
    print(f"i.e. the correlation it structurally cannot reach.  ->  {OUTPUT}")


if __name__ == "__main__":
    main()
