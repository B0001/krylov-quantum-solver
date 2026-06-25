#!/usr/bin/env python3
"""
Materials driver: CIF -> PySCF CASCI -> validated quantum Krylov solver.

    python run_real_material_simulation.py <path_to_cif> [--active_space E,O] [--krylov_dim M] [--shots N]

This is the CIF/materials front end for the validated pipeline. It reuses the same classical
pre-processing as the main CLI (hybrid_quantum_solver.chemistry_gateway.load_and_compute_integrals:
ground-state spin scan, smart ECP basis, CASCI active-space truncation) and feeds the active-space
integrals into hybrid_quantum_solver.pipeline.run_from_integrals.

The retired path (EnterprisePipelineOrchestrator with its qDRIFT/QCIVET core, and this file's old
AO-basis-integral shortcut) produced energies hundreds of Ha below the true ground state. See
REFACTOR_PLAN.md.

SCIENTIFIC CAVEAT: a crystalline CIF is treated as a FINITE molecular cluster built from the
unit-cell atoms, with no periodic boundary conditions. This is not a calculation of the periodic
solid; treat transition-metal/materials numbers as a research probe, not a validated result
(REFACTOR_PLAN.md, Phase 4).
"""

import argparse

from hybrid_quantum_solver.chemistry_gateway import load_and_compute_integrals
from hybrid_quantum_solver.pipeline import run_from_integrals


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("cif", help="Path to the .cif structure file.")
    parser.add_argument("--active_space", "-a", default="8,8",
                        help='CASCI active space as "electrons,orbitals" (default: 8,8).')
    parser.add_argument("--krylov_dim", "-d", type=int, default=8,
                        help="Krylov subspace dimension M (default: 8).")
    parser.add_argument("--shots", "-s", type=int, default=None,
                        help="Model finite-sampling shot noise (omit for exact statevector).")
    parser.add_argument("--seed", type=int, default=None, help="RNG seed for shot noise.")
    args = parser.parse_args()

    try:
        cas_elec, cas_orb = map(int, args.active_space.split(","))
    except ValueError:
        parser.error('--active_space must be "electrons,orbitals", e.g. "8,8"')

    print("=" * 80)
    print(f"[MATERIALS] {args.cif}  |  CAS({cas_elec},{cas_orb}) -> {cas_orb * 2} qubits  |  M={args.krylov_dim}")
    print("=" * 80)

    h1, eri, _n_orb, casci_total, e_core, nelecas = load_and_compute_integrals(
        args.cif, cas_electrons=cas_elec, cas_orbitals=cas_orb
    )
    print(f"[CLASSICAL] CASCI total (active-space FCI target): {casci_total:.6f} Ha | nelecas={nelecas}")

    result = run_from_integrals(
        h1, eri, num_particles=nelecas, e_core=float(e_core),
        krylov_dim=args.krylov_dim, shots=args.shots, seed=args.seed,
        reference_energy=float(casci_total), track_convergence=True,
    )

    print("\n" + "=" * 80)
    print(f"[RESULT] Ground-state estimate: {result.computed_energy:.6f} Ha")
    print(f"         Hartree-Fock:          {result.hf_energy:.6f} Ha")
    print(f"         CASCI (active FCI):     {casci_total:.6f} Ha")
    print(f"         |error| vs CASCI:       {result.error_vs_reference:.3e} Ha")
    print(f"         qubits={result.n_qubits} rank={result.rank} dt={result.dt:.4f} shots={args.shots}")
    print("=" * 80)
    print("NOTE: finite-cluster approximation of a periodic crystal (see module docstring).")


if __name__ == "__main__":
    main()
