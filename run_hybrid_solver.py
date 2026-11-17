#!/usr/bin/env python3
"""
Krylov Quantum Solver -- command-line entry point.

Pipeline:  CIF -> PySCF CASCI (classical)  ->  vetted Jordan-Wigner qubit Hamiltonian
           ->  real-time quantum Krylov subspace diagonalisation.

This replaces the previous entry point, which drove the broken
``EnterprisePipelineOrchestrator`` (incorrect mapping, near-identity Krylov basis,
asymmetric-Gaussian "noise", QCIVET stamp). See REFACTOR_PLAN.md.
"""
import argparse
import csv
import os
import time

from hybrid_quantum_solver.chemistry_gateway import load_and_compute_integrals
from hybrid_quantum_solver.pipeline import run_from_integrals


def main():
    parser = argparse.ArgumentParser(
        description="Run a hybrid quantum-classical ground-state estimate via real-time quantum Krylov.")
    parser.add_argument("--input_file", "-i", required=True,
                        help="Path to the molecular .cif structure file (e.g. data/.../NbN.cif).")
    parser.add_argument("--active_space", "-a", default="8,8",
                        help='CASCI active space as "electrons,orbitals" (default: 8,8).')
    parser.add_argument("--krylov_dim", "-d", type=int, default=8,
                        help="Krylov subspace dimension M (default: 8).")
    parser.add_argument("--shots", "-s", type=int, default=None,
                        help="If set, model finite-sampling shot noise (sigma ~ 1/sqrt(shots)). "
                             "Omit for the exact statevector result.")
    parser.add_argument("--seed", type=int, default=None, help="RNG seed for shot noise.")
    parser.add_argument("--output", "-o", default="krylov_telemetry.csv",
                        help="Output CSV filename for telemetry (default: krylov_telemetry.csv).")
    parser.add_argument("--uncertified-estimates-only", action="store_true",
                        help="Skip certification: no checker run, and the result is an UNVERIFIED "
                             "estimate. One checker run is about half an hour at 16 qubits.")
    args = parser.parse_args()

    if not os.path.isfile(args.input_file):
        parser.error(f"input_file does not exist: {args.input_file}")
    try:
        cas_elec, cas_orb = map(int, args.active_space.split(","))
    except ValueError:
        parser.error('active_space must be "electrons,orbitals", e.g. "8,8"')

    print("=" * 80)
    print("[INIT] Krylov Quantum Solver")
    print("=" * 80)
    print(f"-> Structure:     {os.path.basename(args.input_file)}")
    print(f"-> Active space:  CAS({cas_elec},{cas_orb})  ->  {cas_orb * 2} qubits")
    print(f"-> Krylov dim M:  {args.krylov_dim}")
    print(f"-> Shots:         {args.shots if args.shots else 'exact (noiseless)'}")

    print("\n[PHASE 1] Classical CASCI (PySCF)...")
    h1, eri, n_orb, casci_total, e_core, nelecas = load_and_compute_integrals(
        args.input_file, cas_electrons=cas_elec, cas_orbitals=cas_orb
    )
    print(f"   CASCI total energy (active-space FCI target): {casci_total:.6f} Ha "
          f"| nelecas={nelecas}")

    print("\n[PHASE 2] Quantum Krylov subspace diagonalisation...")
    if not args.uncertified_estimates_only:
        print("   Certifying the result out of process: seconds up to 12 qubits, about half an "
              "hour at 16.\n   --uncertified-estimates-only skips it and returns an "
              "UNVERIFIED estimate.")
    start = time.time()
    result = run_from_integrals(
        h1, eri, num_particles=nelecas, e_core=float(e_core),
        krylov_dim=args.krylov_dim, shots=args.shots, seed=args.seed,
        reference_energy=float(casci_total), track_convergence=True,
        uncertified_estimates_only=args.uncertified_estimates_only,
    )
    elapsed = time.time() - start

    # Certificate-or-abstention: the number is an energy only if the checker VERIFIED it.
    # Either way it is printed with the checker's own line, never as a bare result.
    estimate, certificate = result.unverified_energy, result.certificate
    error = abs(estimate - casci_total)
    label = "Ground-state energy:  " if certificate.ok else "UNVERIFIED estimate:  "
    print("\n" + "=" * 80)
    print(f"[RESULT] {label} {estimate:.6f} Ha")
    print(f"         Certificate:           {certificate.line}")
    print(f"         Hartree-Fock:          {result.hf_energy:.6f} Ha")
    print(f"         CASCI (active FCI):     {casci_total:.6f} Ha")
    print(f"         |error| vs CASCI:       {error:.3e} Ha")
    print(f"         qubits={result.n_qubits} rank={result.rank} dt={result.dt:.4f} "
          f"time={elapsed:.1f}s")
    print("=" * 80)

    fieldnames = ["molecule", "qubits", "krylov_dim", "rank", "shots", "computed_energy",
                  "hf_energy", "casci_total", "error_vs_casci", "time_s", "certificate"]
    file_exists = os.path.isfile(args.output)
    if file_exists:
        with open(args.output, newline="") as f:
            if next(csv.reader(f), None) != fieldnames:
                parser.error(f"{args.output} has other columns (it predates the certificate "
                             "column); pass another --output.")
    with open(args.output, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if not file_exists:
            writer.writeheader()
        writer.writerow({
            "molecule": os.path.splitext(os.path.basename(args.input_file))[0],
            "qubits": result.n_qubits, "krylov_dim": result.krylov_dim,
            "rank": result.rank, "shots": args.shots,
            "computed_energy": estimate, "hf_energy": result.hf_energy,
            "casci_total": casci_total, "error_vs_casci": error,
            "time_s": round(elapsed, 3),
            # VERIFIED, or why computed_energy is only an estimate.
            "certificate": "VERIFIED" if certificate.ok else certificate.line,
        })
    print(f"Telemetry appended to {args.output}")


if __name__ == "__main__":
    main()
