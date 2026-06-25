#!/usr/bin/env python3
"""
Parameter sweep harness for the validated quantum Krylov solver.

Runs in two stages so you can inspect Stage 1 before spending compute on Stage 2:

    python sweep_hybrid_solver.py 1     # active-space convergence
    # ... inspect data/sweep_stage1_active_space.csv, pick CHOSEN_ACTIVE_SPACE ...
    python sweep_hybrid_solver.py 2     # krylov_dim x shots grid

Stage 1 - active-space convergence. Varies (cas_elec, cas_orb) at a fixed krylov_dim with
          exact (noiseless) evaluation, looking for where the computed energy stops moving.
          This is the expensive axis (PySCF CASCI + the O(N^4) Jordan-Wigner build scale with
          cas_orb), so keep the list short and commit to ONE active space before Stage 2.

Stage 2 - krylov_dim x shots grid at the single active space chosen from Stage 1. Each finite
          shot count is repeated over several seeds, since one shot-noisy run is one draw.

This drives the validated path (hybrid_quantum_solver.molecular_hamiltonian +
quantum_krylov_solver), which reproduces FCI/CASCI to sub-mHa and is variationally bounded
(no estimate can fall below the true ground state). The old EnterprisePipelineOrchestrator /
qDRIFT / QCIVET core it used to call is retired (see REFACTOR_PLAN.md).

Caching mirrors the genuinely expensive operations:
  - integrals  cached per (cas_elec, cas_orb)   -> reused across every krylov_dim / shots / seed.
  - Hamiltonian cached per (cas_elec, cas_orb)   -> the Jordan-Wigner build runs once per active
                                                    space; only the (cheap) Krylov solve repeats.

SCIENTIFIC CAVEAT: for a crystalline CIF, load_and_compute_integrals builds a finite molecular
cluster from the unit-cell atoms with no periodic boundary conditions -- it is not the periodic
solid. Treat materials numbers accordingly (REFACTOR_PLAN.md, Phase 4).
"""

import concurrent.futures
import csv
import os
import sys
import time
import numpy as np

from hybrid_quantum_solver.chemistry_gateway import load_and_compute_integrals
from hybrid_quantum_solver.molecular_hamiltonian import build_hamiltonian_from_integrals
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from hybrid_quantum_solver.noise import shot_noise_sigma

INPUT_FILE = "data/nb_structures/Nb3I8_mp-27772_R-3m_stable.cif"

# ---- Stage 1: active-space convergence ------------------------------------
ACTIVE_SPACES = [(6, 6), (8, 8), (8, 10), (10, 10)]   # (electrons, orbitals)

# Stage 1 only needs to see where the energy plateau is as a function of active space.
# A large krylov_dim is irrelevant here; use the smallest M that still gives a converged
# Krylov estimate.  Stage 2 sweeps M properly.
STAGE1_KRYLOV_DIM = 8
STAGE1_OUTPUT = "data/sweep_stage1_active_space.csv"

# ---- Stage 2: krylov_dim x shots grid -------------------------------------
CHOSEN_ACTIVE_SPACE = (8, 8)          # <- set this after reviewing Stage 1
KRYLOV_DIMS = [2, 4, 6, 8, 12]
# None = exact statevector expectation (noiseless). Finite ints model shot noise (sigma ~ 1/sqrt(shots)).
SHOTS = [None, 100_000, 10_000, 1_000]
SEEDS_PER_SHOTS = 3                   # repeats for finite-shot points

STAGE2_OUTPUT = "data/sweep_stage2_krylov_shots.csv"

# Per-point wall-clock timeout in seconds.  Points that exceed this are written to CSV with
# audit_state="TIMEOUT" and skipped.  Set to None to disable.
MAX_RUNTIME_PER_POINT_S = 300


def _run_with_timeout(fn, timeout_s, *args, **kwargs):
    """Call fn(*args, **kwargs) in a daemon thread; raise TimeoutError past timeout_s seconds.

    The thread itself is NOT killed -- it keeps running until the process exits.  That is
    acceptable here because each point is independent and the next point uses cached inputs.
    """
    if timeout_s is None:
        return fn(*args, **kwargs)
    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
        future = pool.submit(fn, *args, **kwargs)
        try:
            return future.result(timeout=timeout_s)
        except concurrent.futures.TimeoutError:
            raise TimeoutError(
                f"Exceeded {timeout_s}s wall-clock limit. "
                f"Try lowering STAGE1_KRYLOV_DIM or MAX_RUNTIME_PER_POINT_S."
            )


def get_integrals(cas_elec, cas_orb, integral_cache):
    """CACHE LEVEL 1: PySCF CASCI extraction, keyed only by active space."""
    key = (cas_elec, cas_orb)
    if key not in integral_cache:
        integral_cache[key] = load_and_compute_integrals(
            INPUT_FILE, cas_electrons=cas_elec, cas_orbitals=cas_orb
        )
    return integral_cache[key]


def get_hamiltonian(cas_elec, cas_orb, integral_cache, ham_cache):
    """CACHE LEVEL 2: qubit Hamiltonian (Jordan-Wigner build), keyed by active space.

    Returns (MolecularHamiltonian, casci_total_energy).  casci_total is the active-space FCI
    target the Krylov estimate should approach -- and the energy frame is already correct
    because build_hamiltonian_from_integrals carries e_core as the energy_offset.
    """
    key = (cas_elec, cas_orb)
    if key not in ham_cache:
        h1, eri, _n_orb, casci_total, e_core, nelecas = get_integrals(
            cas_elec, cas_orb, integral_cache
        )
        mh = build_hamiltonian_from_integrals(
            h1, eri, num_particles=nelecas, energy_offset=float(e_core)
        )
        ham_cache[key] = (mh, casci_total)
    return ham_cache[key]


def run_one(cas_elec, cas_orb, krylov_dim, shots, seed, integral_cache, ham_cache):
    mh, casci_total = get_hamiltonian(cas_elec, cas_orb, integral_cache, ham_cache)

    sigma = shot_noise_sigma(shots) if shots else 0.0
    solver = QuantumKrylovSolver(mh, dt=None, noise_sigma=sigma, seed=seed)

    start = time.time()
    step = _run_with_timeout(solver.solve, MAX_RUNTIME_PER_POINT_S, krylov_dim)
    elapsed = time.time() - start

    computed_energy = step.energy
    error_vs_casci = abs(computed_energy - casci_total)

    return {
        "cas_elec": cas_elec,
        "cas_orb": cas_orb,
        "qubits": mh.num_qubits,
        "krylov_dim": krylov_dim,
        "rank": step.rank,
        "shots": shots,
        "seed": seed,
        "dt": round(solver.dt, 6),
        "computed_energy": computed_energy,
        "casci_total": casci_total,      # active-space FCI target (same frame)
        "error_vs_casci": error_vs_casci,
        "time_s": round(elapsed, 4),
        "audit_state": "OK",
    }


def write_row(path, row):
    """Append a single row immediately after each run so a crash mid-sweep keeps partials."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    file_exists = os.path.isfile(path)
    with open(path, mode="a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)


def _timeout_row(cas_elec, cas_orb, krylov_dim, shots, seed):
    return {
        "cas_elec": cas_elec, "cas_orb": cas_orb, "qubits": float("nan"),
        "krylov_dim": krylov_dim, "rank": float("nan"), "shots": shots, "seed": seed,
        "dt": float("nan"), "computed_energy": float("nan"), "casci_total": float("nan"),
        "error_vs_casci": float("nan"), "time_s": MAX_RUNTIME_PER_POINT_S,
        "audit_state": "TIMEOUT",
    }


def run_stage1():
    print(f"[Stage 1] Active-space convergence sweep on {INPUT_FILE}")
    integral_cache, ham_cache = {}, {}
    valid_rows = []
    for cas_elec, cas_orb in ACTIVE_SPACES:
        print(f"  -> CAS({cas_elec},{cas_orb})  [{cas_orb * 2} qubits]")
        print("=" * 80)
        try:
            row = run_one(cas_elec, cas_orb, STAGE1_KRYLOV_DIM, None, None,
                          integral_cache, ham_cache)
        except TimeoutError as exc:
            print(f"     TIMEOUT ({exc})")
            write_row(STAGE1_OUTPUT, _timeout_row(cas_elec, cas_orb, STAGE1_KRYLOV_DIM, None, None))
            continue
        except Exception as exc:
            print(f"     SKIPPED ({type(exc).__name__}: {exc})")
            continue

        print(f"     E = {row['computed_energy']:.6f} Ha  "
              f"CASCI = {row['casci_total']:.6f} Ha  "
              f"|err| = {row['error_vs_casci']:.2e} Ha  "
              f"rank = {row['rank']}  t = {row['time_s']}s")
        write_row(STAGE1_OUTPUT, row)
        valid_rows.append(row)

    if len(valid_rows) > 1:
        print("\n  Successive energy shifts (look for the plateau):")
        for prev, cur in zip(valid_rows, valid_rows[1:]):
            shift = abs(cur["computed_energy"] - prev["computed_energy"])
            print(f"    CAS({prev['cas_elec']},{prev['cas_orb']}) -> "
                  f"CAS({cur['cas_elec']},{cur['cas_orb']}): |dE| = {shift:.2e} Ha")
    else:
        print("\n  [WARN] Fewer than 2 results -- cannot identify a plateau.")

    print(f"\n[Stage 1] done -> {STAGE1_OUTPUT}")
    print("Set CHOSEN_ACTIVE_SPACE in this script, then run stage 2.")


def run_stage2():
    cas_elec, cas_orb = CHOSEN_ACTIVE_SPACE
    print(f"[Stage 2] krylov_dim x shots grid at CAS({cas_elec},{cas_orb}) on {INPUT_FILE}")
    integral_cache, ham_cache = {}, {}
    for krylov_dim in KRYLOV_DIMS:
        for shots in SHOTS:
            seeds = [None] if shots is None else range(SEEDS_PER_SHOTS)
            for seed in seeds:
                print(f"  -> M={krylov_dim:<3} shots={str(shots):<8} seed={seed}")
                try:
                    row = run_one(cas_elec, cas_orb, krylov_dim, shots, seed,
                                  integral_cache, ham_cache)
                except TimeoutError as exc:
                    print(f"     TIMEOUT ({exc})")
                    write_row(STAGE2_OUTPUT, _timeout_row(cas_elec, cas_orb, krylov_dim, shots, seed))
                    continue
                except Exception as exc:
                    print(f"     SKIPPED ({type(exc).__name__}: {exc})")
                    continue
                print(f"     E = {row['computed_energy']:.6f} Ha   "
                      f"|err| vs CASCI = {row['error_vs_casci']:.2e} Ha")
                write_row(STAGE2_OUTPUT, row)

    print(f"\n[Stage 2] done -> {STAGE2_OUTPUT}")


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "1"
    if stage == "1":
        run_stage1()
    elif stage == "2":
        run_stage2()
    else:
        print("Usage: python sweep_hybrid_solver.py [1|2]")
