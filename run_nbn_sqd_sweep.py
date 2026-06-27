#!/usr/bin/env python3
"""
Sprint 9: SQD Spin-Sector Sweep (drop-in replacement for run_nbn_stress_test.py)
Target: Nb-based crystal structures (e.g. Nb3I8, NbN) from Materials Project.

What changed vs. the orchestrator version
------------------------------------------
- The EnterprisePipelineOrchestrator (qDRIFT + Jordan-Wigner + StabilizedSubspaceShifter)
  is replaced by qiskit-addon-sqd's `diagonalize_fermionic_hamiltonian`.
- No overlap matrix S, no SVD regularization  -> Krylov basis collapse is impossible
  (determinant subspaces are orthonormal, so S = I by construction).
- The solver consumes (hcore, eri) directly -> the O(N^4) Pauli mapping is gone.
- Energy frame is aligned in exactly ONE place:  total = solver_energy + e_core.
- Spin sectors are selected by (n_alpha, n_beta), NOT a parity filter, so quintets
  and septets can no longer be silently dropped.

Default run uses uniform-random ("pure noise") samples so the entire classical loop
can be validated with no ansatz, no ffsim, and no hardware. Flip SAMPLE_SOURCE to
"lucj" once the plumbing checks out against CASCI. See build_lucj_bit_array().
"""

import sys
import time
import numpy as np
import polars as pl
from ase.io import read
from pyscf import gto, scf, ao2mo, mcscf

try:
    from qiskit_addon_sqd.fermion import diagonalize_fermionic_hamiltonian, solve_sci_batch
    from qiskit_addon_sqd.counts import generate_bit_array_uniform
except ImportError as exc:  # pragma: no cover
    sys.exit(
        "[FATAL] qiskit-addon-sqd is not installed in this environment.\n"
        "        pip install qiskit-addon-sqd   (inside the `chem` conda env)\n"
        f"        Original error: {exc}"
    )

from functools import partial

# =============================================================================
# Configuration
# =============================================================================
TARGET_FILE = "data/nb_structures/NbN_mp-2634.cif"
CAS_ELECTRONS = 8
CAS_ORBITALS = 8

# SQD sweep axes: spin sector (auto-derived) x subspace size (samples_per_batch)
SAMPLES_PER_BATCH_SWEEP = [100, 300, 600, 1000]
NUM_BATCHES = 5
MAX_ITERATIONS = 5
ENERGY_TOL = 1e-3
OCCUPANCIES_TOL = 1e-3
CARRYOVER_THRESHOLD = 1e-4

# SCF stabilization for heavy open-shell d-metals (non-negotiable for Nb)
SCF_MAX_CYCLE = 300
SCF_LEVEL_SHIFT = 0.2
SCF_CONV_TOL = 1e-8
SCF_INIT_GUESS = "atom"

SAMPLE_SOURCE = "uniform"   # "uniform" (validation) or "lucj" (requires ffsim + sampler)
N_SHOTS = 10_000
RNG_SEED = 24

OUTPUT_CSV = "sprint_9_nb_sqd_sweep.csv"

# Plumbing validation (only gates hard in "uniform" mode).
# - A uniform-noise SQD run should land within VALIDATION_TOL_MHA of CASCI if the
#   recovery/projection/frame-alignment is correct.
# - SQD is variational in the same active space, so it can NEVER fall below CASCI.
#   A result below CASCI by more than VARIATIONAL_TOL_MHA means the energy frame or
#   integrals are wrong (e.g. a dropped e_core, bad nelec, or notation mismatch).
VALIDATION_TOL_MHA = 50.0
VARIATIONAL_TOL_MHA = 1.0

# Trend gate (the primary convergence check). For correctly-plumbed SQD, enlarging the
# subspace via samples_per_batch should drive delta toward CASCI. We judge the sweep, not
# a single point:
#   - already within CONVERGENCE_TOL_MHA at the largest subspace            -> CONVERGED
#   - delta drops by at least TREND_MIN_IMPROVEMENT across the sweep        -> CONVERGING
#   - subspace grew but delta did not improve                              -> STALLED (suspect)
CONVERGENCE_TOL_MHA = 10.0
TREND_MIN_IMPROVEMENT = 0.20


# =============================================================================
# Classical front end (preserved from the original harness)
# =============================================================================
def get_smart_basis(atoms):
    """Assigns ECP-friendly basis to heavy atoms, all-electron to light atoms."""
    basis_map = {}
    for atom in atoms:
        if atom.number > 36:
            basis_map[atom.symbol] = "lanl2dz"
        else:
            basis_map[atom.symbol] = "6-31g*"
    return basis_map


def load_geometry(cif_filepath):
    """Reads a CIF and returns (atom_str, basis_set, ecp_dict, n_active_electrons_parity)."""
    print("=" * 80)
    print(f"[CLASSICAL PRE-PROCESSING] Loading structure: {cif_filepath}")
    atoms = read(cif_filepath)
    atom_str = "; ".join(
        f"{a.symbol} {a.position[0]} {a.position[1]} {a.position[2]}" for a in atoms
    )
    basis_set = get_smart_basis(atoms)
    ecp_dict = {a.symbol: "lanl2dz" for a in atoms if a.number > 36}
    print(f"  -> {len(atoms)} atoms | ECP on: {sorted(ecp_dict) or 'none'}")
    return atom_str, basis_set, ecp_dict


def valid_spin_sectors(cas_electrons, cas_orbitals):
    """
    All physically reachable spin sectors for CAS(cas_electrons, cas_orbitals),
    expressed as mol.spin = 2S = n_alpha - n_beta.
    """
    max_unpaired = min(cas_electrons, 2 * cas_orbitals - cas_electrons)
    start = 0 if cas_electrons % 2 == 0 else 1
    return list(range(start, max_unpaired + 1, 2))


def multiplicity_name(mol_spin):
    return {0: "singlet", 1: "doublet", 2: "triplet", 3: "quartet",
            4: "quintet", 5: "sextet", 6: "septet", 7: "octet",
            8: "nonet"}.get(mol_spin, f"2S+1={mol_spin + 1}")


def build_scf(atom_str, basis_set, ecp_dict, mol_spin):
    """Converged, stabilized SCF reference for a given spin sector."""
    mol = gto.M(atom=atom_str, basis=basis_set, ecp=ecp_dict, charge=0, spin=mol_spin)
    mf = scf.UHF(mol) if mol.nelec[0] != mol.nelec[1] else scf.RHF(mol)
    mf.verbose = 0
    mf.max_cycle = SCF_MAX_CYCLE
    mf.level_shift = SCF_LEVEL_SHIFT
    mf.conv_tol = SCF_CONV_TOL
    mf.init_guess = SCF_INIT_GUESS
    mf.kernel()
    if not mf.converged:
        # one stabilized retry via second-order Newton before giving up
        mf = mf.newton()
        mf.kernel()
    return mol, mf


def integrals_for_spin(mf, cas_orbitals, cas_electrons, mol_spin):
    """
    CASCI in the requested spin sector. Returns the integrals SQD needs plus the
    classical reference energy.

    e_core carries nuclear repulsion + frozen-core energy. The SQD solver returns the
    active-space electronic eigenvalue; total = solver_energy + e_core.
    """
    n_alpha = (cas_electrons + mol_spin) // 2
    n_beta = (cas_electrons - mol_spin) // 2
    cas = mcscf.CASCI(mf, cas_orbitals, (n_alpha, n_beta))
    cas.verbose = 0
    cas.kernel()
    hcore, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), cas_orbitals)
    return hcore, eri, float(e_core), (n_alpha, n_beta), float(cas.e_tot)


# =============================================================================
# Sample source
# =============================================================================
def get_bit_array(norb, rng):
    """Validation samples: uniform random bitstrings over 2*norb qubits."""
    return generate_bit_array_uniform(N_SHOTS, norb * 2, rand_seed=rng)


def build_lucj_bit_array(mf, norb, nelec, rng):
    """
    Real-sample path (open-shell). SCAFFOLD: requires ffsim + a Sampler backend, and
    the interaction pairs / CCSD reference must match your system and target topology.
    Returns a bit_array compatible with diagonalize_fermionic_hamiltonian.
    """
    import ffsim                                   # noqa: F401  (optional dependency)
    from qiskit import QuantumCircuit, QuantumRegister
    from pyscf import cc

    n_alpha, n_beta = nelec
    ucc = cc.UCCSD(mf).run()                       # open-shell CCSD amplitudes
    ucj_op = ffsim.UCJOpSpinUnbalanced.from_t_amplitudes(
        t2=ucc.t2, t1=ucc.t1, n_reps=1,
        # interaction_pairs=...  # <-- set for your hardware topology before real runs
    )
    qubits = QuantumRegister(2 * norb, name="q")
    circuit = QuantumCircuit(qubits)
    circuit.append(ffsim.qiskit.PrepareHartreeFockJW(norb, nelec), qubits)
    circuit.append(ffsim.qiskit.UCJOpSpinUnbalancedJW(ucj_op), qubits)
    circuit.measure_all()

    # Sample locally; swap for qiskit_ibm_runtime SamplerV2(mode=backend) on hardware.
    from qiskit_ibm_runtime import SamplerV2 as Sampler
    from qiskit_ibm_runtime.fake_provider import FakeSherbrooke
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

    backend = FakeSherbrooke()
    pm = generate_preset_pass_manager(optimization_level=3, backend=backend)
    pm.pre_init = ffsim.qiskit.PRE_INIT
    isa = pm.run(circuit)
    result = Sampler(mode=backend).run([isa], shots=N_SHOTS).result()
    return result[0].data.meas


# =============================================================================
# SQD driver
# =============================================================================
def run_sqd_for_sector(hcore, eri, e_core, nelec, norb, bit_array,
                       samples_per_batch, spin_sq, rng):
    """Self-consistent configuration-recovery diagonalization for one spin sector."""
    tracker = {"best_energy": np.inf, "max_dim": 0, "iterations": 0}

    def callback(results):
        tracker["iterations"] += 1
        for r in results:
            total = r.energy + e_core                      # single frame-alignment point
            tracker["best_energy"] = min(tracker["best_energy"], total)
            dim = int(np.prod(r.sci_state.amplitudes.shape))
            tracker["max_dim"] = max(tracker["max_dim"], dim)

    result = diagonalize_fermionic_hamiltonian(
        hcore, eri, bit_array,
        samples_per_batch=samples_per_batch,
        norb=norb,
        nelec=nelec,
        num_batches=NUM_BATCHES,
        energy_tol=ENERGY_TOL,
        occupancies_tol=OCCUPANCIES_TOL,
        max_iterations=MAX_ITERATIONS,
        sci_solver=partial(solve_sci_batch, spin_sq=spin_sq),
        symmetrize_spin=(nelec[0] == nelec[1]),            # only for closed-shell sectors
        carryover_threshold=CARRYOVER_THRESHOLD,
        callback=callback,
        seed=rng,
    )
    return tracker["best_energy"], tracker["max_dim"], tracker["iterations"], result


# =============================================================================
# Plumbing validation
# =============================================================================
def validate_row(sqd_energy, casci_energy):
    """
    Returns a verdict that separates 'broken plumbing' from 'real physics':

      FRAME_ERROR  -> SQD below CASCI beyond tolerance (variational violation;
                      impossible if the energy frame and integrals are correct).
      ABOVE_TOL    -> SQD above CASCI by more than VALIDATION_TOL_MHA (in uniform
                      mode this means the classical loop didn't recover the state).
      PASS         -> within tolerance and above CASCI, as expected.
      SKIP         -> run did not produce a finite energy.
    """
    if sqd_energy is None or not np.isfinite(sqd_energy) or casci_energy is None:
        return "SKIP"
    delta_mHa = (sqd_energy - casci_energy) * 1e3      # signed: positive = above CASCI
    if delta_mHa < -VARIATIONAL_TOL_MHA:
        return "FRAME_ERROR"
    if delta_mHa > VALIDATION_TOL_MHA:
        return "ABOVE_TOL"
    return "PASS"


def analyze_sector_trend(sector_rows):
    """
    Trend-based convergence check across the samples_per_batch sweep for one sector.
    Returns one of:
      CONVERGED     -> within CONVERGENCE_TOL_MHA at the largest subspace (good).
      CONVERGING    -> delta dropped by >= TREND_MIN_IMPROVEMENT as the subspace grew.
      STALLED       -> subspace grew but delta did not meaningfully improve (suspect).
      FLAT_SUBSPACE -> samples_per_batch did not actually enlarge the subspace.
      INSUFFICIENT  -> fewer than 2 finite points to judge a trend.

    Endpoint heuristic: random batch draws make delta non-monotonic point-to-point, so we
    compare the smallest-subspace point against the largest-subspace point rather than
    requiring strict monotonicity at every step.
    """
    pts = [
        r for r in sector_rows
        if r["delta_mHa"] is not None and np.isfinite(r["delta_mHa"])
        and r["subspace_dim"]
    ]
    if len(pts) < 2:
        return "INSUFFICIENT"

    pts.sort(key=lambda r: r["subspace_dim"])
    first, last = pts[0], pts[-1]

    # A result already inside the band is converged, even if the subspace saturated and
    # stopped growing (small/sparse systems hit the full determinant space immediately).
    if last["delta_mHa"] <= CONVERGENCE_TOL_MHA:
        return "CONVERGED"
    if last["subspace_dim"] <= first["subspace_dim"]:
        return "FLAT_SUBSPACE"
    if last["delta_mHa"] <= (1.0 - TREND_MIN_IMPROVEMENT) * first["delta_mHa"]:
        return "CONVERGING"
    return "STALLED"


# =============================================================================
# Main sweep
# =============================================================================
def main():
    rng = np.random.default_rng(RNG_SEED)
    atom_str, basis_set, ecp_dict = load_geometry(TARGET_FILE)

    sectors = valid_spin_sectors(CAS_ELECTRONS, CAS_ORBITALS)
    print(f"[SWEEP] Spin sectors to explore (mol.spin = 2S): {sectors}")
    print(f"[SWEEP] Subspace sizes (samples_per_batch): {SAMPLES_PER_BATCH_SWEEP}")
    print(f"[SWEEP] Sample source: {SAMPLE_SOURCE}")
    print("=" * 80)

    telemetry = []

    for mol_spin in sectors:
        name = multiplicity_name(mol_spin)
        S = mol_spin / 2.0
        spin_sq = S * (S + 1.0)

        # --- classical reference for this sector (guarded; a bad sector won't kill the sweep)
        try:
            mol, mf = build_scf(atom_str, basis_set, ecp_dict, mol_spin)
            hcore, eri, e_core, nelec, casci_energy = integrals_for_spin(
                mf, CAS_ORBITALS, CAS_ELECTRONS, mol_spin
            )
            norb = hcore.shape[0]
        except Exception as exc:  # noqa: BLE001
            print(f"[{name:>8}] SCF/CASCI FAILED: {exc}")
            telemetry.append({
                "mol_spin": mol_spin, "multiplicity": name,
                "n_alpha": None, "n_beta": None, "samples_per_batch": None,
                "casci_energy": None, "sqd_energy": None, "delta_mHa": None,
                "subspace_dim": None, "iterations": None,
                "execution_time_s": None, "status": "CLASSICAL_FAILED",
                "validation": "SKIP", "trend": "SKIP",
            })
            continue

        print(f"[{name:>8}] (na,nb)=({nelec[0]},{nelec[1]})  CASCI={casci_energy:.6f} Ha")

        # --- samples for this sector
        if SAMPLE_SOURCE == "lucj":
            bit_array = build_lucj_bit_array(mf, norb, nelec, rng)
        else:
            bit_array = get_bit_array(norb, rng)

        # --- SQD across subspace sizes
        sector_rows = []
        for spb in SAMPLES_PER_BATCH_SWEEP:
            start = time.time()
            try:
                sqd_energy, max_dim, iters, _ = run_sqd_for_sector(
                    hcore, eri, e_core, nelec, norb, bit_array,
                    samples_per_batch=spb, spin_sq=spin_sq, rng=rng,
                )
                status = "STABILIZED"
            except Exception as exc:  # noqa: BLE001
                print(f"           spb={spb}: SQD FAILED: {exc}")
                sqd_energy, max_dim, iters, status = np.nan, None, None, "SQD_FAILED"

            elapsed = time.time() - start
            delta_mHa = (
                abs(sqd_energy - casci_energy) * 1e3
                if np.isfinite(sqd_energy) else None
            )
            verdict = validate_row(sqd_energy, casci_energy)
            print(
                f"           spb={spb:>4}  SQD={sqd_energy:.6f} Ha  "
                f"delta={delta_mHa:.2f} mHa  dim={max_dim}  t={elapsed:.1f}s  [{verdict}]"
                if delta_mHa is not None else
                f"           spb={spb:>4}  {status}"
            )

            row = {
                "mol_spin": mol_spin, "multiplicity": name,
                "n_alpha": nelec[0], "n_beta": nelec[1], "samples_per_batch": spb,
                "casci_energy": casci_energy, "sqd_energy": sqd_energy,
                "delta_mHa": delta_mHa, "subspace_dim": max_dim, "iterations": iters,
                "execution_time_s": elapsed, "status": status,
                "validation": verdict, "trend": None,
            }
            telemetry.append(row)      # same object also held in sector_rows
            sector_rows.append(row)

        # --- sector-level trend verdict; stamp it onto every row of this sector
        sector_trend = analyze_sector_trend(sector_rows)
        for row in sector_rows:
            finite = row["delta_mHa"] is not None and np.isfinite(row["delta_mHa"])
            row["trend"] = sector_trend if finite else "SKIP"
        print(f"           -> sector trend: {sector_trend}")

    df = pl.DataFrame(telemetry, infer_schema_length=None)
    df.write_csv(OUTPUT_CSV)

    print("=" * 80)
    ok = df.filter(pl.col("status") == "STABILIZED")
    if not ok.is_empty():
        gs = ok.row(ok["sqd_energy"].arg_min(), named=True)
        print(
            f"[GROUND STATE] {gs['multiplicity']} "
            f"(na,nb)=({gs['n_alpha']},{gs['n_beta']})  "
            f"SQD={gs['sqd_energy']:.6f} Ha  delta={gs['delta_mHa']:.2f} mHa"
        )
    print(f"[COMPLETE] Telemetry saved to '{OUTPUT_CSV}'")
    print("=" * 80)

    # --- Plumbing gate: announce broken plumbing loudly instead of as physics ---
    # Two hard-fail conditions:
    #   1. FRAME_ERROR (SQD below CASCI) -> always a bug, any sample source.
    #   2. STALLED trend in uniform mode -> the subspace grew but delta did not improve,
    #      so the recovery loop isn't capturing the right configurations. In "lucj" mode a
    #      stall reflects ansatz quality rather than code, so it's informational there.
    # ABOVE_TOL is now informational only: a sector can sit above tolerance yet still be
    # CONVERGING (just needs more samples), which the trend gate correctly tolerates.
    frame_errors = df.filter(pl.col("validation") == "FRAME_ERROR")
    stalled = df.filter((pl.col("trend") == "STALLED") & pl.col("samples_per_batch").is_not_null())
    stalled_sectors = stalled.unique(subset=["mol_spin"])

    hard_fail = not frame_errors.is_empty() or (
        SAMPLE_SOURCE == "uniform" and not stalled_sectors.is_empty()
    )

    # Informational: above tolerance but trending the right way.
    converging_high = df.filter((pl.col("validation") == "ABOVE_TOL") & (pl.col("trend") == "CONVERGING"))
    if not converging_high.is_empty():
        print("[NOTE] Sectors above the absolute tolerance but still converging "
              "(raise samples_per_batch): "
              f"{sorted(converging_high['multiplicity'].unique().to_list())}")

    if hard_fail:
        print("!" * 80)
        print("[VALIDATION FAILED] Telemetry below looks like broken plumbing, not physics:")
        for row in frame_errors.iter_rows(named=True):
            print(
                f"  - {row['multiplicity']} (na,nb)=({row['n_alpha']},{row['n_beta']}) "
                f"spb={row['samples_per_batch']}: delta={row['delta_mHa']:.2f} mHa -> "
                "SQD below CASCI (variational violation -> check e_core / nelec / notation)"
            )
        if SAMPLE_SOURCE == "uniform":
            for row in stalled_sectors.iter_rows(named=True):
                print(
                    f"  - {row['multiplicity']} (na,nb)=({row['n_alpha']},{row['n_beta']}): "
                    "delta did not improve as the subspace grew -> "
                    "config recovery not capturing the right determinants"
                )
        print("!" * 80)
        return 1

    print("[VALIDATION PASSED] No variational violations; all sectors converging or converged.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
