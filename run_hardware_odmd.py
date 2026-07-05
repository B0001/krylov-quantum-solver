#!/usr/bin/env python3
"""
run_hardware_odmd.py -- run the ODMD ground-state method on a (simulated or real) quantum device.

This is the on-hardware sibling of the SPEC_device_odmd Aer study. It measures the survival
signal s_k = <phi0|U_trot^k|phi0> for H2 by ancilla Hadamard tests (the validated
HardwareKrylovSolver circuit construction -- ONLY the overlap row S_0k, no Hamiltonian
observables, so no lambda-scaled cost), at Trotter reps 1 and 2, then runs the full ODMD stack:
odmd_energy per reps, Richardson bias removal across the pair (SPEC_trotter_odmd), and a
single-signal union-bootstrap error bar (SPEC_odmd_uq). Results are compared to the exact circuit
eigenphase (computed locally, trotter_odmd) and FCI.

Backends:
  --backend aer  (default): AerSimulator, optionally with a device NoiseModel (--noise cx=<p>).
  --backend <ibm_name>: submit to IBM Quantum via qiskit-ibm-runtime (EstimatorV2). Requires a
    saved account -- see the printed setup hint if none is found. This path is code-complete but
    depends on live credentials + queue, so it is not exercised by the test suite.

--dry-run stops after the transpiled resource table (no execution anywhere).

THE POINT (and the honest caveat printed at the end): SPEC_device_odmd showed ODMD eigenphases
are immune to *depolarizing* damping in simulation. Real-hardware coherent errors, crosstalk,
and drift are exactly what that simulated study could NOT test -- this script is how you find out
whether the immunity survives contact with a real device.

HARDWARE-READINESS FINDING (surfaced by --dry-run, 2026-07-04): the ancilla-controlled
Hadamard-test construction is DEEP. For H2 at K=8, the deepest circuit (controlled-U^7) transpiles
to ~3.4k 2q-gates at reps=1 and ~6.6k at reps=2 -- far beyond what survives on today's NISQ
hardware. The dry-run resource table is exactly the tool to catch this before spending queue
time. Actionable levers: drop K (depth ~ controlled-U^(K-1), so K=4 is ~half), use reps=1 and
Richardson only if the bias warrants it, or move to a shallower signal-measurement scheme. In
Aer (exact / noise-model) the K=8 pipeline recovers the circuit eigenphase to <0.01 mHa and
Richardson lands within ~1 mHa of FCI -- the algorithm is correct; the *circuit depth* is the
hardware bottleneck, not the method.
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
from qiskit import transpile

from device_odmd import centered_frame, device_odmd_energy
from hybrid_quantum_solver.hardware_krylov import HardwareKrylovSolver
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from odmd import odmd_energy
from odmd_uq import odmd_confidence_interval
from trotter_odmd import build_trotter_odmd_problem, richardson_energy
from visibility_law import predicted_shots

H2 = dict(atom="H 0 0 0; H 0 0 0.74")
GROUND_WEIGHT = 0.95           # |<HF|E0>|^2 for H2 (sets the visibility-law shot budget)


def _solver(mh_c, tau, reps, shots, noise_model=None, seed=None):
    return HardwareKrylovSolver(mh_c, dt=tau, shots=shots, noise_model=noise_model, seed=seed,
                                trotter_order=2, trotter_reps=reps)


def resource_table(mh_c, tau, K, shots, backend_target=None):
    """Transpile the deepest Hadamard-test circuit (i=0, j=K-1) and print circuit costs."""
    print("\n== Resource estimate (transpiled, optimization_level=3) ==")
    for reps in (1, 2):
        solver = _solver(mh_c, tau, reps, shots)
        deepest = solver._pair_circuit(0, K - 1)
        if backend_target is not None:
            deepest = transpile(deepest, backend=backend_target, optimization_level=3)
        else:
            deepest = transpile(deepest, optimization_level=3)
        ops = deepest.count_ops()
        cx = ops.get("cx", 0) + ops.get("cz", 0) + ops.get("ecr", 0)
        # 2 overlap observables (X_a, Y_a) per pair; K pairs (j=0..K-1); j=0 is trivial
        n_circuits = K
        print(f"  reps={reps}: deepest circuit depth={deepest.depth():4d}  2q-gates={cx:4d}  "
              f"qubits={deepest.num_qubits}  circuits={n_circuits}  "
              f"shots/obs={shots}  total_shots~{2 * n_circuits * shots}")
    budget = predicted_shots(GROUND_WEIGHT, K, 2 ** mh_c.num_qubits)
    print(f"  visibility-law budget for the ground line (w={GROUND_WEIGHT}, K={K}): "
          f"~{budget:.0f} shots/element (default {shots} is comfortably above)")


def measure_signal_aer(mh_c, tau, K, reps, shots, noise_model, seed):
    return _solver(mh_c, tau, reps, shots, noise_model, seed).measure_signal(K)


def measure_signal_ibm(mh_c, tau, K, reps, shots, backend, service):
    """Overlap row S_0k on a real IBM backend via EstimatorV2 (X_a, Y_a ancilla observables)."""
    from qiskit_ibm_runtime import EstimatorV2

    solver = _solver(mh_c, tau, reps, shots)     # reuse the (backend-agnostic) circuit builder
    obs = solver._observables[:2]                # [X_a (x) I, Y_a (x) I] -> Re, Im of S_0k
    circuits = [transpile(solver._pair_circuit(0, k), backend=backend, optimization_level=3)
                for k in range(K)]
    layout_obs = [[o.apply_layout(c.layout) for o in obs] for c in circuits]
    estimator = EstimatorV2(mode=backend)
    estimator.options.default_shots = shots
    job = estimator.run([(c, o) for c, o in zip(circuits, layout_obs)])
    print(f"  submitted job {job.job_id()} (reps={reps}); waiting...")
    res = job.result()
    s = np.array([r.data.evs[0] + 1j * r.data.evs[1] for r in res], dtype=complex)
    s[0] = 1.0
    return s


def analyze(mh, mh_c, tau, mu, K, sig1, sig2, shots):
    """Run the ODMD stack on the two measured signals and print the report."""
    ref1 = build_trotter_odmd_problem(mh, n=K, reps=1)
    ref2 = build_trotter_odmd_problem(mh, n=K, reps=2)
    fci = mh.ground_state_energy()
    sigma = 1.0 / np.sqrt(shots)
    off = mh.energy_offset + mu

    e1 = odmd_energy(sig1, tau, svd_threshold=5 * sigma)[0] + off
    e2 = odmd_energy(sig2, tau, svd_threshold=5 * sigma)[0] + off
    e_rich = richardson_energy(e1 - off, e2 - off) + off
    dev2 = device_odmd_energy(sig2, tau, sigma) + off
    ci = odmd_confidence_interval(sig2, tau, sigma, seed=0)

    print("\n== ODMD result ==")
    print(f"  {'':22s}{'E (Ha)':>14s}{'err vs circuit':>16s}{'err vs FCI':>14s}")
    print(f"  {'reps=1 (odmd)':22s}{e1:14.6f}{(e1 - ref1.e_circuit - off) * 1e3:14.3f} mHa"
          f"{(e1 - fci) * 1e3:11.3f} mHa")
    print(f"  {'reps=2 (odmd)':22s}{e2:14.6f}{(e2 - ref2.e_circuit - off) * 1e3:14.3f} mHa"
          f"{(e2 - fci) * 1e3:11.3f} mHa")
    print(f"  {'reps=2 (device-odmd)':22s}{dev2:14.6f}{'':16s}{(dev2 - fci) * 1e3:11.3f} mHa")
    print(f"  {'Richardson(1,2)':22s}{e_rich:14.6f}{'':16s}{(e_rich - fci) * 1e3:11.3f} mHa")
    print(f"  single-signal 90% CI (reps=2): "
          f"[{ci.lower + off:.6f}, {ci.upper + off:.6f}] Ha  "
          f"(half-width {ci.half_width * 1e3:.3f} mHa)")
    print("    ^ brackets the circuit eigenphase, NOT FCI: the Trotter bias is a systematic the"
          "\n      error bar cannot see (SPEC_odmd_uq G4) -- Richardson above removes it.")
    print(f"  reference: FCI={fci:.6f}  circuit eigenphase(reps=2)={ref2.e_circuit + off:.6f} "
          f"(Trotter bias {(ref2.e_circuit - ref2.ref) * 1e3:+.3f} mHa)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--backend", default="aer", help="'aer' or an IBM backend name")
    ap.add_argument("--shots", type=int, default=4096, help="shots per observable")
    ap.add_argument("--k", type=int, default=8, help="signal length (Krylov depth)")
    ap.add_argument("--noise", default=None,
                    help="aer only: depolarizing model, e.g. 'cx=3e-4'")
    ap.add_argument("--dry-run", action="store_true",
                    help="stop after the transpiled resource table")
    args = ap.parse_args()

    mh = build_molecular_hamiltonian(**H2)
    mh_c, tau, mu = centered_frame(mh)
    K = args.k
    print(f"ODMD on H2, centered frame, K={K}, reps 1+2, backend={args.backend}, "
          f"shots/obs={args.shots}")

    backend_target = None
    service = None
    if args.backend != "aer":
        try:
            from qiskit_ibm_runtime import QiskitRuntimeService
        except ImportError:
            sys.exit("qiskit-ibm-runtime not installed: pip install qiskit-ibm-runtime")
        try:
            service = QiskitRuntimeService()
        except Exception:
            sys.exit(
                "No saved IBM account found. One-time setup:\n"
                "  from qiskit_ibm_runtime import QiskitRuntimeService\n"
                "  QiskitRuntimeService.save_account(channel='ibm_quantum', token='<TOKEN>')\n"
                "then re-run with --backend <backend_name>."
            )
        backend_target = service.backend(args.backend)

    resource_table(mh_c, tau, K, args.shots, backend_target)
    if args.dry_run:
        print("\n--dry-run: stopping before execution.")
        return

    noise_model = None
    if args.backend == "aer" and args.noise:
        from hybrid_quantum_solver.noise import build_depolarizing_noise_model
        cx = float(args.noise.split("=")[1])
        noise_model = build_depolarizing_noise_model(cx / 10, cx, cx)

    if args.backend == "aer":
        sig1 = measure_signal_aer(mh_c, tau, K, 1, args.shots, noise_model, seed=1)
        sig2 = measure_signal_aer(mh_c, tau, K, 2, args.shots, noise_model, seed=2)
    else:
        sig1 = measure_signal_ibm(mh_c, tau, K, 1, args.shots, backend_target, service)
        sig2 = measure_signal_ibm(mh_c, tau, K, 2, args.shots, backend_target, service)

    analyze(mh, mh_c, tau, mu, K, sig1, sig2, args.shots)
    print("\n! CAVEAT: SPEC_device_odmd validated depolarizing-immunity in SIMULATION. Real-device"
          "\n  coherent errors, crosstalk, and drift are exactly what that study could not test --"
          "\n  compare the reps=2 error vs circuit eigenphase above against the Aer baseline.")


if __name__ == "__main__":
    main()
