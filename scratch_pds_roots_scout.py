#!/usr/bin/env python3
"""Scout probe for chem-pc1: higher roots of P_K(E) as excited-state estimates.

Not part of the test suite -- exploratory only, to pre-register the gate thresholds in
specs/SPEC_pds_excited_roots.md before writing tests/test_pds_excited_roots_spec.py.
"""
import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from moment_expansion import hamiltonian_moments, pds_roots


def reachable_spectrum(mh, overlap_tol=1e-8):
    w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    hf = np.asarray(mh.hf_state().data, dtype=complex)
    overlaps = np.abs(V.conj().T @ hf) ** 2
    reachable = w[overlaps > overlap_tol].real + mh.energy_offset
    return np.sort(reachable)


print("=== H4 equilibrium-ish (1,2,3 A spacing, as in SPEC_moment_pds) ===")
mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0")
reachable = reachable_spectrum(mh)
print("reachable spectrum (first 6):", reachable[:6])
mu, off = hamiltonian_moments(mh, 16)
for K in range(2, 9):
    try:
        roots = pds_roots(mu, K, off)
    except ValueError as e:
        print(f"K={K}: FAILED {e}")
        continue
    n = min(len(roots), len(reachable))
    viol = [(j, roots[j], reachable[j]) for j in range(n) if roots[j] < reachable[j] - 1e-9]
    err0 = (roots[0] - reachable[0]) * 1e3
    err1 = (roots[1] - reachable[1]) * 1e3 if len(roots) > 1 else float("nan")
    print(f"K={K}: n_real_roots={len(roots)} err_root0={err0:+.4f} mHa err_root1={err1:+.4f} mHa "
          f"violations={viol}")

print()
print("=== H4 stretched / near-degenerate (2,4,6 A spacing) ===")
mh2 = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 2.0; H 0 0 4.0; H 0 0 6.0")
reachable2 = reachable_spectrum(mh2)
print("reachable spectrum (first 8):", reachable2[:8])
mu2, off2 = hamiltonian_moments(mh2, 16)
for K in range(2, 9):
    try:
        roots = pds_roots(mu2, K, off2)
    except ValueError as e:
        print(f"K={K}: FAILED {e}")
        continue
    n = min(len(roots), len(reachable2))
    viol = [(j, roots[j], reachable2[j]) for j in range(n) if roots[j] < reachable2[j] - 1e-9]
    print(f"K={K}: n_real_roots={len(roots)} roots={np.round(roots[:6],6)} "
          f"reach={np.round(reachable2[:6],6)} violations={viol}")

print()
print("=== index-matching diagnostic (stretched H4, nearest-neighbor assignment) ===")
for K in range(4, 9):
    roots = pds_roots(mu2, K, off2)
    n_show = min(len(roots), 6)
    assigned = []
    for j in range(n_show):
        nearest = np.argmin(np.abs(reachable2[:8] - roots[j]))
        assigned.append(nearest)
    skipped = sorted(set(range(n_show)) - set(assigned))
    dup = len(assigned) != len(set(assigned))
    print(f"K={K}: roots[:6]={np.round(roots[:n_show],4)}")
    print(f"       nearest-ref-index per root = {assigned}  duplicate_match={dup}  "
          f"skipped_ref_indices(<{n_show})={skipped}")

print()
print("=== conditioning of M at high K (equilibrium H4) ===")
for K in range(2, 9):
    Mmat = np.array([[mu[2 * K - i - j] for j in range(1, K + 1)] for i in range(1, K + 1)])
    print(f"K={K}: cond(M)={np.linalg.cond(Mmat):.3e}")
