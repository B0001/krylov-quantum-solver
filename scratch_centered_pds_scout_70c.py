#!/usr/bin/env python3
"""Scout probe for chem-70c: does device_odmd.centered_frame remove PDS's det(M)->0 accuracy floor?

Not part of the test suite -- exploratory only, to pre-register gate thresholds before writing
tests/test_centered_pds_spec.py (or extending test_moment_pds_spec.py).

Claim under test (specs/BACKLOG.md "Method rungs"):
  1. centered-frame PDS(K) == raw-frame PDS(K) to <1e-6 Ha wherever raw cond(M) < 1e10 (same
     functional, just better conditioned).
  2. centering drops cond(M) by >= 8 orders of magnitude at K=8.
  3. raw PDS is non-monotone past K=4 on N2 CAS(6,6) (0.044 mHa at K=6 -> 0.049 at K=7); centered
     frame stays monotone.
"""
import numpy as np

from device_odmd import centered_frame
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from moment_expansion import hamiltonian_moments, pds_energy

MAX_K = 8
MOMENT_ORDER = 2 * MAX_K  # need up to 2K-1 for PDS(K); pad by 1


def cond_at_K(moments, K):
    M = np.array([[moments[2 * K - i - j] for j in range(1, K + 1)] for i in range(1, K + 1)])
    return np.linalg.cond(M)


def scout(name, mh):
    fci = mh.ground_state_energy()
    mh_c, tau, mu_shift = centered_frame(mh)

    mu_raw, off_raw = hamiltonian_moments(mh, MOMENT_ORDER)
    mu_c, off_c = hamiltonian_moments(mh_c, MOMENT_ORDER)

    print(f"=== {name} (n_qubits={mh.num_qubits}, FCI={fci:.6f}, mu_shift={mu_shift:.4f} Ha, "
          f"tau={tau:.4f}) ===")
    header = (f"{'K':>2} {'cond_raw':>12} {'cond_cent':>12} {'orders_drop':>11} "
              f"{'E_raw':>12} {'E_cent':>12} {'diff(mHa)':>10} "
              f"{'err_raw(mHa)':>12} {'err_cent(mHa)':>13}")
    print(header)
    print("-" * len(header))

    raw_energies, cent_energies, raw_conds = [], [], []
    for K in range(1, MAX_K + 1):
        try:
            e_raw = pds_energy(mu_raw, K, off_raw)
            c_raw = cond_at_K(mu_raw, K)
        except (ValueError, np.linalg.LinAlgError) as e:
            print(f"{K:2d} raw FAILED: {e}")
            e_raw, c_raw = float("nan"), float("nan")
        try:
            e_cent = pds_energy(mu_c, K, off_c)
            c_cent = cond_at_K(mu_c, K)
        except (ValueError, np.linalg.LinAlgError) as e:
            print(f"{K:2d} centered FAILED: {e}")
            e_cent, c_cent = float("nan"), float("nan")

        raw_energies.append(e_raw)
        cent_energies.append(e_cent)
        raw_conds.append(c_raw)

        drop = np.log10(c_raw / c_cent) if c_cent > 0 and c_raw > 0 else float("nan")
        diff_mha = (e_cent - e_raw) * 1e3
        err_raw = (e_raw - fci) * 1e3
        err_cent = (e_cent - fci) * 1e3
        print(f"{K:2d} {c_raw:12.3e} {c_cent:12.3e} {drop:11.2f} "
              f"{e_raw:12.6f} {e_cent:12.6f} {diff_mha:10.4f} "
              f"{err_raw:12.4f} {err_cent:13.4f}")

    print()
    print("Monotonicity (E(K+1) <= E(K) + 1e-9):")
    for K in range(1, MAX_K):
        raw_mono = cent_mono = "n/a"
        if not (np.isnan(raw_energies[K - 1]) or np.isnan(raw_energies[K])):
            raw_mono = raw_energies[K] <= raw_energies[K - 1] + 1e-9
        if not (np.isnan(cent_energies[K - 1]) or np.isnan(cent_energies[K])):
            cent_mono = cent_energies[K] <= cent_energies[K - 1] + 1e-9
        print(f"  K={K}->{K + 1}: raw_monotone={raw_mono}  centered_monotone={cent_mono}")

    print()
    print("Agreement check (centered vs raw, only where raw cond(M) < 1e10):")
    for K in range(1, MAX_K + 1):
        if raw_conds[K - 1] < 1e10:
            diff = abs(cent_energies[K - 1] - raw_energies[K - 1])
            print(f"  K={K}: cond_raw={raw_conds[K - 1]:.3e} diff={diff:.3e} Ha "
                  f"{'OK' if diff < 1e-6 else 'VIOLATION'}")
    print()
    return raw_energies, cent_energies, raw_conds


if __name__ == "__main__":
    print("=== H4 equilibrium-ish (sanity check, matches SPEC_moment_pds geometry) ===")
    mh_h4 = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0")
    scout("H4", mh_h4)

    print("=== N2 CAS(6,6), R=1.10 A (equilibrium-ish, matches benchmark_n2.py) ===")
    mh_n2 = build_molecular_hamiltonian(
        atom="N 0 0 0; N 0 0 1.10", basis="sto3g",
        active_electrons=6, active_orbitals=6,
    )
    scout("N2_CAS(6,6)", mh_n2)
