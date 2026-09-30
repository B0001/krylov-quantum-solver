"""One-off probe for chem-ayr: does tightening PySCFDriver's conv_tol (now threaded through
build_molecular_hamiltonian) move any number actually consumed by a shipped gate?

For every geometry used by the certified-arc gates that reference REACHABLE_TOL_CERTIFIED / the
1e-8 family (SPEC_reachability_tolerance's enumerated sites, plus test_chained_overlap_spec.py's
DIRECT_SURVIVES/DIRECT_VACUOUS set), build the Hamiltonian at the driver default (1e-9) and at
1e-13, and diff: ground-state energy, HF energy, and the full HF-population spectrum. Not a pytest
gate -- a pre-registered, throwaway measurement to decide the killable check before touching any
shipped file.
"""
import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

GEOMETRIES = {
    # test_chained_overlap_spec.py DIRECT_SURVIVES + DIRECT_VACUOUS
    "H2_EQ": ("H 0 0 0; H 0 0 0.74", {}),
    "H2_STRETCHED": ("H 0 0 0; H 0 0 2.0", {}),
    "H4_LINEAR_chained": ("H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0", {}),
    "H4_SQUARE_105": ("H 0 0 0; H 1.05 0 0; H 1.05 1.05 0; H 0 1.05 0", {}),
    "H6_LINEAR": ("H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0; H 0 0 4.0; H 0 0 5.0", {}),
    # test_hf_overlap_certificate_spec.py
    "H4_certificate": ("H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0", {}),
    # test_hf_overlap_subspace_spec.py square-H4 sweep (non-artifact sides)
    "H4_SQUARE_140": ("H 0 0 0; H 1.4 0 0; H 1.4 1.4 0; H 0 1.4 0", {}),
    "H4_SQUARE_120": ("H 0 0 0; H 1.2 0 0; H 1.2 1.2 0; H 0 1.2 0", {}),
    "H4_SQUARE_100": ("H 0 0 0; H 1.0 0 0; H 1.0 1.0 0; H 0 1.0 0", {}),
    # test_certified_gaps_spec.py
    "gaps_H4": ("H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7", {}),
    "gaps_LiH": ("Li 0 0 0; H 0 0 1.6", dict(active_electrons=2, active_orbitals=5)),
    "gaps_N2": ("N 0 0 0; N 0 0 1.1", dict(active_electrons=6, active_orbitals=6)),
    # test_certified_dipole_spec.py
    "dipole_HeHp": ("He 0 0 0; H 0 0 0.772", dict(charge=1)),
    "dipole_LiH": ("Li 0 0 0; H 0 0 1.6", {}),
    "dipole_square_H4_140": ("H 0 0 0; H 1.4 0 0; H 1.4 1.4 0; H 0 1.4 0", {}),
    # test_certified_noise_spec.py
    "noise_H4": ("H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7", {}),
    # the two geometries excluded from SPEC_chained_overlap -- these SHOULD move (positive control)
    "EXCLUDED_square_110": ("H 0 0 0; H 1.10 0 0; H 1.10 1.10 0; H 0 1.10 0", {}),
    "EXCLUDED_square_135": ("H 0 0 0; H 1.35 0 0; H 1.35 1.35 0; H 0 1.35 0", {}),
}


def hf_population_spectrum(mh):
    H = mh.qubit_hamiltonian.to_matrix()
    u = np.asarray(mh.hf_state().data, dtype=complex)
    w, V = np.linalg.eigh(H)
    return w, np.abs(V.conj().T @ u) ** 2


def main():
    print(f"{'geometry':24s} {'dE_ground':>12s} {'dE_HF':>12s} {'max|dp0|':>12s} {'moved?':>8s}")
    for name, (atom, kwargs) in GEOMETRIES.items():
        mh_loose = build_molecular_hamiltonian(atom=atom, conv_tol=1e-9, **kwargs)
        mh_tight = build_molecular_hamiltonian(atom=atom, conv_tol=1e-13, **kwargs)

        dE_ground = abs(mh_loose.ground_state_energy() - mh_tight.ground_state_energy())
        dE_hf = abs(mh_loose.hf_energy - mh_tight.hf_energy)

        w_l, p_l = hf_population_spectrum(mh_loose)
        w_t, p_t = hf_population_spectrum(mh_tight)
        # Compare eigenvalues first (they must line up before comparing populations index-by-index).
        dE_spectrum = float(np.max(np.abs(w_l - w_t)))
        dp0 = float(np.max(np.abs(p_l - p_t)))

        moved = dE_ground > 1e-8 or dp0 > 1e-6
        print(f"{name:24s} {dE_ground:12.3e} {dE_hf:12.3e} {dp0:12.3e} "
              f"{'YES' if moved else 'no':>8s}  (max|dE_spectrum|={dE_spectrum:.3e})")


if __name__ == "__main__":
    main()
