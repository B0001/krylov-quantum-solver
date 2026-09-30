import numpy as np

from certified_dipole import certified_dipole
from certified_gaps import gap_bracket, reachable_gap
from hf_overlap_certificate import certify_hf_overlap, exact_reachable_overlap
from hf_overlap_subspace import certify_hf_subspace_overlap, exact_hf_subspace_overlap
from hybrid_quantum_solver.molecular_hamiltonian import build_dipole_operators, build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from reachability import reachable_eigenpairs

GEOMS = {
    "persistent-escape H4 (0.9/2.0/2.0)": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 2.9; H 0 0 4.9"),
    "persistent-escape H4 (0.9/2.5/2.0)": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 3.4; H 0 0 5.4"),
}
DIMS = (6, 8, 12, 16)

for name, spec in GEOMS.items():
    mh = build_molecular_hamiltonian(**spec)
    solver = QuantumKrylovSolver(mh)
    gap = reachable_gap(mh)
    print("=" * 100)
    print(f"{name}   exact reachable gap = {gap:.6f} Ha")
    print("  -- upstream gap_bracket (self) --")
    for m in DIMS:
        br = gap_bracket(mh, m, solver=solver)
        print(f"    M={m:2d} Delta_lo={br.gap_lower:+.6f} Delta_hi={br.gap_upper:+.6f} "
              f"lower_escape={gap < br.gap_lower}")

    Az = build_dipole_operators(**spec)[2].to_matrix(sparse=True)
    psi_ex = reachable_eigenpairs(mh)[1][:, 0]
    mu_exact = float((psi_ex.conj() @ (Az @ psi_ex)).real)
    print(f"  -- downstream certified_dipole (exact mu_z = {mu_exact:+.6f}) --")
    for m in DIMS:
        cd = certified_dipole(mh, Az, m, solver=solver)
        if cd.finite:
            inside = cd.mu - cd.half_width <= mu_exact <= cd.mu + cd.half_width
            print(f"    M={m:2d} mu={cd.mu:+.6f} half_width={cd.half_width:.6f} "
                  f"gap_lower_used={cd.gap_lower:+.6f} inside={inside}")
        else:
            print(f"    M={m:2d} VACUOUS gap_lower_used={cd.gap_lower:+.6f}")

    exact1 = exact_reachable_overlap(mh)
    print(f"  -- downstream certify_hf_overlap d=1 (exact = {exact1:.6f}) --")
    for m in DIMS:
        c1 = certify_hf_overlap(mh, m, solver=solver)
        if c1.vacuous:
            print(f"    M={m:2d} VACUOUS  reason={c1.vacuous_reason}")
        else:
            ok = c1.gamma_min <= exact1 + 1e-9
            print(f"    M={m:2d} gamma_min={c1.gamma_min:.6f} valid_lower_bound={ok}")

    exact2 = exact_hf_subspace_overlap(mh, 2)
    print(f"  -- downstream certify_hf_subspace_overlap d=2 (exact = {exact2:.6f}) --")
    for m in DIMS:
        try:
            c2 = certify_hf_subspace_overlap(mh, 2, m=m, solver=solver)
        except ValueError as e:
            print(f"    M={m:2d} RAISED {e}")
            continue
        if c2.vacuous:
            print(f"    M={m:2d} VACUOUS  reason={c2.vacuous_reason}")
        else:
            ok = c2.gamma_min <= exact2 + 1e-9
            print(f"    M={m:2d} gamma_min={c2.gamma_min:.6f} valid_lower_bound={ok}")
