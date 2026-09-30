import numpy as np
import traceback

from certified_dipole import certified_dipole
from certified_gaps import gap_bracket
from hf_overlap_certificate import certify_hf_overlap
from hf_overlap_subspace import certify_hf_subspace_overlap
from hybrid_quantum_solver.molecular_hamiltonian import build_dipole_operators, build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

for R in (0.5, 0.6, 0.772, 0.9, 1.1, 1.4, 1.8, 2.5):
    spec = dict(atom=f"He 0 0 0; H 0 0 {R}", charge=1)
    try:
        mh = build_molecular_hamiltonian(**spec)
    except Exception as e:
        print(f"R={R} build RAISED {type(e).__name__}: {e}")
        continue
    solver = QuantumKrylovSolver(mh)
    Az = build_dipole_operators(**spec)[2].to_matrix(sparse=True)
    for m in (6, 8, 12, 16, 20):
        try:
            cd = certified_dipole(mh, Az, m, solver=solver)
            dipole_note = f"hw={cd.half_width:.3e}"
        except Exception as e:
            dipole_note = f"RAISED {type(e).__name__}: {e}"
        try:
            c1 = certify_hf_overlap(mh, m, solver=solver)
            hf1_note = "VACUOUS" if c1.vacuous else f"g={c1.gamma_min:.6f}"
        except Exception as e:
            hf1_note = f"RAISED {type(e).__name__}: {e}"
        try:
            c2 = certify_hf_subspace_overlap(mh, 2, m=m, solver=solver)
            hf2_note = "VACUOUS" if c2.vacuous else f"g={c2.gamma_min:.6f}"
        except Exception as e:
            hf2_note = f"RAISED {type(e).__name__}: {e}"
        print(f"R={R:.3f} M={m:2d}  dipole[{dipole_note}]  hf1[{hf1_note}]  hf2[{hf2_note}]")
