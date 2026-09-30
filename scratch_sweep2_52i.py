import itertools
import numpy as np

from certified_dipole import certified_dipole
from certified_gaps import gap_bracket, reachable_gap
from hf_overlap_certificate import certify_hf_overlap, exact_reachable_overlap
from hf_overlap_subspace import certify_hf_subspace_overlap, exact_hf_subspace_overlap
from hybrid_quantum_solver.molecular_hamiltonian import build_dipole_operators, build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from reachability import reachable_eigenpairs

DIMS = (6, 8, 12, 16)


def check_reproducible(atom, n=2):
    vals = []
    for _ in range(n):
        mh = build_molecular_hamiltonian(atom=atom)
        vals.append(reachable_gap(mh))
    return max(vals) - min(vals) < 1e-6 * max(1.0, abs(vals[0])), vals[0]


# systematic asymmetric spacings, kept under 8.5 Ang max separation to dodge the SCF
# near-degeneracy instability found for R>~9 Ang.
grid = []
for d1 in (0.9, 1.0, 1.1):
    for d2 in (1.5, 2.0, 2.5):
        for d3 in (2.0, 3.0, 4.0):
            x1, x2, x3 = d1, d1 + d2, d1 + d2 + d3
            if x3 > 8.3:
                continue
            grid.append((d1, d2, d3, f"H 0 0 0; H 0 0 {x1}; H 0 0 {x2}; H 0 0 {x3}"))

print(f"{len(grid)} candidate geometries")
escapers = []
for d1, d2, d3, atom in grid:
    stable, gap = check_reproducible(atom)
    if not stable:
        continue
    mh = build_molecular_hamiltonian(atom=atom)
    solver = QuantumKrylovSolver(mh)
    for m in DIMS:
        br = gap_bracket(mh, m, solver=solver)
        if gap < br.gap_lower:  # lower-certificate escape
            escapers.append((d1, d2, d3, atom, m, gap, br.gap_lower, br.gap_upper))

print(f"{len(escapers)} (geometry, M) lower-escapes found among stable geometries")
for row in escapers[:40]:
    print(row)
