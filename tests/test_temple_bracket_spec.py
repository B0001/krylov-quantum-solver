"""
Acceptance gates G1-G4 for specs/SPEC_temple_bracket.md (certified two-sided energy brackets).

Test-first: ``temple_bounds`` does not exist yet, so this file is RED until the spec is
implemented. Claim: the QKSD ground Ritz eigenstate + ONE extra expectation <Psi0|H^2|Psi0> gives
a rigorous Temple lower bound, so every Krylov solve carries a certified bracket
[E_Temple, E_Ritz] containing the exact reachable-sector ground energy, closing as M grows.
Weinstein is the looser fallback; the oracle-free mode feeds Temple eps = theta1 - sigma1 from
the same Krylov data and is valid only once the subspace resolves the excited state (M >= 6 here
-- the recorded boundary).

Noiseless by design: the claim is about RIGOR (one containment escape kills the spec), not
sampling. References: exact reachable spectrum by dense diagonalization, as
SPEC_qksd_excited.md. PySCF/qiskit, no block2; `make gates` runs it in its own process.
"""
import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from temple_bounds import bracket_ladder, krylov_bracket

SYSTEMS = {
    "h2": dict(atom="H 0 0 0; H 0 0 0.74"),
    "h4": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
    "lih": dict(atom="Li 0 0 0; H 0 0 1.6"),
    "n2": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
}
DIMS = (2, 4, 6, 8, 12, 16, 20, 24)
_CACHE = {}


def _case(key):
    """(mh, shared solver, exact reachable E0/E1 as total energies)."""
    if key not in _CACHE:
        mh = build_molecular_hamiltonian(**SYSTEMS[key])
        w_eig, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
        psi0 = np.asarray(mh.hf_state().data, dtype=complex)
        pops = np.abs(V.conj().T @ psi0) ** 2
        reach = w_eig[pops > 1e-8].real + mh.energy_offset
        _CACHE[key] = (mh, QuantumKrylovSolver(mh), float(reach[0]), float(reach[1]))
    return _CACHE[key]


def test_G1_containment_no_exceptions():
    """DEFINITION OF DONE: oracle-mode bracket contains the exact energy at EVERY system x M."""
    for key in SYSTEMS:
        mh, solver, e0, e1 = _case(key)
        for br in bracket_ladder(mh, DIMS, eps=e1, solver=solver):
            assert br.lower <= e0 + 1e-9, (key, br.m, br.lower - e0)
            assert br.upper >= e0 - 1e-9, (key, br.m, br.upper - e0)


def test_G2_bracket_closes():
    """width(M=16) < width(M=4) everywhere; micro-Hartree widths at convergence."""
    for key, (m_tight, tol) in {"h4": (16, 1e-5), "n2": (16, 1e-5), "lih": (24, 1e-4)}.items():
        mh, solver, _, e1 = _case(key)
        w4 = krylov_bracket(mh, 4, eps=e1, solver=solver).width
        w_tight = krylov_bracket(mh, m_tight, eps=e1, solver=solver).width
        assert w_tight < w4, (key, w_tight, w4)
        assert w_tight < tol, (key, w_tight)


def test_G3_certification_overhead_small():
    """At mid-convergence the certified width is < 5x the uncertified Ritz error (measured
    ~2.7x): rigor costs a small constant factor, not orders of magnitude."""
    for key, m in (("h4", 8), ("n2", 12)):
        mh, solver, e0, e1 = _case(key)
        br = krylov_bracket(mh, m, eps=e1, solver=solver)
        ritz_err = br.upper - e0
        assert br.width < 5.0 * ritz_err, (key, br.width, ritz_err)


def test_G4_self_consistent_validity_region():
    """(a) eps = theta1 - sigma1 (no oracle) stays a valid lower bound at every M >= 6 and is
    tight (< 1e-4 Ha) by M=24 on all systems; (b) at M=4 its premise eps <= E1 FAILS on H4 and
    N2 -- small-M self-certification is not rigorous (the recorded boundary)."""
    for key in SYSTEMS:
        mh, solver, e0, _ = _case(key)
        for m in (6, 8, 12, 16, 20, 24):
            br = krylov_bracket(mh, m, solver=solver)      # self mode
            assert br.eps_source == "self"
            assert br.lower <= e0 + 1e-9, (key, m, br.lower - e0)
        assert krylov_bracket(mh, 24, solver=solver).width < 1e-4, key
    for key in ("h4", "n2"):
        mh, solver, _, e1 = _case(key)
        br = krylov_bracket(mh, 4, solver=solver)
        assert br.eps > e1, (key, br.eps, e1)              # premise violated at small M
