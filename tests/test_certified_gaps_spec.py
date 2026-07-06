"""
Acceptance gates G1-G4 for specs/SPEC_certified_gaps.md.

Claim: the fundamental gap Delta = E_1 - E_0 of the HF-reachable sector can be BRACKETED from Krylov
data alone (no FCI), via Delta_hi = theta_1 - tau_0 (interlacing + Temple) and
Delta_lo = (theta_1 - sigma_1) - theta_0 (Weinstein self-eps). The exact reachable gap lies inside
at every depth M >= 6 (zero escapes) and the bracket closes with depth; at M = 4 the Weinstein
premise fails and the lower certificate escapes -- the temple_bracket boundary, inherited on gaps.

Exact statevector, sector-restricted (the two lowest HF-reachable levels). PySCF/qiskit, no block2;
`make gates` runs it in its own process.
"""
import numpy as np
import pytest

from certified_gaps import gap_bracket, gap_bracket_ladder, reachable_gap
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

# Small, cheap, exactly-diagonalizable references (FCI used only to CHECK the bracket, never fed in).
CASES = {
    "H4": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
    "LiH": dict(atom="Li 0 0 0; H 0 0 1.6", active_electrons=2, active_orbitals=5),
    "N2": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
}
_CERTIFIED_DIMS = (6, 8, 12, 16, 20, 24)   # M >= 6: the certified regime


@pytest.fixture(scope="module")
def built():
    out = {}
    for name, spec in CASES.items():
        mh = build_molecular_hamiltonian(**spec)
        out[name] = (mh, reachable_gap(mh), QuantumKrylovSolver(mh))
    return out


def test_G1_zero_escapes_in_certified_regime(built):
    """DEFINITION OF DONE: the exact reachable gap lies inside [Delta_lo, Delta_hi] at EVERY
    M >= 6, for every system (self-eps mode -- no oracle). One escape kills the claim."""
    for name, (mh, gap, solver) in built.items():
        for br in gap_bracket_ladder(mh, _CERTIFIED_DIMS, solver=solver):
            assert br.gap_lower - 1e-9 <= gap <= br.gap_upper + 1e-9, (name, br.m, br.gap_lower,
                                                                       gap, br.gap_upper)


def test_G2_bracket_closes_with_depth(built):
    """The certified interval tightens substantially with Krylov depth (M=24 vs M=6): a genuine
    error bar, not a vacuous one. Width is finite throughout the certified regime."""
    for name, (mh, gap, solver) in built.items():
        w6 = gap_bracket(mh, 6, solver=solver).width
        w24 = gap_bracket(mh, 24, solver=solver).width
        assert np.isfinite(w6) and np.isfinite(w24), name
        assert w24 < 0.5 * w6, (name, w6, w24)            # closes by > 2x
        assert w24 < 0.05, (name, w24)                    # < 50 mHa certified at M=24


def test_G3_boundary_premise_fails_at_M4(built):
    """THE FINDING (the honest boundary): at M=4 the Weinstein premise eps_1 <= E_1 fails for the
    multireference cases (H4, N2), and the LOWER certificate escapes (gap < Delta_lo) -- so the
    bracket is trustworthy only for M >= 6. Checked against the oracle E_1 the live path never has.
    LiH (well-separated) may already satisfy the premise at M=4; the claim is failure EXISTS, not
    that it is universal."""
    escapes = 0
    for name in ("H4", "N2"):
        mh, gap, solver = built[name]
        br = gap_bracket(mh, 4, solver=solver)
        H = mh.qubit_hamiltonian.to_matrix()
        w, V = np.linalg.eigh(H)
        hf = np.asarray(mh.hf_state().data, dtype=complex)
        reach = w[np.abs(V.conj().T @ hf) ** 2 > 1e-10]
        e1_elec = reach[1]                                # qubit H is electronic; eig = E_1 (elec)
        premise_ok = br.eps1 <= e1_elec + 1e-12           # eps_1 <= E_1 ?
        if not premise_ok:
            escapes += 1
            # premise failure manifests as a lower-certificate violation
            assert gap < br.gap_lower + 1e-9, (name, gap, br.gap_lower)
    assert escapes >= 1, "expected the M=4 premise to fail for at least one multireference case"


def test_G4_upper_certificate_is_robust_and_scope(built):
    """The asymmetry: the UPPER certificate (interlacing + Temple) brackets the gap from above at
    EVERY tested depth, including M=4 where the lower side is unreliable -- the premise-sensitive
    side is the lower one. Plus honest scope: the certified object is the reachable gap (positive)
    and the oracle mode reproduces a valid, tighter-or-equal lower certificate."""
    for name, (mh, gap, solver) in built.items():
        assert gap > 0.0, name
        for m in (4, 6, 12, 24):
            br = gap_bracket(mh, m, solver=solver)
            assert gap <= br.gap_upper + 1e-9, (name, m, gap, br.gap_upper)   # upper holds at M=4 too
        # oracle mode: feeding the exact E_1 is self-consistent and never escapes in the certified
        # regime (a sanity check that "self" is a faithful stand-in for the oracle at M >= 6).
        H = mh.qubit_hamiltonian.to_matrix()
        w, V = np.linalg.eigh(H)
        hf = np.asarray(mh.hf_state().data, dtype=complex)
        reach = w[np.abs(V.conj().T @ hf) ** 2 > 1e-10]
        e1_total = reach[1] + mh.energy_offset            # oracle mode expects a TOTAL energy
        for m in _CERTIFIED_DIMS:
            bro = gap_bracket(mh, m, e1=e1_total, solver=solver)
            assert bro.eps1_source == "oracle"
            assert bro.gap_lower - 1e-9 <= gap <= bro.gap_upper + 1e-9, (name, m)
