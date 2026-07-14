"""
Acceptance gates G1-G4 for specs/SPEC_qubitization_spectrum.md (the qubitization walk operator
recovers EVERY Hamiltonian eigenvalue, exactly twice -- not just "every recovered value is valid",
the one-directional check `qubitization_blueprint.verify_qubitization` already does).

Deliberately no new library code: the bidirectional/counting check lives entirely here, reusing
`qubitization_blueprint.py`'s existing `build_qubit_hamiltonian`/`pauli_decompose`/
`build_walk_operator` unmodified, so it is a genuine external verification rather than a round-trip
through `verify_qubitization`'s own one-directional logic.
"""
import numpy as np
import pytest
from pyscf import ao2mo, gto, mcscf, scf

from qubitization_blueprint import build_qubit_hamiltonian, build_walk_operator, pauli_decompose

TOL = 1e-6

SYSTEMS = ("H2", "LiH")


def _reference(label):
    if label == "H2":
        atom, norb, ne = "H 0 0 0; H 0 0 0.74", 2, 2
    else:
        atom, norb, ne = "Li 0 0 0; H 0 0 1.6", 2, 2
    mol = gto.M(atom=atom, basis="sto-3g")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    cas = mcscf.CASCI(mf, norb, ne)
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), norb)
    return h1, eri, norb


def _build(label):
    h1, eri, norb = _reference(label)
    H, n = build_qubit_hamiltonian(h1, eri, norb)
    eigH = np.linalg.eigvalsh(H).real          # NOT deduplicated -- one entry per degenerate slot
    terms = pauli_decompose(H, n)
    W, lam, L, a = build_walk_operator(terms, n)
    theta = np.angle(np.linalg.eigvals(W))
    recovered_all = lam * np.cos(theta)         # NOT deduplicated -- one entry per W eigenvalue
    return eigH, recovered_all, lam, terms, W, n, a


@pytest.mark.parametrize("label", SYSTEMS)
def test_G1_every_true_eigenvalue_is_recovered(label):
    """The direction verify_qubitization does NOT test: every eigenvalue of H has a matching
    W-recovered value nearby. A dropped eigenvalue (e.g. the ground state) would fail this."""
    eigH, recovered_all, *_ = _build(label)
    max_err = max(np.min(np.abs(recovered_all - e)) for e in eigH)
    assert max_err < TOL, (label, max_err)


@pytest.mark.parametrize("label", SYSTEMS)
def test_G2_exactly_two_recovered_phases_per_eigenvalue_slot(label):
    """THE FINDING / definition of done: counted from the RECOVERED (W) side -- how many of W's
    (non-deduplicated) eigenphases land within tolerance of SOME true eigenvalue -- the total is
    exactly 2 * dim(H): every eigenvalue slot contributes exactly one theta/-theta pair, no
    eigenvalue missing, no extra collision. (Summing local match-counts from the eigH side instead
    double-counts whenever two eigenvalues are within TOL of each other -- caught while writing
    this gate: an early version did exactly that and over-counted by 60 vs the correct 32.)"""
    eigH, recovered_all, *_ = _build(label)
    matched = sum(1 for r in recovered_all if np.min(np.abs(r - eigH)) < TOL)
    assert matched == 2 * len(eigH), (label, matched, len(eigH))


@pytest.mark.parametrize("label", SYSTEMS)
def test_G3_theta_minus_theta_pairing_never_odd(label):
    """Every individual eigenvalue's multiplicity count is even -- the theta/-theta structural
    pairing holds without exception, never an unpaired phase."""
    eigH, recovered_all, *_ = _build(label)
    counts = [int(np.sum(np.abs(recovered_all - e) < TOL)) for e in np.unique(np.round(eigH, 8))]
    assert all(c % 2 == 0 for c in counts), (label, counts)
    assert all(c >= 2 for c in counts), (label, counts)


@pytest.mark.parametrize("label", SYSTEMS)
def test_G4_construction_bookkeeping(label):
    """lambda equals the independently-recomputed sum of |Pauli coefficients|, and W's dimension
    is exactly ancilla_dim * system_dim -- pins the basic bookkeeping, not just the spectral
    relation."""
    _eigH, _recovered_all, lam, terms, W, n, a = _build(label)
    lam_recomputed = sum(abs(c) for _, c in terms)
    assert abs(lam - lam_recomputed) < 1e-10, (label, lam, lam_recomputed)
    assert W.shape[0] == (2 ** a) * (2 ** n), (label, W.shape[0], a, n)
