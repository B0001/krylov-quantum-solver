"""
Acceptance gates G1-G4 for specs/SPEC_classical_shadows.md (classical shadows).

Test-first: ``classical_shadows`` does not exist yet, so this file is RED until the spec is
implemented. Random-Pauli classical shadows estimate <psi|H|psi> from randomized single-qubit
measurements; we check the estimator is unbiased (HF and FCI states), converges as 1/sqrt(shots),
and that its single-shot variance is bounded by the HKP shadow norm sum_k |c_k|^2 3^{w_k} -- whose
3^{weight} growth is the honest sampling-cost finding. Reference: exact <psi|H|psi>.

Seeded (deterministic). PySCF/qiskit only (no block2); `make gates` runs it in its own process.
"""
import numpy as np

from classical_shadows import (
    collect_classical_shadow,
    shadow_energy_samples,
    shadow_norm,
)
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian


def _h2():
    return build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.74")


def _exact(op, psi):
    return float(np.real(psi.conj() @ op.to_matrix() @ psi))


def test_G1_unbiased_on_hf_and_fci_states():
    """The shadow energy mean is within 4 stderr of the exact expectation, HF and FCI states."""
    mh = _h2()
    op = mh.qubit_hamiltonian
    hf = np.asarray(mh.hf_state().data, dtype=complex)
    _, V = np.linalg.eigh(op.to_matrix())
    fci = V[:, 0]
    for psi, seed in ((hf, 1), (fci, 4)):
        bases, signs = collect_classical_shadow(psi, mh.num_qubits, 16000, seed=seed)
        s = shadow_energy_samples(bases, signs, op)
        stderr = s.std() / np.sqrt(s.size)
        assert abs(s.mean() - _exact(op, psi)) < 4.0 * stderr, (s.mean(), _exact(op, psi), stderr)


def test_G2_converges_with_shots():
    """DEFINITION OF DONE: more shots -> smaller |estimate - exact| (1/sqrt(shots))."""
    mh = _h2()
    op = mh.qubit_hamiltonian
    psi = np.asarray(mh.hf_state().data, dtype=complex)
    exact = _exact(op, psi)
    errs = []
    for shots in (500, 16000):
        bases, signs = collect_classical_shadow(psi, mh.num_qubits, shots, seed=7)
        errs.append(abs(shadow_energy_samples(bases, signs, op).mean() - exact))
    assert errs[1] < errs[0], errs


def test_G3_variance_within_shadow_norm():
    """The empirical single-shot variance is bounded by the HKP shadow norm."""
    mh = _h2()
    op = mh.qubit_hamiltonian
    psi = np.asarray(mh.hf_state().data, dtype=complex)
    bases, signs = collect_classical_shadow(psi, mh.num_qubits, 8000, seed=2)
    var = shadow_energy_samples(bases, signs, op).var()
    norm = shadow_norm(op)
    assert 0 < var <= norm, (var, norm)


def test_G4_high_weight_cost():
    """THE FINDING: the shadow norm is inflated by 3^{weight}; high-weight terms carry real cost."""
    mh = _h2()
    op = mh.qubit_hamiltonian
    weights = np.array([int(np.sum(p.x | p.z)) for p in op.paulis])
    c2 = np.abs(op.coeffs) ** 2
    norm = shadow_norm(op)
    capped = float(np.sum(c2 * 3.0 ** np.minimum(weights, 1)))   # weight capped at 1
    assert norm > capped, (norm, capped)                         # 3^w grows past weight 1
    hi_share = float(np.sum((c2 * 3.0 ** weights)[weights >= 3])) / norm
    assert hi_share > 0.05 and weights.max() >= 3, (hi_share, weights.max())
