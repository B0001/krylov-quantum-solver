"""
Spec gates for SPEC_subspace_floor_resolvability: the self-mode subspace E_d floor is a HEURISTIC,
not rigorous, for d >= 2, and the disjoint-Weinstein-interval guard is a pre-filter that is
demonstrably INSUFFICIENT (a parallel falsification sweep found guard-passing invalid floors).

  G1  the floor bug (regression, killable): raw theta_d - sigma_d > true E_d on H6 R=1.2 d=3.
  G2  pre-filter rejects gross failure: guarded self-mode is VACUOUS on H6 R=1.2 d=3.
  G3  no over-rejection on resolved cases: square-H4 / linear-H4 d=2 pass, non-vacuous & valid.
  G4  THE FALSIFICATION (killable): a guard-PASSING case with an invalid floor exists
      (linear H6 R=1.0 d=3 M=16) -- the guard is not sound.
  G4b the rigorous path (zero-tol): oracle-mode certificates are valid across the sweep.

Each unique molecule is built once (module fixture) with its dense eigendecomposition cached, so
exact ||P_S u|| and the true E_d floor are derived without repeated O(2^n) eighs.
"""

import numpy as np
import pytest

from hf_overlap_subspace import _weinstein_intervals_disjoint, certify_hf_subspace_overlap
from hybrid_quantum_solver.certified_overlap import residual_norm
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver


def _chain(n, R):
    return "; ".join(f"H 0 0 {i * R}" for i in range(n))


def _square(a):
    return f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0"


_ATOMS = {
    "h6_r12": _chain(6, 1.2),   # the floor-bug case (d=3)
    "h6_r10": _chain(6, 1.0),   # the guard-insufficiency witness (d=3, M=16)
    "sq12": _square(1.2), "sq10": _square(1.0), "linh4": _chain(4, 1.0),  # resolved d=2
}


@pytest.fixture(scope="module")
def B():
    """Build each molecule once; cache the dense eigendecomposition and the reachable indices."""
    out = {}
    for key, atom in _ATOMS.items():
        mh = build_molecular_hamiltonian(atom=atom)
        w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
        u = np.asarray(mh.hf_state().data, dtype=complex)
        reach = np.where(np.abs(V.conj().T @ u) ** 2 > 1e-8)[0]
        out[key] = dict(mh=mh, solver=QuantumKrylovSolver(mh), off=mh.energy_offset,
                        Hs=mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc(),
                        w=w, V=V, u=u, reach=reach)
    return out


def _exact_sub(b, d):
    return float(np.linalg.norm(b["V"][:, b["reach"][:d]].conj().T @ b["u"]))


def _e_d_total(b, d):
    return float(b["w"][b["reach"][d]]) + b["off"]


def _ritz(b, d, m):
    en, st = b["solver"].eigenstates(m, n_states=d + 1)
    centers = [en[k] - b["off"] for k in range(d + 1)]
    sig = [residual_norm(b["Hs"], st[k], centers[k]) for k in range(d + 1)]
    return centers, sig


def _raw_floor_total(b, d, m):
    centers, sig = _ritz(b, d, m)
    return centers[d] - sig[d] + b["off"]


@pytest.mark.parametrize("m", [8, 12])
def test_G1_raw_floor_unsound_on_h6_d3(B, m):
    """G1 (the floor bug): the raw self-mode floor exceeds the true reachable E_d on H6 R=1.2 d=3."""
    b = B["h6_r12"]
    assert _raw_floor_total(b, 3, m) > _e_d_total(b, 3) + 1e-6, (
        f"G1 premise gone at M={m}: raw floor became valid -- revisit the spec"
    )


@pytest.mark.parametrize("m", [8, 12])
def test_G2_prefilter_rejects_gross_failure(B, m):
    """G2: the guard rejects the gross non-resolution -> VACUOUS on H6 R=1.2 d=3."""
    b = B["h6_r12"]
    c = certify_hf_subspace_overlap(b["mh"], 3, m=m, solver=b["solver"])
    assert c.vacuous and c.gamma_min == 0.0
    assert "unresolved" in (c.vacuous_reason or "").lower()


def test_G2_oracle_mode_still_valid_on_h6_d3(B):
    """Oracle mode (rigorous) returns a valid bound where self-mode is rejected."""
    b = B["h6_r12"]
    c = certify_hf_subspace_overlap(b["mh"], 3, m=8, e_d=_e_d_total(b, 3), solver=b["solver"])
    floor = 0.0 if c.vacuous else c.gamma_min
    assert floor <= _exact_sub(b, 3) + 1e-12


@pytest.mark.parametrize("key", ["sq12", "sq10", "linh4"])
def test_G3_no_over_rejection_on_resolved_d2(B, key):
    """G3: resolved d=2 cases pass the guard, non-vacuous and valid (regression-guards PR #20)."""
    b = B[key]
    c = certify_hf_subspace_overlap(b["mh"], 2, m=8, solver=b["solver"])
    assert not c.vacuous, f"guard over-rejected a resolved case: {key}"
    assert c.gamma_min <= _exact_sub(b, 2) + 1e-12


def test_G4_guard_is_insufficient_witness(B):
    """G4 (the falsification, killable): a guard-PASSING case with an INVALID floor exists.

    Deterministic witness: linear H6, R=1.0 A, d=3, M=16. The Weinstein intervals are disjoint
    (guard passes) yet beta_self > true E_d (invalid floor) -- a reachable level of tiny HF
    amplitude near the cluster boundary is localized as a higher level, and more Krylov dimension
    does not fix it. So the guard is a heuristic pre-filter, NOT a soundness proof; oracle mode is
    the only rigorous path. If a future change makes this case sound, this test flips -> revisit.
    """
    b = B["h6_r10"]
    centers, sig = _ritz(b, 3, 16)
    assert _weinstein_intervals_disjoint(centers, sig), "witness premise: guard should PASS here"
    beta_self = centers[3] - sig[3] + b["off"]
    assert beta_self > _e_d_total(b, 3) + 1e-6, (
        "witness premise: floor should be INVALID here; if it became valid the guard may now be "
        "sound on this case -- revisit SPEC_subspace_floor_resolvability"
    )


@pytest.mark.parametrize("key,d", [("sq12", 2), ("linh4", 2), ("h6_r12", 3), ("h6_r10", 3)])
def test_G4b_oracle_mode_is_sound(B, key, d):
    """G4b (the rigorous path, zero-tol): oracle-mode certificates are valid across the sweep."""
    b = B[key]
    c = certify_hf_subspace_overlap(b["mh"], d, m=8, e_d=_e_d_total(b, d), solver=b["solver"])
    floor = 0.0 if c.vacuous else c.gamma_min
    assert floor <= _exact_sub(b, d) + 1e-12


def test_guard_unit():
    """_weinstein_intervals_disjoint: disjoint -> True, overlapping -> False."""
    assert _weinstein_intervals_disjoint([0.0, 1.0, 2.0], [0.1, 0.1, 0.1])
    assert not _weinstein_intervals_disjoint([0.0, 1.0, 2.0], [0.6, 0.6, 0.1])
    assert _weinstein_intervals_disjoint([0.0], [0.5])


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
