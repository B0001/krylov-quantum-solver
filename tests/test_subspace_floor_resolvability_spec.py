"""
Spec gates for SPEC_subspace_floor_resolvability: the self-mode subspace E_d floor is not
rigorous for d >= 3, and the disjoint-Weinstein-interval guard makes it fail-safe.

THE FINDING: certify_hf_subspace_overlap's self-mode floor theta_d - sigma_d can EXCEED the true
reachable E_d (an unsound floor) when the Krylov space has not resolved the cluster -- demonstrated
on linear H6 (R=1.2, d=3). THE FIX: the guard returns VACUOUS when the d+1 lowest Weinstein
intervals overlap, so the library never emits the unsound positive.

  G1 the bug (regression, killable): raw self-mode floor theta_d - sigma_d > true E_d on H6 d=3.
  G2 guard catches it: guarded certificate is VACUOUS on H6 d=3 (M in {8,12}).
  G3 guard does not over-reject: square-H4 / linear-H4 d=2 pass, non-vacuous and valid.
  G4 empirical soundness (zero-tol): guard-passes => valid floor AND valid certificate, no escape.
"""

import pytest

from hf_overlap_subspace import (
    _reachable_e_d_total,
    _weinstein_intervals_disjoint,
    certify_hf_subspace_overlap,
    exact_hf_subspace_overlap,
)
from hybrid_quantum_solver.certified_overlap import residual_norm
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver


def _chain(n, R):
    return "; ".join(f"H 0 0 {i * R}" for i in range(n))


def _square(a):
    return f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0"


def _raw_self_floor_total(mh, d, m, solver):
    """The UNGUARDED self-mode floor theta_d - sigma_d as a total energy (what the bug produces)."""
    off = mh.energy_offset
    Hs = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    en, st = solver.eigenstates(m, n_states=d + 1)
    th_e = en[d] - off
    return th_e - residual_norm(Hs, st[d], th_e) + off


# ---- H6 d=3 is expensive (12 qubits); build once ----
@pytest.fixture(scope="module")
def h6_r12():
    mh = build_molecular_hamiltonian(atom=_chain(6, 1.2))
    return {"mh": mh, "solver": QuantumKrylovSolver(mh),
            "true_Ed": _reachable_e_d_total(mh, 3), "exact": exact_hf_subspace_overlap(mh, 3)}


@pytest.mark.parametrize("m", [8, 12])
def test_G1_raw_self_floor_is_unsound_on_h6_d3(h6_r12, m):
    """G1 (the bug): the raw self-mode floor exceeds the true reachable E_d on H6 d=3."""
    raw = _raw_self_floor_total(h6_r12["mh"], 3, m, h6_r12["solver"])
    assert raw > h6_r12["true_Ed"] + 1e-6, (
        f"G1 premise gone: raw self floor {raw} <= true E_d {h6_r12['true_Ed']} at M={m} "
        "-- if the raw floor became rigorous here, revisit the spec"
    )


@pytest.mark.parametrize("m", [8, 12])
def test_G2_guard_returns_vacuous_on_h6_d3(h6_r12, m):
    """G2: the guarded certificate is VACUOUS on H6 d=3 -- never the unsound positive floor."""
    c = certify_hf_subspace_overlap(h6_r12["mh"], 3, m=m, solver=h6_r12["solver"])
    assert c.vacuous, f"G2: expected VACUOUS on unresolved H6 d=3 at M={m}, got gamma={c.gamma_min}"
    assert c.gamma_min == 0.0
    assert "unresolved" in (c.vacuous_reason or "").lower()


def test_G2_oracle_mode_still_works_on_h6_d3(h6_r12):
    """The guard is self-mode only: oracle mode (rigorous) still returns a valid bound."""
    c = certify_hf_subspace_overlap(h6_r12["mh"], 3, m=8, e_d=h6_r12["true_Ed"],
                                    solver=h6_r12["solver"])
    floor = 0.0 if c.vacuous else c.gamma_min
    assert floor <= h6_r12["exact"] + 1e-12


@pytest.mark.parametrize("atom,d", [(_square(1.4), 2), (_square(1.2), 2), (_square(1.0), 2),
                                    (_chain(4, 1.0), 2)])
def test_G3_guard_does_not_over_reject(atom, d):
    """G3: on the known-good (resolved) cases the guard passes -- non-vacuous and valid."""
    mh = build_molecular_hamiltonian(atom=atom)
    solver = QuantumKrylovSolver(mh)
    c = certify_hf_subspace_overlap(mh, d, m=8, solver=solver)
    assert not c.vacuous, f"G3: guard over-rejected a resolved case {atom} d={d}"
    assert c.gamma_min <= exact_hf_subspace_overlap(mh, d) + 1e-12


# G4: empirical soundness sweep -- guard-passes => valid floor and valid certificate, zero escapes.
_SWEEP = [
    (_square(1.4), 2), (_square(1.2), 2), (_square(1.0), 2),
    (_chain(4, 1.0), 2), (_chain(4, 1.5), 2),
    (_chain(6, 1.2), 2), (_chain(6, 1.2), 3), (_chain(6, 1.8), 2), (_chain(6, 1.8), 3),
]


@pytest.mark.parametrize("atom,d", _SWEEP)
@pytest.mark.parametrize("m", [6, 8, 12])
def test_G4_guard_soundness_no_escape(atom, d, m):
    """G4 (zero-tol): every case the guard PASSES has a valid self-mode floor and certificate."""
    mh = build_molecular_hamiltonian(atom=atom)
    solver = QuantumKrylovSolver(mh)
    off = mh.energy_offset
    Hs = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    en, st = solver.eigenstates(m, n_states=d + 1)
    centers = [en[k] - off for k in range(d + 1)]
    sigmas = [residual_norm(Hs, st[k], centers[k]) for k in range(d + 1)]
    if not _weinstein_intervals_disjoint(centers, sigmas):
        pytest.skip("guard rejects this case (unresolved) -- soundness claim is over passes only")
    # Guard passed: the self-mode floor must be a valid lower bound on the true reachable E_d,
    beta_self = centers[d] - sigmas[d] + off
    assert beta_self <= _reachable_e_d_total(mh, d) + 1e-9, (
        f"G4 ESCAPE {atom} d={d} M={m}: guard passed but floor {beta_self} > true E_d"
    )
    # and the resulting certificate must be valid.
    c = certify_hf_subspace_overlap(mh, d, m=m, solver=solver)
    floor = 0.0 if c.vacuous else c.gamma_min
    assert floor <= exact_hf_subspace_overlap(mh, d) + 1e-12


def test_guard_unit():
    """_weinstein_intervals_disjoint: disjoint -> True, overlapping -> False."""
    assert _weinstein_intervals_disjoint([0.0, 1.0, 2.0], [0.1, 0.1, 0.1])
    assert not _weinstein_intervals_disjoint([0.0, 1.0, 2.0], [0.6, 0.6, 0.1])  # I0,I1 overlap
    assert _weinstein_intervals_disjoint([0.0], [0.5])  # single interval trivially disjoint


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
