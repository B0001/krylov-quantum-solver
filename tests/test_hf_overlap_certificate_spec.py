"""
Spec gates for SPEC_hf_overlap_certificate: certified HF guiding-state overlap on molecules.

Feeds the repo's premise-gated Krylov E1 floor (certified_gaps.gap_bracket, self mode) into
the SPEC-21 Davis-Kahan machinery and certifies gamma_min <= |<HF|psi_0_reachable>| from Krylov
data alone -- the quantity the guided-LH literature assumes but never certifies.

  G1 validity (killable): gamma_min <= exact reachable overlap, zero tolerance.
  G2 usefulness: non-vacuous with gamma_min >= 1/n_qubits at equilibrium.
  G3 ordering: self-mode <= oracle-mode <= exact.
  G4 premise boundary: self mode at M < 6 raises.
  G5 Krylov refinement: self-mode floor non-decreasing in M, -> oracle (H4); self == oracle (H2).
"""

import numpy as np
import pytest

from hf_overlap_certificate import certify_hf_overlap, exact_reachable_overlap
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

# Built once per geometry; the Krylov basis is reused across M via a shared solver.
_GEOMETRIES = {
    "H2_eq": "H 0 0 0; H 0 0 0.74",
    "H2_stretched": "H 0 0 0; H 0 0 2.0",
    "H4": "H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0",
}


def _oracle_e1_total(mh):
    """Exact reachable E1 as a TOTAL energy, for oracle-mode comparison."""
    Hd = mh.qubit_hamiltonian.to_matrix()
    w, V = np.linalg.eigh(Hd)
    hf = np.asarray(mh.hf_state().data, dtype=complex)
    reach = np.where(np.abs(V.conj().T @ hf) ** 2 > 1e-10)[0]
    return float(w[reach[1]]) + mh.energy_offset


@pytest.fixture(scope="module")
def bundles():
    out = {}
    for key, atom in _GEOMETRIES.items():
        mh = build_molecular_hamiltonian(atom=atom)
        out[key] = {
            "mh": mh,
            "solver": QuantumKrylovSolver(mh),
            "exact": exact_reachable_overlap(mh),
            "e1": _oracle_e1_total(mh),
        }
    return out


@pytest.mark.parametrize("key", list(_GEOMETRIES))
@pytest.mark.parametrize("m", [6, 8, 12])
def test_G1_validity_self_and_oracle(bundles, key, m):
    """G1: certified floor never exceeds the exact reachable overlap. Zero tolerance."""
    b = bundles[key]
    c_self = certify_hf_overlap(b["mh"], m, solver=b["solver"])
    c_orac = certify_hf_overlap(b["mh"], m, e1=b["e1"], solver=b["solver"])
    for c, mode in ((c_self, "self"), (c_orac, "oracle")):
        floor = 0.0 if c.vacuous else c.gamma_min
        assert floor <= b["exact"] + 1e-12, (
            f"G1 VIOLATION {key} M={m} {mode}: gamma_min={floor} > exact={b['exact']}"
        )


@pytest.mark.parametrize("key", ["H2_eq", "H4"])
def test_G2_useful_at_equilibrium(bundles, key):
    """G2: at equilibrium the self-mode floor is non-vacuous and clears the 1/n_qubits threshold."""
    b = bundles[key]
    c = certify_hf_overlap(b["mh"], m=6, solver=b["solver"])
    assert not c.vacuous, f"G2: {key} self-mode certificate is vacuous at equilibrium"
    threshold = 1.0 / b["mh"].qubit_hamiltonian.num_qubits
    assert c.gamma_min >= threshold, (
        f"G2: {key} gamma_min={c.gamma_min} below 1/n={threshold}"
    )


@pytest.mark.parametrize("key", list(_GEOMETRIES))
@pytest.mark.parametrize("m", [6, 8, 12])
def test_G3_ordering_self_le_oracle_le_exact(bundles, key, m):
    """G3: self-mode floor <= oracle-mode floor <= exact overlap."""
    b = bundles[key]
    c_self = certify_hf_overlap(b["mh"], m, solver=b["solver"])
    c_orac = certify_hf_overlap(b["mh"], m, e1=b["e1"], solver=b["solver"])
    g_self = 0.0 if c_self.vacuous else c_self.gamma_min
    g_orac = 0.0 if c_orac.vacuous else c_orac.gamma_min
    assert g_self <= g_orac + 1e-12, f"G3: {key} M={m} self {g_self} > oracle {g_orac}"
    assert g_orac <= b["exact"] + 1e-12, f"G3: {key} M={m} oracle {g_orac} > exact {b['exact']}"


def test_G4_self_mode_below_M6_raises(bundles):
    """G4: the self-mode premise boundary (M >= 6) is inherited as a loud raise."""
    b = bundles["H2_eq"]
    with pytest.raises(ValueError, match="m >= 6"):
        certify_hf_overlap(b["mh"], m=4, solver=b["solver"])
    # oracle mode has no such premise -- must NOT raise below M=6
    c = certify_hf_overlap(b["mh"], m=4, e1=b["e1"], solver=b["solver"])
    assert c.gamma_min <= b["exact"] + 1e-12


def test_G5_krylov_refinement_h4_monotone(bundles):
    """G5: on H4 the self-mode floor is non-decreasing in M and converges up to the oracle floor."""
    b = bundles["H4"]
    floors = []
    for m in (6, 8, 12):
        c = certify_hf_overlap(b["mh"], m, solver=b["solver"])
        assert not c.vacuous, f"G5: H4 self-mode vacuous at M={m}"
        floors.append(c.gamma_min)
    assert floors[0] <= floors[1] + 1e-9 <= floors[2] + 1e-9, (
        f"G5: H4 self-mode floors not non-decreasing in M: {floors}"
    )
    oracle = certify_hf_overlap(b["mh"], 12, e1=b["e1"], solver=b["solver"]).gamma_min
    assert floors[-1] <= oracle + 1e-12, f"G5: H4 self floor {floors[-1]} exceeds oracle {oracle}"


def test_G5_h2_self_equals_oracle(bundles):
    """G5: H2's reachable sector is 2-dimensional (sigma_1 = 0), so self mode == oracle mode."""
    b = bundles["H2_eq"]
    c_self = certify_hf_overlap(b["mh"], 8, solver=b["solver"])
    c_orac = certify_hf_overlap(b["mh"], 8, e1=b["e1"], solver=b["solver"])
    assert c_self.gamma_min == pytest.approx(c_orac.gamma_min, abs=1e-9)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
