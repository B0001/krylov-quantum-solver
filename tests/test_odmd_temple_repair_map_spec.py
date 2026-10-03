"""
Acceptance gates G-a, G-b, G-c for specs/SPEC_odmd_temple_repair_map.md (chem-z3h): where ODMD's
E_1 estimate repairs -- and where it silently breaks -- the Temple bracket's oracle-free premise.

Claim under test: substituting ODMD's E_1 (first row of the solver's own S) for self mode's
theta_1 - sigma_1 repairs the Temple premise at shallow M with no extra measurement. These gates
pin the measured K-vs-validity map instead: where self mode actually fails (G-a), that ODMD's
overshoot is located and recorded INVALID by containment -- with the oracle-free
`premise_refuted` flag's sound but partial coverage (G-b), and that the free K <= M region does
not repair (G-c).

Noiseless; reference = exact reachable E_0/E_1 by dense diagonalisation (an oracle, validation
only). PySCF/qiskit, no block2; two 12-qubit dense references dominate the run time.
"""
from functools import lru_cache

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from odmd_temple_repair_map import SYSTEMS, TOL, k_vs_validity_map, summarize
from temple_bounds import krylov_bracket


@lru_cache(maxsize=None)
def _map(key):
    mp = k_vs_validity_map(key)
    return mp, summarize(mp)


def test_Ga_self_mode_violation_at_m4_is_not_universal():
    """KILLS 'violated at M=4 on all three systems': violated on H4, N2 and full-space LiH; not on
    LiH CAS(2,5); vacuous on H2 (its reachable sector is exhausted by M=4)."""
    for key in ("h4", "n2", "lih"):
        assert _map(key)[1]["self_m4_mha"] > TOL * 1e3, key
    assert _map("lih_cas")[1]["self_m4_mha"] < -TOL * 1e3
    mp, s = _map("h2")
    assert mp["n_reachable"] <= 4 and s["self_m4_mha"] <= TOL * 1e3, (mp["n_reachable"], s)


def test_Gb_overshoot_located_and_recorded_invalid():
    """(b1) an overshooting K exists, >100 mHa at K=4 on H4/N2; (b2) genuine escapes (>1e-6 Ha),
    worst >1 mHa; (b3) premise violation without escape exists -- so INVALID is decided by
    containment, not by the premise. (b2)'s pre-registered "every escape at K <= 8" was KILLED by
    full-space LiH (M=6, K=12-16; spec section 5): pinned both ways so neither half moves silently."""
    for key in ("h4", "n2", "lih"):
        mp, s = _map(key)
        assert s["overshoot_k"], key
        assert s["escapes"] and s["worst_mha"] > 1e-3, (key, s["worst_mha"])
    for key in ("h4", "n2", "lih_cas"):
        assert all(k <= 8 for _, k in _map(key)[1]["escapes"]), key
    assert any(k > 8 for _, k in _map("lih")[1]["escapes"]), _map("lih")[1]["escapes"]
    for key in ("h4", "n2"):
        mp, s = _map(key)
        assert (mp["odmd"][4] - mp["e1"]) * 1e3 > 100.0, (key, mp["odmd"][4] - mp["e1"])
        assert s["degraded"], key
    assert max(_map(key)[1]["worst_mha"] for key in ("h4", "n2", "lih")) > 1.0


def test_Gb_oracle_free_flag_sound_and_partial():
    """(b4) premise_refuted (eps > theta_1(M)) never fires unless eps > E_1 + TOL, never in self
    mode. The pre-registered "fires on every M >= 4 escape, on no M = 2 escape" was KILLED both
    ways (spec section 5): coverage is set by eps - theta_1(M), not by M, and on every system with
    escapes the flag both catches and misses some -- a silent flag never certifies the bracket."""
    for key in SYSTEMS:
        mp, s = _map(key)
        for cell, br in mp["cells"].items():
            assert not br.premise_refuted or br.eps > mp["e1"] + TOL, (key, cell)
        assert not any(br.premise_refuted for br in mp["self"].values()), key
        caught = [mp["cells"][c].premise_refuted for c in s["escapes"]]
        assert not caught or (any(caught) and not all(caught)), (key, s["escapes"], caught)
    mh = build_molecular_hamiltonian(**SYSTEMS["h2"])   # shot noise: theta_1 is no bound -> off
    noisy = QuantumKrylovSolver(mh, noise_sigma=1e-3, seed=0)
    assert krylov_bracket(mh, 2, eps=10.0).premise_refuted
    assert not krylov_bracket(mh, 2, eps=10.0, solver=noisy).premise_refuted


def test_Gc_free_region_does_not_repair():
    """At M=4 on H4/N2 the only free ODMD depth (K=4 <= M) escapes where self mode still contains
    E_0, and the first repairing K exceeds 2M = 8."""
    for key in ("h4", "n2"):
        mp, s = _map(key)
        assert (4, 4) in s["escapes"] and s["self_m4_contains"], (key, s["escapes"])
        assert s["repair_k_m4"] and min(s["repair_k_m4"]) > 8, (key, s["repair_k_m4"])
