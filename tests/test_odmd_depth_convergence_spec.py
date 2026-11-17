"""
Acceptance gates G1-G3 for specs/SPEC_odmd_depth_convergence.md (chem-36m, follow-up to chem-z3h).

Claim under test: (1) certified_gaps.GapBracket.premise_refuted is a sound oracle-free proof that
the Weinstein/self premise eps1 <= E_1 failed, mirroring temple_bounds.EnergyBracket's flag; (2) a
depth-convergence check on ODMD's eps(K) -- stability over a trailing window of already-measured K,
no look-ahead -- never green-lights a K that the gated (M, K) map shows produces a containment
escape, but is NOT sufficient: it can be satisfied at a K where the premise is still violated
(a stable-but-wrong fixed point).

Reuses the already-gated map from tests/test_odmd_temple_repair_map_spec.py's module
(odmd_temple_repair_map.k_vs_validity_map) -- no new Krylov/ODMD machinery, same reference (exact
reachable E_0/E_1, an oracle used only to check the flag/the check, never to compute them).

Noiseless; PySCF/qiskit, no block2; `make gates` runs this in its own process.
"""
from functools import lru_cache

import numpy as np

from certified_gaps import gap_bracket
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from odmd_depth_convergence import TOL, converged_depths, depth_converged, safe_eps1
from odmd_temple_repair_map import K_DEPTHS, SYSTEMS, k_vs_validity_map, summarize

# The four non-trivial systems from SPEC_odmd_temple_repair_map (H2 has <= 4 reachable levels and
# no escapes anywhere on the grid -- nothing for G2 to check).
SYSTEM_KEYS = ("h4", "n2", "lih_cas", "lih")


@lru_cache(maxsize=None)
def _map(key):
    mp = k_vs_validity_map(key)
    return mp, summarize(mp)


def test_G1_gap_bracket_premise_refuted_sound():
    """GapBracket.premise_refuted mirrors EnergyBracket's: never fires without eps1 > E1 + TOL
    against the oracle; never fires in self mode; a deliberately-wrong oracle e1 (set above
    theta1) does fire it; under shot noise theta1 is left at +inf so the flag stays False even
    though eps1 is then unverifiable (soundness over usefulness)."""
    for key in SYSTEM_KEYS:
        mp, _ = _map(key)
        e1 = mp["e1"]
        mh = build_molecular_hamiltonian(**SYSTEMS[key])
        solver = QuantumKrylovSolver(mh)
        for m in (2, 4, 6, 8):
            br = gap_bracket(mh, m, solver=solver)
            assert not br.premise_refuted or br.eps1 > e1 + TOL, (key, m, br.eps1, e1)
            if br.eps1_source == "self":
                assert not br.premise_refuted, (key, m, "self mode must never trip the flag")

    # A deliberately-wrong oracle e1, set strictly above theta1, MUST be caught.
    mh = build_molecular_hamiltonian(**SYSTEMS["h4"])
    solver = QuantumKrylovSolver(mh)
    br0 = gap_bracket(mh, 8, solver=solver)
    assert np.isfinite(br0.theta1)
    bad_e1_total = (br0.theta1 + 1.0) + mh.energy_offset
    br_bad = gap_bracket(mh, 8, e1=bad_e1_total, solver=solver)
    assert br_bad.premise_refuted

    # Under shot noise theta1 is no longer a variational bound -> left at +inf -> flag off, even
    # though eps1 (the same bad oracle value) is now unverifiable.
    noisy = QuantumKrylovSolver(mh, noise_sigma=1e-3, seed=0)
    br_noisy = gap_bracket(mh, 8, e1=bad_e1_total, solver=noisy)
    assert not np.isfinite(br_noisy.theta1) and br_noisy.theta1 > 0
    assert not br_noisy.premise_refuted


def test_G2_depth_convergence_never_coincides_with_an_escape():
    """DEFINITION OF DONE: no K that depth_converged declares True, for any system in
    SYSTEM_KEYS, coincides with any (M, K) the gated map records as a containment escape (lower
    > E0 + TOL) at that K, for ANY base depth M on the grid. If this fails, requiring
    depth-convergence before substitution buys nothing over using eps(K) raw."""
    any_escapes_checked = False
    for key in SYSTEM_KEYS:
        mp, s = _map(key)
        escape_ks = {k for _, k in s["escapes"]}
        assert escape_ks, key  # sanity: every one of these four systems has >= 1 escape K
        any_escapes_checked = True
        conv = set(converged_depths(mp["odmd"], K_DEPTHS))
        overlap = conv & escape_ks
        assert not overlap, (key, "converged-but-escaping K's", overlap, "converged:", conv,
                             "escaping:", escape_ks)
    assert any_escapes_checked


def test_G3_convergence_is_necessary_not_sufficient():
    """The honest boundary: depth_converged can be True at a K where the premise is still
    violated against the oracle -- a stable-but-wrong fixed point, not a transient. Measured on
    N2 at K in {26, 28} (eps(K) - E1 ~ +34...+35 mHa, stable within ~1 mHa over the trailing
    window). Convergence alone must never be read as a correctness proof."""
    mp, s = _map("n2")
    e1 = mp["e1"]
    conv = converged_depths(mp["odmd"], K_DEPTHS)
    violated_and_converged = [k for k in conv if mp["odmd"][k] is not None
                              and mp["odmd"][k] > e1 + TOL]
    assert violated_and_converged, (conv, {k: mp["odmd"][k] - e1 for k in conv})


def test_safe_eps1_abstains_unless_converged():
    """safe_eps1 returns None (abstain) at a shallow, unconverged K, and a finite value once the
    trailing window has stabilized -- the function other code should call instead of
    odmd_temple_repair_map.odmd_eps1 directly when substituting into krylov_bracket/gap_bracket."""
    mh = build_molecular_hamiltonian(**SYSTEMS["h4"])
    solver = QuantumKrylovSolver(mh)
    assert safe_eps1(mh, solver, 4) is None          # shallow: too few points in the window
    deep_val = safe_eps1(mh, solver, 22)              # H4's one genuinely stable repair point
    assert deep_val is not None and np.isfinite(deep_val)


def test_depth_converged_requires_full_resolved_window():
    """depth_converged is False when fewer than `window` grid points are available at or below
    k, and False when any point in the trailing window is an unresolved ODMD fit (None) --
    both must hold regardless of how close the resolved values happen to be."""
    odmd = {4: -1.0, 5: None, 6: -1.0005, 7: -1.0003, 8: -1.0001}
    assert not depth_converged(odmd, (4, 5, 6, 7, 8), 4, window=3, tol=2e-3)   # too few points
    assert not depth_converged(odmd, (4, 5, 6, 7, 8), 6, window=3, tol=2e-3)   # window hits None
    assert depth_converged(odmd, (4, 5, 6, 7, 8), 8, window=3, tol=2e-3)       # {6,7,8}: resolved, tight
