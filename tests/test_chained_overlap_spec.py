"""Gates for specs/SPEC_chained_overlap.md.

`krylov_refine.refine_via_lanczos` was a NotImplementedError stub since 2026-07-17. It chains the
Davis-Kahan bound through the Krylov ground Ritz vector, so the certificate uses v's residual rather
than the (much larger) HF residual that makes the direct SPEC-21 bound go vacuous on exactly the
multireference systems it is wanted for.

All references are built at the builder's default conv_tol (1e-9) -- NOT 1e-13, which an earlier
version of this docstring claimed. The gated systems are insensitive to that choice
(specs/SPEC_scf_conv_tol.md G1: on the <= 8-qubit members |dE0| <= 3.4e-14 Ha and |d overlap| <=
6.8e-9 between 1e-9 and 1e-13; linear H6 was not re-measured). a=1.10/1.35 are excluded because at
the default an SCF residue contaminates the reachable sector (specs/SPEC_reachability_tolerance.md)
and the certified target is not well defined.
"""
import sys

import numpy as np
import pytest

from certified_gaps import gap_bracket
from hf_overlap_certificate import REACHABLE_TOL_CERTIFIED, certify_hf_overlap
from hybrid_quantum_solver.certified_overlap.krylov_refine import (
    SATURATION_SLACK,
    refine_via_lanczos,
)
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

H2_EQ = "H 0 0 0; H 0 0 0.74"
H2_STRETCHED = "H 0 0 0; H 0 0 2.0"
H4_LINEAR = "H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"
# square H4 ONLY at its symmetric-SCF geometries -- elsewhere RHF breaks symmetry and there is no
# well-defined target (specs/SPEC_symmetry_reachability.md).
H4_SQUARE_105 = "H 0 0 0; H 1.05 0 0; H 1.05 1.05 0; H 0 1.05 0"
H6_LINEAR = "H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0; H 0 0 4.0; H 0 0 5.0"

# NOTE a=1.10 and a=1.35 are deliberately EXCLUDED. They are symmetric-SCF geometries and every
# reference in this file is built at the builder's default conv_tol (1e-9). At the default a=1.35
# carries the SCF-residue artifact (specs/SPEC_reachability_tolerance.md): its "exact reachable
# overlap" is the residue (7.8e-5), not the physical overlap (0.613 at conv_tol=1e-13). a=1.10
# carries it on the Linux freeze (residue 5.07e-10) but not on macOS (~1e-29); see G7 below.
# a=1.05 is clean at the default. The builder DOES accept conv_tol (an earlier version of this note
# said it could not), so re-admitting a=1.35 at TIGHT_SCF_CONV_TOL is possible -- a follow-up,
# specs/SPEC_scf_conv_tol.md section 7.

DIRECT_SURVIVES = (H2_EQ, H2_STRETCHED, H4_LINEAR)
DIRECT_VACUOUS = (H4_SQUARE_105, H6_LINEAR)
DEPTHS = (6, 8, 12)


def _setup(atom):
    mh = build_molecular_hamiltonian(atom=atom)
    H = mh.qubit_hamiltonian.to_matrix()
    u = np.asarray(mh.hf_state().data, dtype=complex)
    w, V = np.linalg.eigh(H)
    amp = np.abs(V.conj().T @ u)
    reach = np.where(amp ** 2 > 1e-10)[0]
    return mh, H, u, float(w[reach[1]]), float(amp[reach[0]])


def _ritz(mh, m, solver):
    v = np.asarray(solver.eigenstates(m, n_states=1)[1][0], dtype=complex)
    return v / np.linalg.norm(v)


def _both_modes(atom, m):
    """(chained_oracle, chained_self, exact) at depth m."""
    mh, H, u, e1_elec, exact = _setup(atom)
    solver = QuantumKrylovSolver(mh)
    v = _ritz(mh, m, solver)
    oracle = refine_via_lanczos(H, u, v, e1_elec)
    self_floor = gap_bracket(mh, m, e1=None, solver=solver).eps1
    self_mode = refine_via_lanczos(H, u, v, self_floor)
    return oracle, self_mode, exact


# --- G1: validity, both modes (DEFINITION OF DONE) ------------------------------------------------

@pytest.mark.parametrize("atom", DIRECT_SURVIVES + DIRECT_VACUOUS)
@pytest.mark.parametrize("m", DEPTHS)
def test_G1_chained_bound_never_exceeds_the_exact_overlap(atom, m):
    """The killable check. Slack is required, not cosmetic: at a machine-converged Ritz vector the
    bound SATURATES (equals the exact overlap) and rounding puts it 1-3 ulp either side."""
    oracle, self_mode, exact = _both_modes(atom, m)
    # Not vacuously satisfiable: an all-None implementation would pass a bare "if g is not None"
    # guard on every cell, so require a bound to actually exist here.
    assert oracle is not None, (atom, m, "oracle mode produced no bound")
    for label, g in (("oracle", oracle), ("self", self_mode)):
        if g is not None:
            assert g <= exact + SATURATION_SLACK, (label, atom, m, g, exact)


# --- G2: it strictly tightens the direct bound ----------------------------------------------------

@pytest.mark.parametrize("atom", DIRECT_SURVIVES)
@pytest.mark.parametrize("m", DEPTHS)
def test_G2_chained_beats_direct_where_direct_survives(atom, m):
    mh, H, u, e1_elec, _ = _setup(atom)
    solver = QuantumKrylovSolver(mh)
    direct = certify_hf_overlap(mh, m, e1=e1_elec + mh.energy_offset, solver=solver)
    assert not direct.vacuous, (atom, m)
    chained = refine_via_lanczos(H, u, _ritz(mh, m, solver), e1_elec)
    assert chained is not None and chained > direct.gamma_min, (atom, m, chained, direct.gamma_min)


# --- G3: it rescues the vacuous cases -- in SELF mode, not just oracle -----------------------------

@pytest.mark.parametrize("atom", DIRECT_VACUOUS)
@pytest.mark.parametrize("m", (8, 12))
def test_G3_chained_rescues_vacuous_direct_in_self_mode(atom, m):
    mh, H, u, e1_elec, exact = _setup(atom)
    solver = QuantumKrylovSolver(mh)
    direct = certify_hf_overlap(mh, m, e1=e1_elec + mh.energy_offset, solver=solver)
    assert direct.vacuous, (atom, m, "direct was expected VACUOUS here")
    _, self_mode, _ = _both_modes(atom, m)
    assert self_mode is not None and self_mode > 0.0, (atom, m)
    assert self_mode <= exact + SATURATION_SLACK, (atom, m, self_mode, exact)


# --- G4: self mode is nearly free, and absorbs a loose floor far better than direct does -----------

@pytest.mark.parametrize("atom", DIRECT_SURVIVES + DIRECT_VACUOUS)
@pytest.mark.parametrize("m", DEPTHS)
def test_G4_self_mode_costs_little_against_oracle(atom, m):
    oracle, self_mode, _ = _both_modes(atom, m)
    if oracle is None or self_mode is None:
        pytest.skip("vacuous in one mode -- covered by G3")
    assert self_mode / oracle >= 0.9, (atom, m, self_mode, oracle)


def test_G4_chained_absorbs_the_self_mode_floor_better_than_direct():
    """THE PRACTICAL ARGUMENT. A loose self-mode E_1 floor costs the direct bound ~38% on linear H4
    at M=6 but costs the chained bound under 1%, because the floor enters only through
    arcsin(r_v/delta_v) with a tiny r_v rather than through the much larger HF residual.
    """
    m = 6
    mh, H, u, e1_elec, _ = _setup(H4_LINEAR)
    solver = QuantumKrylovSolver(mh)
    d_oracle = certify_hf_overlap(mh, m, e1=e1_elec + mh.energy_offset, solver=solver).gamma_min
    d_self = certify_hf_overlap(mh, m, e1=None, solver=solver).gamma_min
    c_oracle, c_self, _ = _both_modes(H4_LINEAR, m)
    direct_loss = 1.0 - d_self / d_oracle
    chained_loss = 1.0 - c_self / c_oracle
    assert direct_loss > 0.30, direct_loss
    assert chained_loss < 0.01, chained_loss


# --- G5: NOT monotone in M ------------------------------------------------------------------------

def test_G5_chained_bound_is_not_monotone_in_krylov_depth():
    """Killed if the bound turns out monotone -- the caveat would be unnecessary and the spec should
    drop it rather than carry a false warning."""
    vals = [_both_modes(H4_SQUARE_105, m)[0] for m in (6, 8, 12)]
    assert all(v is not None for v in vals), vals
    assert vals[1] < vals[0], vals        # M=8 is WORSE than M=6 -- that alone is the claim
    # Deliberately NOT asserting vals[2] > vals[1]: the recovery shape is incidental, no spec claim
    # depends on it, and pinning it would break on unrelated solver changes.


# --- G6: vacuous is None, and inputs are checked --------------------------------------------------

def test_G6_vacuous_returns_none_not_a_number():
    mh, H, u, _, _ = _setup(H4_LINEAR)
    solver = QuantumKrylovSolver(mh)
    v = _ritz(mh, 6, solver)
    lam_v = float(np.real(np.vdot(v, H @ v)))
    assert refine_via_lanczos(H, u, v, lam_v - 1.0) is None      # floor below lambda_v


def test_G6_unnormalized_inputs_raise():
    mh, H, u, e1_elec, _ = _setup(H4_LINEAR)
    v = _ritz(mh, 6, QuantumKrylovSolver(mh))
    with pytest.raises(ValueError, match="normalized"):
        refine_via_lanczos(H, 2.0 * u, v, e1_elec)
    with pytest.raises(ValueError, match="normalized"):
        refine_via_lanczos(H, u, 0.5 * v, e1_elec)
    # The band that a np.isclose(atol=1e-8) guard silently admitted (its default rtol=1e-5 made it
    # ~1000x looser than it read). A 9e-6 denormalization pushes the bound ~1e-5 above the exact
    # overlap -- nine orders past SATURATION_SLACK -- so it must raise, not return a number.
    with pytest.raises(ValueError, match="normalized"):
        refine_via_lanczos(H, u, (1.0 + 9e-6) * v, e1_elec)


# --- G7: R2b (sector restriction inherited in NAME only) is measured, not asserted ----------------

# leak^2 is the term Davis-Kahan-on-v neglects. At 1e-12 the neglected term is 1e-24, ten orders
# below SATURATION_SLACK, so the bound cannot be affected. This is NOT a rounding tolerance -- the
# excluded geometries below sit eight orders on the other side of it.
LEAK_TOL = 1e-12


def _leakage(atom, m):
    """(||P_unreach v||, #unreachable levels below E1_reach) -- u's unreachable components vanish by
    construction, v's only numerically, which is the whole content of R2b."""
    mh = build_molecular_hamiltonian(atom=atom)
    H = mh.qubit_hamiltonian.to_matrix()
    u = np.asarray(mh.hf_state().data, dtype=complex)
    w, V = np.linalg.eigh(H)
    amp2 = np.abs(V.conj().T @ u) ** 2
    reach = np.where(amp2 > REACHABLE_TOL_CERTIFIED)[0]
    unreach = np.where(amp2 <= REACHABLE_TOL_CERTIFIED)[0]
    v = _ritz(mh, m, QuantumKrylovSolver(mh))
    leak = float(np.linalg.norm(V[:, unreach].conj().T @ v))
    return leak, int(np.sum(w[unreach] < float(w[reach[1]])))


@pytest.mark.parametrize("atom", DIRECT_SURVIVES + DIRECT_VACUOUS)
def test_G7_unreachable_leakage_is_negligible_on_the_gated_set(atom):
    leak, n_below = _leakage(atom, 6)
    # Non-vacuity: if no unreachable level sat below E1_reach there would be nothing to neglect and
    # R2b should be deleted rather than carried.
    assert n_below > 0, (atom, "no unreachable level below E1_reach -- R2b would be vacuous")
    assert leak < LEAK_TOL, (atom, leak)


# a=1.10 pins a platform-specific SCF stopping point, so it is Linux-reference-only. Measured on
# macOS 27 (Apple M3, scipy 1.15.3/Accelerate), default conv_tol=1e-9: the a=1.10 forbidden-level
# residue is ~1e-29 (Linux freeze: 5.07e-10), so the Ritz leakage is 8.9e-15 < LEAK_TOL and the
# positive control cannot fire; a=1.35 still does (1.6e-6). That is the finding, not a failure to fix:
# specs/SPEC_reachability_tolerance.md section 10. The a=1.35 case carries the claim everywhere.
_A110_LINUX_ONLY = pytest.mark.skipif(
    sys.platform != "linux",
    reason="a=1.10 residue is platform-dependent: 5.07e-10 on the Linux freeze, ~1e-29 on macOS 27 "
           "(leakage 8.9e-15 < LEAK_TOL=1e-12 there); a=1.35 still fires (1.6e-6) on every platform. "
           "specs/SPEC_reachability_tolerance.md section 10.",
)


@pytest.mark.parametrize("a", (pytest.param(1.10, marks=_A110_LINUX_ONLY), 1.35))
def test_G7_excluded_geometries_are_the_positive_control(a):
    """The exclusion is load-bearing for R2b, not only for R3. These geometries must BREACH the
    threshold -- if they passed it, excluding them would be unnecessary and R2b overstated."""
    atom = f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0"
    leak, _ = _leakage(atom, 6)
    assert leak > LEAK_TOL, (a, leak, "excluded geometry no longer breaches -- revisit R2b/R3")
    # And it breaches by enough to actually matter: leak^2 above SATURATION_SLACK.
    assert leak ** 2 > SATURATION_SLACK, (a, leak ** 2)
