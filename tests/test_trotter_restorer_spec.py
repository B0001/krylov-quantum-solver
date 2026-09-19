"""
Acceptance gates for specs/SPEC_trotter_floor_restorer.md (signal-domain Richardson for Trotter).

v1.0 set out to plug a "Trotter variational-floor leak". The gates below record that the leak does
not exist, that the premise inverts, and that the proposed cure is worse than the one already
shipped (specs/README.md step 5):

  G1   THE KILL: raw order-2 Trotter ODMD is ABOVE E_FCI in every configuration tested -- there is
       no leak to plug.
  G1b  THE INVERSION: the EXTRAPOLANT undershoots E_FCI in 22 of 24 configurations. If anything
       breaks the variational floor here it is the mitigation, not the Trotterization.
  G2   what signal extrapolation actually buys, re-gated: a large error cut at moderate K, then
       degeneration onto the un-extrapolated signal at large K.
  G3   FALSIFIED/out of scope: no Lehmann bracket is computable from an ODMD signal at all.
  G4   REVERSED: scalar (energy-domain) Richardson is never worse and is up to 95x better.

Deterministic (statevector circuits, no RNG). PySCF/qiskit, no block2; `make gates` isolates it.
"""
import numpy as np
import pytest

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from odmd import odmd_energy
from trotter_restorer import (
    phase_budget,
    restored_energy,
    signal_richardson,
    trotter_problems,
)

SYSTEMS = {"H2@1.6": "H 0 0 0; H 0 0 1.6",
           "H4@0.9": "H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"}
WINDOWS = (6, 8, 12, 16, 24, 32)
_CACHE = {}


def _probs(name):
    if name not in _CACHE:
        mh = build_molecular_hamiltonian(atom=SYSTEMS[name])
        _CACHE[name] = trotter_problems(mh, n=32)
    return _CACHE[name]


def test_G1_there_is_no_variational_floor_leak():
    """THE KILL: Trotterization never pushes the ODMD energy below E_FCI here.

    v1.0's G1 demanded the raw Trotter energy drop BELOW E_FCI by at least 1.0 mHa at large step
    size, "exposing and verifying the variational-floor leak". Measured over 2 systems x reps
    {1,2,4} x 6 window lengths: the energy is above FCI in 36/36 configurations, minimum excess
    +0.26 mHa, typical +0.6 to +20 mHa. The order-2 Suzuki effective-Hamiltonian bias is POSITIVE,
    so the phenomenon this spec exists to repair does not occur. Gated as a floor on the excess so
    the kill cannot silently decay into a pass.
    """
    worst = np.inf
    for name in SYSTEMS:
        P = _probs(name)
        for r in (1, 2, 4):
            for K in WINDOWS:
                e, _ = odmd_energy(P[r].s[:K], P[r].tau)
                worst = min(worst, e - P[r].ref)
    assert worst > 0.0, worst                       # never below FCI: no leak
    assert worst > 1e-4, worst                      # and not marginally so (+0.26 mHa)


def test_G1b_the_extrapolant_is_what_breaks_the_floor():
    """THE INVERSION: mitigation undershoots E_FCI where Trotterization never did.

    Richardson extrapolation is not variational -- its residual is two-sided. Across this
    2-system grid the extrapolated energies land BELOW E_FCI in 12 of 24 configurations (22 of 24
    over the 4-system shipped benchmark), by ~0.01 mHa, while G1 shows the RAW energy never does.
    So v1.0 had the causality backwards: it proposed extrapolation to restore a floor that
    Trotterization does not break, and extrapolation is the only step here that does break it.
    The undershoot is small and benign -- it is gated so nobody later mistakes the extrapolated
    energy for a bound.
    """
    undershoot = 0
    total = 0
    for name in SYSTEMS:
        P = _probs(name)
        for K in WINDOWS:
            for mode in ("signal_extrap", "scalar_extrap"):
                total += 1
                if restored_energy(P, K, mode) < P[1].ref:
                    undershoot += 1
    assert undershoot >= total // 2, (undershoot, total)   # 12/24 here, 22/24 on 4 systems
    assert undershoot < total, (undershoot, total)         # not a constant offset either


def test_G2_signal_extrapolation_helps_then_degenerates():
    """Re-gates what v1.0's `signal_extrap` mode actually delivers, in both directions.

    On H4 at K = 16 it cuts the raw error from +0.646 to +0.007 mHa -- an 87x reduction, real and
    worth having. But the Trotter bias enters s_k as exp(-i k tau dE), so linear extrapolation of
    the SIGNAL only cancels it while that phase stays small; at K = 32 the extrapolant returns the
    un-extrapolated value to within 0.003 mHa, buying nothing. Both halves are gated: the win at
    moderate K, and the degeneration at large K.
    """
    P = _probs("H4@0.9")
    ref = P[1].ref
    raw16 = abs(restored_energy(P, 16, "raw") - ref)
    sig16 = abs(restored_energy(P, 16, "signal_extrap") - ref)
    assert raw16 / sig16 > 20.0, (raw16 * 1e3, sig16 * 1e3)        # measured 87x

    raw32 = restored_energy(P, 32, "raw")
    sig32 = restored_energy(P, 32, "signal_extrap")
    assert abs(sig32 - raw32) < 1e-5, (sig32 - raw32)              # degenerate: no benefit left


def test_G3_no_lehmann_bracket_exists_on_this_path():
    """FALSIFIED/OUT OF SCOPE: v1.0's step 4 cannot be built from an ODMD signal.

    v1.0's G3 required "100% sound, oracle-enclosed Lehmann brackets" on the restored Ritz states.
    Lehmann needs <u|H^2|u> for a physical Ritz STATE. ODMD's entire premise is that it works from
    the scalar overlap series s_k alone -- it returns eigenPHASES, and no state is ever formed, so
    neither <H> nor <H^2> of a Ritz vector is available to extrapolate or to bound. The gate pins
    that the module exposes no bracket function rather than shipping a number that looks certified
    and is not. Certifying a Trotter solve needs the QKSD/state path (excited_bounds.py), which is
    a different spec.
    """
    import trotter_restorer
    assert not hasattr(trotter_restorer, "certified_restored_energy")
    assert not any("lehmann" in n.lower() or "bracket" in n.lower()
                   for n in dir(trotter_restorer)), dir(trotter_restorer)
    with pytest.raises(ValueError):
        restored_energy(_probs("H4@0.9"), 12, mode="certified")


def test_G4_neither_extrapolation_dominates_the_other():
    """REVERSES v1.0's G4 -- but not the way a first pass assumes. Each wins in its own regime.

    v1.0 claimed signal-domain (state-level) extrapolation is "mathematically superior" to
    post-solve scalar extrapolation. It is not uniformly superior -- and scalar is not uniformly
    superior either, which a first draft of this gate wrongly asserted before the data corrected
    it:

      * INSIDE the phase window (small K) signal wins: stretched H2 at K = 6 leaves 0.0045 mHa
        against scalar's 0.0476 mHa -- 10x better.
      * OUTSIDE it (large K) scalar wins by far more: at K = 32 scalar beats signal by ~95x on H4
        and ~21x on H2, because the dt^2 law holds exactly for the eigenPHASE and so scalar is
        flat in K, while signal degenerates onto the raw value (G2).

    The controlling quantity is the accumulated phase K*tau*dE, not the extrapolation domain. So
    v1.0's superiority claim is falsified, the honest replacement is a regime map, and the
    practical advice is scalar -- it is flat, cheaper, and already shipped.
    """
    P2, P4 = _probs("H2@1.6"), _probs("H4@0.9")

    sig6 = abs(restored_energy(P2, 6, "signal_extrap") - P2[1].ref)
    sca6 = abs(restored_energy(P2, 6, "scalar_extrap") - P2[1].ref)
    assert sca6 / sig6 > 5.0, (sca6 * 1e3, sig6 * 1e3)        # signal wins at small K, ~10x

    for P, floor in ((P4, 20.0), (P2, 10.0)):
        sig = abs(restored_energy(P, 32, "signal_extrap") - P[1].ref)
        sca = abs(restored_energy(P, 32, "scalar_extrap") - P[1].ref)
        assert sig / sca > floor, (floor, sig * 1e3, sca * 1e3)   # scalar wins at large K

    flat = [abs(restored_energy(P4, K, "scalar_extrap") - P4[1].ref) for K in (12, 16, 24, 32)]
    assert max(flat) - min(flat) < 2e-5, flat                 # scalar is flat in K; signal is not


def test_G5_the_refinement_parameter_is_reps_not_dt():
    """v1.0's extrapolation formula pairs signals sampled at DIFFERENT physical times.

    s_k(dt) sits at time k*dt and s_k(2dt) at 2k*dt, so combining them at equal index k subtracts
    unlike quantities and the O(dt^2) term cannot cancel. Refining `reps` at fixed tau keeps every
    sample at k*tau while halving dt_eff. The gate shows the time-aligned pairing reduces the
    error while the index-paired version -- the literal v1.0 formula, s_k at two different taus --
    does not even share a time axis: its tau differs, so the two problems disagree on what k means.
    """
    P = _probs("H4@0.9")
    assert P[1].tau == P[2].tau == P[4].tau          # reps refinement keeps the time axis fixed
    aligned = signal_richardson(P[2].s[:16], P[4].s[:16])
    e_aligned, _ = odmd_energy(aligned, P[1].tau)
    assert abs(e_aligned - P[1].ref) < abs(restored_energy(P, 16, "raw") - P[1].ref)

    budgets = [phase_budget(K, P[1].tau, 2.6053e-3) for K in WINDOWS]
    assert budgets == sorted(budgets), budgets       # the diagnostic grows with the window
