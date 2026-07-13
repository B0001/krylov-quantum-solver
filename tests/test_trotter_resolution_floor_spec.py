"""
Acceptance gates G1-G4 for specs/SPEC_trotter_resolution_floor.md.

Claim 1 (the fix): Trotter synthesis was nondeterministic -- SuzukiTrotter products depend on
Pauli term order, which varied with hash randomization. build_trotter_step now orders terms
canonically (largest |coeff| first), pinning every derived number.

Claim 2 (the law): a circuit eigenphase is resolvable only if the reference population on it
exceeds the leakage floor ||U_trot - U_exact||^2. Below the floor, the extracted branch is an
artifact of the ordering -- probed by re-synthesizing under different orderings.

Origin: the nb3x8_device_gap G2 coin-flip flake (see BACKLOG). PySCF/qiskit, no block2;
`make gates` runs this in its own process.
"""
import numpy as np

from nb3x8_device_gap import circuit_gap, exact_gap, sector_models
from nb3x8_gaps import NB3X8_LT_BULK
from trotter_odmd import select_ground_eigenphase
from trotter_resolution_floor import (
    is_resolvable,
    leakage_floor,
    ordering_spread,
    reference_population,
)

F8 = NB3X8_LT_BULK["Nb3F8"]


def _sector2():
    return sector_models(**F8)[2].to_hamiltonian()


# --- G1: determinism -----------------------------------------------------------------------


def test_G1_trotter_step_is_order_canonical():
    """build_trotter_step yields the SAME unitary for shuffled input term orders, and the
    Nb3F8 sector-2 reps=1 deviation is pinned at the coeff-desc value (5.96e-3), not the
    1.04e-2 alternative the hash seed used to deal."""
    from qiskit.quantum_info import Operator, SparsePauliOp

    from hybrid_quantum_solver.trotter_krylov import build_trotter_step

    mh = _sector2()
    op = mh.qubit_hamiltonian
    rng = np.random.default_rng(0)
    terms = [(p.to_label(), complex(c)) for p, c in zip(op.paulis, op.coeffs)]
    ref = None
    for _ in range(3):
        rng.shuffle(terms)
        shuffled = SparsePauliOp.from_list(terms)
        U = Operator(build_trotter_step(shuffled, 0.1, order=2, reps=1)).data
        if ref is None:
            ref = U
        assert np.linalg.norm(U - ref, 2) < 1e-12

    dev = np.sqrt(leakage_floor(mh, reps=1))
    assert abs(dev - 5.96e-3) < 5e-4, dev  # pinned to coeff-desc, not the 1.04e-2 ordering


def test_G1_floor_shrinks_as_reps_squared():
    """The leakage floor follows the order-2 law: floor(reps) ~ reps^-4 (dev ~ reps^-2)."""
    mh = _sector2()
    f1, f2 = leakage_floor(mh, reps=1), leakage_floor(mh, reps=2)
    ratio = f1 / f2  # dev ratio ~4 -> floor ratio ~16
    assert 10.0 < ratio < 40.0, ratio


# --- G2: the floor separates the flaky assertion from the stable ones ----------------------


def test_G2_floor_criterion_matches_the_flake():
    """Nb3F8 sector-2: BELOW the floor at reps=1 (the assertion that coin-flipped), ABOVE it
    at reps=2 (the assertions that never flaked). The population is genuine (1.4e-5), the
    floor moves."""
    mh = _sector2()
    pop = reference_population(mh)
    assert 1e-6 < pop < 1e-4, pop  # the genuine 1.36e-5, not the leakage-corrupted 1e-13

    assert not is_resolvable(mh, reps=1)
    assert is_resolvable(mh, reps=2)
    assert leakage_floor(mh, reps=1) > pop > leakage_floor(mh, reps=2)


def test_G2_ordering_probe_flips_only_below_the_floor():
    """The probe: re-synthesize the F8 gap under different canonical orderings. Below the
    floor (reps=1) the branch flips -- spread ~ the 2pi/tau wrap quantum (>1000 meV). Above
    it (reps=2) orderings differ by ordinary Trotter bias (<10 meV)."""
    assert ordering_spread(F8, reps=1) > 1000.0
    assert ordering_spread(F8, reps=2) < 10.0


# --- G3: the capstone gate is deterministically green ---------------------------------------


def test_G3_f8_circuit_gap_is_deterministic_and_correct():
    """circuit_gap(F8, reps=1) was a coin flip (-1171.1 on ~half the runs); with canonical
    ordering it returns the historical green-run value +2580.70 every time. In-process
    repeats catch nondeterministic synthesis; cross-process nondeterminism is G1's job."""
    ref = exact_gap(**F8)
    vals = [circuit_gap(**F8, reps=1) for _ in range(5)]
    assert max(vals) - min(vals) < 1e-9, vals
    assert abs(vals[0] - ref) < 1.0, (vals[0], ref)


# --- G4: no committed Trotter number moves --------------------------------------------------


def test_G4_device_gap_recorded_numbers_reproduced():
    """The canonical order reproduces the committed device-gap numbers exactly -- the fix
    pins the green-run history rather than moving it: I8 reps=1/2 biases keep the recorded
    order-2 ratio in (3.3, 5.5), and F8 reps=1 lands < 1 meV from exact (the original G2
    assertions, now deterministic)."""
    i8 = NB3X8_LT_BULK["Nb3I8"]
    ref = exact_gap(**i8)
    b1 = circuit_gap(**i8, reps=1) - ref
    b2 = circuit_gap(**i8, reps=2) - ref
    assert abs(b1) > 50.0, b1
    assert 3.3 < b1 / b2 < 5.5, b1 / b2
    assert abs(circuit_gap(**F8, reps=1) - exact_gap(**F8)) < 1.0


# --- G5: unphysical-branch hardening (backlog latent item) ---------------------------------


def test_G5_selection_rejects_out_of_band_branch_even_at_higher_population():
    """The latent bug this hardens: `min(-angle(lam)/tau)` over the population cut alone would
    prefer a periodic image outside the physical band [-width/2, width/2] if it ever had lower
    (i.e. more negative) angle AND cleared the population cut -- exactly the F8 sector-2
    scenario noted in BACKLOG (an e = -1657.45 meV image, |e| > width/2). Synthetic here because
    it is untriggered by every committed system: a genuine in-band candidate at pop just above
    cut must beat a lower, higher-population out-of-band decoy."""
    width = 100.0
    angles = np.array([-200.0, -40.0, 30.0, 55.0])  # first is a decoy outside +-width/2 = 50
    pops = np.array([0.9, 1e-7, 0.5, 1e-6])
    e = select_ground_eigenphase(angles, pops, width, pop_cut=1e-8)
    assert e == -40.0, e  # not -200.0: the decoy is excluded despite pop=0.9 > cut


def test_G5_selection_keeps_the_closed_boundary_despite_roundoff():
    """A candidate that IS the band edge (|e| == width/2 up to float roundoff) must still be
    selected -- the tolerance in select_ground_eigenphase exists so ULP noise from the
    eig/angle chain cannot silently drop a genuine boundary branch (this is what a naive strict
    `<= width/2` filter breaks on I8 sector-1 reps=2, see BACKLOG)."""
    width = 436.4  # width/2 = 218.2
    angles = np.array([-218.20000000000002, 218.19999999999996])
    pops = np.array([0.5, 0.5])
    e = select_ground_eigenphase(angles, pops, width, pop_cut=1e-8)
    assert e == -218.20000000000002, e


def test_G5_f8_i8_eigenphases_unchanged_under_physical_band_selection():
    """The hardening is pure -- every committed circuit_gap number reproduces exactly through
    select_ground_eigenphase, for both the systems that previously exercised the coin-flip
    (F8 sector-2 reps=1/2) and the deterministic ones (I8)."""
    i8 = NB3X8_LT_BULK["Nb3I8"]
    assert abs(circuit_gap(**F8, reps=1) - 2580.70) < 0.5
    assert abs(circuit_gap(**i8, reps=1) - exact_gap(**i8)) > 50.0
    ref = exact_gap(**i8)
    b1 = circuit_gap(**i8, reps=1) - ref
    b2 = circuit_gap(**i8, reps=2) - ref
    assert 3.3 < b1 / b2 < 5.5, b1 / b2
