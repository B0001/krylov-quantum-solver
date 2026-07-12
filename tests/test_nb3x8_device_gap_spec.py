"""
Acceptance gates G1-G4 for specs/SPEC_nb3x8_device_gap.md (capstone composition study).

Test-first: ``nb3x8_device_gap`` does not exist yet, so this file is RED until the spec is
implemented. Claim: the validated materials thread (Nb3X8 downfolded dimer clusters, exact ED
gaps) and the device-validated ODMD stack compose end-to-end -- a material cluster's charge gap
Delta = E(3) + E(1) - 2 E(2) measured through genuinely-Trotterized Hadamard-test circuits under
an Aer device noise model, one depolarizing-immune ground-state ODMD run per particle sector,
with the Trotter bias Richardson-removed.

Recorded findings gated here: sector Trotter biases do NOT cancel in the gap (reps=1 leaves 35%
on Nb3I8, while near-commuting Nb3F8 is at 0.1 meV); and the bias-vs-noise crossover -- at
cx=1e-4 Richardson pays, at cx=1e-3 it stops (noise floor above the reps=2 bias).

REVISED (specs/SPEC_trotter_resolution_floor.md): the original recorded values (12% bias,
crossover at cx=3e-4) were measured under hash-ordering-nondeterministic Trotter synthesis and
are unreproducible under ANY canonical ordering. With deterministic (coeff-desc) ordering the
I8 reps=1 bias is -292.9 meV (35%), ratio 4.29, and Richardson still pays at cx=3e-4 and 6e-4
(16.7 < 60.7, 39.0 < 43.8 meV) -- it stops at cx=1e-3 (154.1 > 39.8).

References: fixed_filling_energy sector FCI + SPEC_nb3x8_gaps' recorded 842.44 meV. Aer is
seeded and Trotter synthesis is canonically ordered, so the medians are deterministic. Energies
in meV. PySCF/qiskit-aer, no block2; `make gates` runs it in its own process.
"""
import numpy as np

from hybrid_quantum_solver.noise import build_depolarizing_noise_model
from nb3x8_gaps import NB3X8_LT_BULK
from nb3x8_device_gap import (
    circuit_gap,
    device_gap,
    device_gap_richardson,
    exact_gap,
    statevector_gap,
)

I8 = NB3X8_LT_BULK["Nb3I8"]
SHOTS = 32768


def _noise(cx):
    return build_depolarizing_noise_model(cx / 10, cx, cx)


def test_G1_statevector_pipeline_is_exact():
    """ODMD sector gaps == sector FCI (< 0.01 meV) on all four LT-bulk materials, and the Nb3I8
    reference reproduces SPEC_nb3x8_gaps' recorded 842.44 meV (the cross-spec pin)."""
    for name, p in NB3X8_LT_BULK.items():
        ref = exact_gap(**p)
        assert abs(statevector_gap(**p) - ref) < 0.01, name
    assert abs(exact_gap(**I8) - 842.44) < 0.5


def test_G2_sector_trotter_biases_do_not_cancel():
    """reps=1 leaves > 50 meV (35% under canonical ordering) on the Nb3I8 gap -- the sector
    biases do NOT cancel -- with the order-2 ratio on reps doubling; near-commuting Nb3F8 sits
    < 1 meV (the contrast). This assertion coin-flipped before canonical Trotter ordering: the
    F8 reps=1 gap sits BELOW the Trotter resolution floor, and the branch depended on the term
    order the hash seed dealt (specs/SPEC_trotter_resolution_floor.md)."""
    ref = exact_gap(**I8)
    b1 = circuit_gap(**I8, reps=1) - ref
    b2 = circuit_gap(**I8, reps=2) - ref
    assert abs(b1) > 50.0, b1
    assert 3.3 < b1 / b2 < 5.5, b1 / b2
    f8 = NB3X8_LT_BULK["Nb3F8"]
    assert abs(circuit_gap(**f8, reps=1) - exact_gap(**f8)) < 1.0


def test_G3_richardson_fixes_the_circuit_gap():
    """Reps-(2,4) Richardson: < 1 meV residual and > 5x below the reps=4 bias; (1,2) < 10 meV."""
    ref = exact_gap(**I8)
    g2, g4 = circuit_gap(**I8, reps=2), circuit_gap(**I8, reps=4)
    rich24 = (4.0 * g4 - g2) / 3.0
    assert abs(rich24 - ref) < 1.0, rich24 - ref
    assert abs(rich24 - ref) < abs(g4 - ref) / 5.0
    g1 = circuit_gap(**I8, reps=1)
    assert abs((4.0 * g2 - g1) / 3.0 - ref) < 10.0


def test_G4_device_measurement_and_crossover():
    """DEFINITION OF DONE: at cx=1e-4 the Richardson device gap lands < 15 meV (1.2%) and > 5x
    below the raw reps=1 device gap; at cx=1e-3 the noise floor exceeds the reps=2 bias and
    Richardson stops paying -- SPEC_trotter_odmd R1 demonstrated on a material.

    REVISED: the crossover was recorded at cx=3e-4 under hash-nondeterministic Trotter ordering.
    Deterministically (canonical ordering) Richardson still pays at 3e-4 AND 6e-4 (16.7 < 60.7,
    39.0 < 43.8 meV); it stops at 1e-3 (154.1 > 39.8). The crossover EXISTS either way -- its
    location was never reproducible before specs/SPEC_trotter_resolution_floor.md."""
    ref = exact_gap(**I8)
    nm = _noise(1e-4)
    rich = [abs(device_gap_richardson(**I8, shots=SHOTS, noise_model=nm, seed=sd) - ref)
            for sd in range(5)]
    plain1 = [abs(device_gap(**I8, shots=SHOTS, noise_model=nm, seed=sd, trotter_reps=1) - ref)
              for sd in range(5)]
    assert np.median(rich) < 15.0, np.median(rich)
    assert np.median(plain1) / np.median(rich) > 5.0, (np.median(plain1), np.median(rich))
    # Richardson still pays in the mid-noise regime under deterministic ordering ...
    nm3 = _noise(3e-4)
    rich3 = [abs(device_gap_richardson(**I8, shots=SHOTS, noise_model=nm3, seed=sd) - ref)
             for sd in range(5)]
    plain3 = [abs(device_gap(**I8, shots=SHOTS, noise_model=nm3, seed=sd, trotter_reps=2) - ref)
              for sd in range(5)]
    assert np.median(rich3) < np.median(plain3), (np.median(rich3), np.median(plain3))
    # ... and stops paying once the noise floor clears the reps=2 bias.
    nm10 = _noise(1e-3)
    rich10 = [abs(device_gap_richardson(**I8, shots=SHOTS, noise_model=nm10, seed=sd) - ref)
              for sd in range(5)]
    plain10 = [abs(device_gap(**I8, shots=SHOTS, noise_model=nm10, seed=sd, trotter_reps=2) - ref)
               for sd in range(5)]
    assert np.median(rich10) > np.median(plain10), (np.median(rich10), np.median(plain10))
