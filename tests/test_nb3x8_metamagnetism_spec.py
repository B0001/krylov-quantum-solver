"""
Acceptance gates G1-G4 for specs/SPEC_nb3x8_metamagnetism.md (field-driven singlet->triplet
crossing).

Test-first: ``nb3x8_metamagnetism`` does not exist yet, so this file is RED until the spec is
implemented. Claim: adding a Zeeman term -h*Sz,tot to the exactly-diagonalizable Nb3X8 dimer
crosses its ground state from the non-magnetic singlet to the fully-polarized triplet member at
an exact closed form h_c = J (the same interlayer exchange gated in SPEC_odmd_spin), reproduced
by DIRECT diagonalization of the full field-augmented Hamiltonian (not assumed from symmetry).
Falsifiable: any numeric crossing field disagreeing with J, or a smeared/fractional magnetization
step, kills the claim. G4 ties the Tesla-scale prediction to two documented pulsed-field records.

PySCF/qiskit, no block2; `make gates` runs it in its own process.
"""
import numpy as np

from nb3x8_gaps import NB3X8_LT_BULK
from nb3x8_metamagnetism import (
    G_MU_B,
    critical_field_numeric,
    critical_field_tesla,
    magnetization,
    zeeman_ground_state,
)
from nb3x8_susceptibility import ionic_singlet_energy, n2_spectrum
from odmd_spin import dimer_exchange_analytic

MAGNETIC = ("Nb3Cl8", "Nb3Br8", "Nb3I8")  # Nb3F8 excluded: J below the model's own noise floor
ALL_MATERIALS = ("Nb3F8", "Nb3Cl8", "Nb3Br8", "Nb3I8")

# Documented pulsed-field records used only as a Tesla-scale feasibility yardstick (see spec
# Sources): non-destructive multi-shot (Los Alamos, March 2012) and destructive indoor
# electromagnetic flux-compression (U. Tokyo, 2018).
NONDESTRUCTIVE_RECORD_T = 100.75
DESTRUCTIVE_RECORD_T = 1200.0


def test_G1_closed_form_matches_direct_field_augmented_ed():
    """DEFINITION OF DONE: bisecting the full H0 - h*Sz matrix for its <Sz> crossing lands on the
    closed-form J to < 1e-6 meV -- no block-diagonal shortcut assumed."""
    for name in MAGNETIC:
        p = NB3X8_LT_BULK[name]
        J = dimer_exchange_analytic(**p)
        h_c = critical_field_numeric(**p)
        assert abs(h_c - J) < 1e-6, (name, h_c, J)


def test_G2_rigid_shift_identity():
    """The full field-augmented ground energy matches min(E_singlet, E_triplet - h) to < 1e-8 meV
    across a grid spanning both sides of the crossing, for every material (incl. Nb3F8)."""
    for name in ALL_MATERIALS:
        p = NB3X8_LT_BULK[name]
        J = dimer_exchange_analytic(**p)
        e_rel, _ = n2_spectrum(**p)          # energies relative to the h=0 ground state
        e0_abs, _ = zeeman_ground_state(**p, h=0.0)
        for h in np.linspace(0.0, 1.5 * J, 13):
            e_num, _ = zeeman_ground_state(**p, h=h)
            e_pred = e0_abs + min(e_rel[0], e_rel[1] - h)
            assert abs(e_num - e_pred) < 1e-8, (name, h, e_num, e_pred)


def test_G3_clean_magnetization_step():
    """<Sz> jumps from ~0 to ~1 within +-0.1% of J around the crossing -- no fractional plateau --
    and the first ionic singlet sits above J for every material, so no other level intervenes."""
    for name in MAGNETIC:
        p = NB3X8_LT_BULK[name]
        J = dimer_exchange_analytic(**p)
        eps = 1e-3 * J
        assert magnetization(**p, h=J - eps) < 0.01, name
        assert magnetization(**p, h=J + eps) > 0.99, name
        assert ionic_singlet_energy(**p) > J, name


def test_G4_tesla_scale_feasibility_boundary():
    """Cl/Br sit above the non-destructive record but below the destructive one; Nb3I8 exceeds
    even the destructive record. Nb3F8 is excluded (its J is below the noise floor, per spec)."""
    b_cl = critical_field_tesla(**NB3X8_LT_BULK["Nb3Cl8"])
    b_br = critical_field_tesla(**NB3X8_LT_BULK["Nb3Br8"])
    b_i = critical_field_tesla(**NB3X8_LT_BULK["Nb3I8"])

    for b in (b_cl, b_br):
        assert NONDESTRUCTIVE_RECORD_T < b < DESTRUCTIVE_RECORD_T, b
    assert b_i > DESTRUCTIVE_RECORD_T, b_i

    # sanity: B_c = J / (g*mu_B) is consistent with the closed-form J directly
    for name in MAGNETIC:
        p = NB3X8_LT_BULK[name]
        J = dimer_exchange_analytic(**p)
        assert abs(critical_field_tesla(**p) - J / G_MU_B) < 1e-9, name
