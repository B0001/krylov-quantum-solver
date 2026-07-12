"""
Acceptance gates G1-G4 for specs/SPEC_lambda_meas_identity.md.

Claim: `precision_cost.measurement_lambda` summed every Pauli coefficient INCLUDING the
identity -- a zero-variance constant that costs zero shots -- so the published near-term shot
counts and crossovers were overstated by the identity fraction (30-46% of the 1-norm). Fixed in
place: identity excluded by default, `include_identity=True` kept for archaeology.

THE FINDING (G3): re-scoring is not a pure rescale -- H2's advisor verdict at rho=1e4 flips
FT -> near-term. The bias was hiding a regime where near-term is the right answer. (The BACKLOG
conjecture that no qualitative verdict would flip was wrong; recorded in the spec.)

PySCF + qiskit, no block2; `make gates` runs this in its own process.
"""
import numpy as np
import pytest
from pyscf import ao2mo, gto, mcscf, scf

from cost_advisor import advise
from hybrid_quantum_solver.molecular_hamiltonian import build_hamiltonian_from_integrals
from precision_cost import (
    crossover_epsilon,
    measurement_lambda,
    qubitization_lambda,
    resource_ratio,
)
from shift_both_sides import shot_lambda

EPS = 1.6e-3

MOLS = {
    "H2": ("H 0 0 0; H 0 0 0.74", 2, 2),
    "H2O(4,3)": ("O 0 0 0.117; H 0 0.757 -0.467; H 0 -0.757 -0.467", 3, 4),
    "N2(6,6)": ("N 0 0 0; N 0 0 1.10", 6, 6),
}

_CACHE: dict = {}


def _mol(name):
    """(mh, h1, eri, nelec, norb) for a molecule, cached across gates."""
    if name not in _CACHE:
        atom, norb, nelec = MOLS[name]
        mf = scf.RHF(gto.M(atom=atom, basis="sto-3g", verbose=0)).run()
        c = mcscf.CASCI(mf, norb, nelec)
        c.kernel()
        h1, e_core = c.get_h1eff()
        eri = ao2mo.restore(1, c.get_h2eff(), norb)
        na = (nelec + nelec % 2) // 2
        h1 = np.asarray(h1)
        ne = (na, nelec - na)
        mh = build_hamiltonian_from_integrals(h1, eri, ne, float(e_core))
        _CACHE[name] = (mh, h1, eri, ne, norb)
    return _CACHE[name]


# --- G1: the metric is honest by default ---------------------------------------------------


@pytest.mark.parametrize("name", list(MOLS))
def test_G1_identity_excluded_by_default(name):
    """measurement_lambda(mh) now excludes the identity (== shot_lambda); include_identity=True
    reproduces the old inflated value; and the identity fraction is material (> 25%) -- this
    was never a rounding error."""
    mh = _mol(name)[0]
    lam = measurement_lambda(mh)
    lam_old = measurement_lambda(mh, include_identity=True)

    assert lam == pytest.approx(shot_lambda(mh), rel=1e-12)
    assert lam_old == pytest.approx(float(np.abs(mh.qubit_hamiltonian.coeffs).sum()), rel=1e-12)
    assert (lam_old - lam) / lam_old > 0.25, (name, lam, lam_old)


# --- G2: the crossovers move by exactly the identity fraction (algebra check) ---------------


@pytest.mark.parametrize("name", list(MOLS))
def test_G2_crossovers_rescale_quadratically(name):
    """eps* and R are quadratic in lambda_meas, so old/honest ratios must equal
    (lam_old/lam_new)^2 -- ties the re-scored table to the algebra."""
    mh, h1, eri, ne, norb = _mol(name)
    lam, lam_old = measurement_lambda(mh), measurement_lambda(mh, include_identity=True)
    lam_df = qubitization_lambda(h1, eri, norb, nelec=ne, shift=True)
    q = (lam_old / lam) ** 2

    assert crossover_epsilon(lam_old, lam_df) / crossover_epsilon(lam, lam_df) == pytest.approx(q, rel=1e-9)
    assert resource_ratio(lam_old, lam_df, EPS) / resource_ratio(lam, lam_df, EPS) == pytest.approx(q, rel=1e-9)


# --- G3 (THE FINDING): one qualitative verdict flips ----------------------------------------


def test_G3_h2_verdict_flips_at_rho_1e4():
    """H2 at rho=1e4, chemical accuracy: the inflated metric said FT is cheaper; the honest
    one says near-term. The identity bias was hiding a regime where the near-term method is
    the right answer. H2O and N2 verdicts are stable across rho = 1..1e6."""
    mh, h1, eri, ne, norb = _mol("H2")
    lam_df = qubitization_lambda(h1, eri, norb, nelec=ne, shift=True)
    lam, lam_old = measurement_lambda(mh), measurement_lambda(mh, include_identity=True)

    assert advise(lam_old, lam_df, EPS, rho=1e4).cheaper == "FT"          # the old story
    assert advise(lam, lam_df, EPS, rho=1e4).cheaper == "near-term"       # the honest one

    for name in ("H2O(4,3)", "N2(6,6)"):
        mh, h1, eri, ne, norb = _mol(name)
        lam_df = qubitization_lambda(h1, eri, norb, nelec=ne, shift=True)
        for rho in (1.0, 1e2, 1e4, 1e6):
            old = advise(measurement_lambda(mh, include_identity=True), lam_df, EPS, rho=rho).cheaper
            new = advise(measurement_lambda(mh), lam_df, EPS, rho=rho).cheaper
            assert old == new, (name, rho, old, new)


# --- G4: the structural findings of SPEC_precision_cost survive -----------------------------


def test_G4_precision_cost_findings_survive_rescoring():
    """Under the honest metric: raw lambda_DF still exceeds lambda_meas for N2 (stronger now);
    the shifted lambda_DF still beats lambda_meas for every molecule; and the margin still
    grows with size -- but it is 1.94/2.95/3.58, NOT the published 2.8/5.4/5.7 (those measured
    identity mass, not FT advantage; the precision_cost G3 bar moves 5.0 -> 3.0)."""
    ratios = {}
    for name in MOLS:
        mh, h1, eri, ne, norb = _mol(name)
        lam = measurement_lambda(mh)
        lam_shift = qubitization_lambda(h1, eri, norb, nelec=ne, shift=True)
        assert lam_shift < lam, name
        ratios[name] = lam / lam_shift

    mh, h1, eri, ne, norb = _mol("N2(6,6)")
    assert qubitization_lambda(h1, eri, norb) > measurement_lambda(mh)  # raw DF > honest meas

    assert ratios["H2"] < ratios["H2O(4,3)"] < ratios["N2(6,6)"], ratios
    assert 3.0 < ratios["N2(6,6)"] < 5.0, ratios["N2(6,6)"]  # the honest margin, not the inflated 5.7
