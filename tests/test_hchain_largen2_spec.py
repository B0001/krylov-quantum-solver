"""
Acceptance gates G1-G3 for specs/SPEC_hchain_largen2.md
(Hn to larger n, done right: adequate-D ramp + bulk per-site estimator).

G1 is pure arithmetic (no DMRG) and always runs. G2-G3 require block2 and MUST run in their own
process (block2 segfaults if it loads after pyscf/aer in the same interpreter -- `make gates`
isolates each spec file). `bulk_per_site_energy` does not exist yet: this file is RED until
Station 4 adds it to hybrid_quantum_solver.dmrg_reference (test-first).
"""
import pytest
from pyscf import gto, scf, ao2mo

from hybrid_quantum_solver.dmrg_reference import (
    DISCARD_WEIGHT_FLOOR,
    dmrg_energy_extrapolated,
    thermodynamic_limit_fit,
    dmrg_available,
)

requires_dmrg = pytest.mark.skipif(not dmrg_available(), reason="block2 not installed")

R_ANG = 1.8 * 0.529177210903   # 1.8 Bohr in Angstrom (~0.9525), near the cohesive minimum
CONV_DIMS = (80, 160, 300)     # well-converged on these small, FCI-tractable chains


def integrals(n):
    """Full-valence (n,n) minimal-basis integrals for an H_n chain."""
    atom = "; ".join(f"H 0 0 {i * R_ANG:.6f}" for i in range(n))
    mol = gto.M(atom=atom, basis="sto-6g", verbose=0)
    mf = scf.RHF(mol).run()
    mo = mf.mo_coeff
    h1 = mo.T @ mf.get_hcore() @ mo
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mo), n)
    ne = (mol.nelectron // 2, mol.nelectron // 2)
    return h1, eri, ne, float(mol.energy_nuc())


def test_G1_bulk_estimator_inverts_surface_term():
    """Synthetic open-chain totals E(n) = n*e_inf + c. The difference quotient must cancel the
    constant surface term c exactly, returning e_inf -- and must match the 1/n fit on ideal data."""
    from hybrid_quantum_solver.dmrg_reference import bulk_per_site_energy  # NEW (Station 4)

    e_inf, c = -0.539967, 0.123
    ns = [8, 10, 12, 16]
    totals = [n * e_inf + c for n in ns]

    e_bulk = bulk_per_site_energy(ns, totals)
    assert abs(e_bulk - e_inf) < 1e-9, e_bulk

    e_fit, _ = thermodynamic_limit_fit(ns, [t / n for t, n in zip(totals, ns)])
    assert abs(e_bulk - e_fit) < 1e-9, (e_bulk, e_fit)


@requires_dmrg
def test_G2_bulk_vs_fit_agree_on_real_chains():
    """The two independent extrapolations -- 1/n surface fit and bulk difference quotient -- must
    be mutually consistent on real H_n at converged bond dimension. At n <= 14 they agree only to
    ~0.5 mHa/atom (the bulk quotient is still on the curved part of E(n)); 0.1 mHa/atom agreement is
    a large-n claim verified by the driver. So this gate is a 1 mHa/atom sanity bound -- see
    SPEC_hchain_largen2 G2 finding + Definition of done."""
    from hybrid_quantum_solver.dmrg_reference import bulk_per_site_energy  # NEW (Station 4)

    ns = [8, 10, 12, 14]
    totals, per_atom = [], []
    for n in ns:
        h1, eri, ne, ec = integrals(n)
        res = dmrg_energy_extrapolated(h1, eri, ne, ec, bond_dims=CONV_DIMS)
        totals.append(res.energy)
        per_atom.append(res.energy / n)

    e_fit, _ = thermodynamic_limit_fit(ns, per_atom)
    e_bulk = bulk_per_site_energy(ns, totals)
    assert abs(e_bulk - e_fit) < 1e-3, (e_bulk, e_fit)   # 1 mHa/atom sanity bound (see G2 finding)


@requires_dmrg
def test_G3_stays_in_discarded_weight_regime():
    """Regime guard: at adequate bond dims the truncation must stay CONTROLLED.

    Was `method == "dweight"`, but that also rejected a converged ladder -- and `method="invD"`
    covers both "converged" and "uncontrolled", so it could not express the intent. The failure that
    killed the cheap spec was the *uncontrolled* one. See specs/SPEC_extrap_regime.md.

    This gate has no independent accuracy assertion, so a bare `!= "uncontrolled"` would accept a
    "converged" verdict on the label alone -- e.g. if block2 ever silently returned near-zero
    weights. Corroborate it against the weights themselves.
    """
    h1, eri, ne, ec = integrals(12)
    res = dmrg_energy_extrapolated(h1, eri, ne, ec, bond_dims=CONV_DIMS)
    assert res.regime != "uncontrolled", res.regime
    dws = [w for _, w, _ in res.per_D]
    assert res.regime != "converged" or max(dws) <= DISCARD_WEIGHT_FLOOR, dws
