"""Gates for specs/SPEC_scf_stability_rule.md (chem-8ls).

Square H4's D2h-symmetric RHF is a saddle point at every side. A free SCF from PySCF's symmetric
minao guess keeps the symmetry up to rounding, so whether it stopped on the saddle or escaped to the
stable broken-symmetry minimum depended on the platform: at the default conv_tol it stopped on the
saddle at a = 1.05/1.10/1.19/1.22/1.35 on macOS arm64 and at 1.15/1.19/1.20/1.25/1.30 on Linux
x86_64. molecular_hamiltonian._driver_problem now follows internal SCF instabilities to a stable
minimum, so every builder lands on the same branch on both.
"""
import numpy as np
import pytest
from pyscf import gto, scf

import hybrid_quantum_solver.molecular_hamiltonian as mham
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from reachability import hf_population_spectrum


def _sq(a):                          # atom ORDER is load-bearing (SPEC_reachability_tolerance R2)
    return f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0"


# side -> (E_HF total in Ha, p0 = lowest-level HF population) of the builder's reference at
# conv_tol=1e-13, recorded with the rule on macOS 27 arm64 and on Linux x86_64 (identical to 1e-10 Ha
# and 8 digits; SPEC section 9).
MINIMUM = {
    1.00: (-1.7610750541, 0.47587995), 1.05: (-1.7751044780, 0.47270840),
    1.10: (-1.7825511826, 0.46922809), 1.15: (-1.7846392221, 0.46541280),
    1.19: (-1.7831249107, 0.46210250), 1.20: (-1.7823613168, 0.46123748),
    1.22: (-1.7804148460, 0.45946120), 1.25: (-1.7765301318, 0.45667855),
    1.30: (-1.7678162175, 0.45171429), 1.35: (-1.7567766045, 0.44632542),
    1.40: (-1.7438767269, 0.44049595),
}
SIDES = tuple(MINIMUM)


def _saddle(a):
    """The symmetry-enforced (D2h) RHF: PySCF's irrep-blocked SCF cannot leave the symmetric
    stationary point, so this is the saddle on every platform."""
    ms = scf.RHF(gto.M(atom=_sq(a), basis="sto3g", symmetry=True, verbose=0))
    ms.conv_tol = 1e-13
    ms.kernel()
    return ms


def _free_rhf(a, conv_tol=1e-13):
    mf = scf.RHF(gto.M(atom=_sq(a), basis="sto3g", verbose=0))
    mf.conv_tol = conv_tol
    return mf


def _lowest_hessian_eigenvalue(mf, mo, occ):
    """Lowest eigenvalue of the full real RHF orbital Hessian at (mo, occ), symmetry-breaking
    rotations included -- dense, so independent of the Davidson solver the rule itself uses."""
    from pyscf.soscf import newton_ah

    g, hop, _ = newton_ah.gen_g_hop_rhf(mf, mo, occ, with_symmetry=False)
    h = np.array([hop(e).real * 2 for e in np.eye(g.size)])
    return float(np.linalg.eigvalsh(0.5 * (h + h.T))[0])


# --- G1: one branch on every platform (DEFINITION OF DONE) ------------------------------------------

@pytest.mark.parametrize("conv_tol", (1e-6, 1e-9, 1e-13))
@pytest.mark.parametrize("a", SIDES)
def test_G1_builder_lands_on_the_recorded_minimum(a, conv_tol):
    """Killed by any side, conv_tol or platform whose build is not the recorded minimum. Before the
    rule this failed at 4-5 sides per platform, a different 4-5 on each."""
    e_ref, p0_ref = MINIMUM[a]
    mh = build_molecular_hamiltonian(atom=_sq(a), conv_tol=conv_tol)
    assert mh.hf_energy == pytest.approx(e_ref, abs=max(1e-8, 10 * conv_tol)), (a, conv_tol)
    assert float(hf_population_spectrum(mh)[0]) == pytest.approx(p0_ref, abs=1e-4), (a, conv_tol)


# --- G2: the reason for the rule -- the symmetric branch is a saddle, the rule's is a minimum -------

@pytest.mark.parametrize("a", SIDES)
def test_G2_symmetric_rhf_is_a_saddle_and_the_reference_a_minimum(a):
    """Killed by a side where the symmetric RHF is itself a minimum: the rule would then be choosing
    between two minima, a physics call it cannot make. Measured: -0.529..-0.672 vs +0.529..+0.661."""
    ms = _saddle(a)
    # STO-3G hydrogen has only s functions, so the symmetric solve's AO coefficients mean the same
    # thing in the unrotated frame its symmetry=True build re-oriented away from.
    assert _lowest_hessian_eigenvalue(_free_rhf(a), ms.mo_coeff, ms.mo_occ) < -0.5, a
    calc = mham._driver_problem(_sq(a), "sto3g", 0, 0, 1e-13, False)[2]._calc
    assert _lowest_hessian_eigenvalue(calc, calc.mo_coeff, calc.mo_occ) > 0.5, a
    assert 0.06 < ms.e_tot - calc.e_tot < 0.09, (a, ms.e_tot - calc.e_tot)


# --- G3: what the rule does, and where it does nothing ----------------------------------------------

@pytest.mark.parametrize("a", SIDES)
def test_G3_from_the_saddle_the_rule_reaches_the_minimum(a):
    """Platform-independent non-vacuity: a free SCF started ON the saddle stays there by itself, and
    the rule takes it to the recorded minimum."""
    ms = _saddle(a)
    mf = _free_rhf(a)
    mf.kernel(dm0=ms.make_rdm1())
    assert mf.e_tot == pytest.approx(ms.e_tot, abs=1e-10), a      # stuck on the saddle
    mham._follow_instabilities(mf)
    assert mf.e_tot == pytest.approx(MINIMUM[a][0], abs=1e-8), a


ORDINARY = {
    "H2": dict(atom="H 0 0 0; H 0 0 0.74"),
    "H2-stretched": dict(atom="H 0 0 0; H 0 0 2.0"),
    "H4-linear": dict(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"),
    "H4-triplet": dict(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0", spin=2),
    "HeH+": dict(atom="He 0 0 0; H 0 0 0.772", charge=1),
    "LiH": dict(atom="Li 0 0 0; H 0 0 1.6"),
    "LiH-CAS(2,5)": dict(atom="Li 0 0 0; H 0 0 1.6", active_electrons=2, active_orbitals=5),
    "N2-CAS(6,6)": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
    "OH-doublet": dict(atom="O 0 0 0; H 0 0 0.97", spin=1),
}


@pytest.mark.parametrize("name", ORDINARY)
def test_G3_stable_references_are_untouched(name, monkeypatch):
    """Where the free SCF is already stable the rule is a read-only check: the operator is
    bit-identical to a build with the rule switched off (the pre-chem-8ls builder)."""
    ruled = build_molecular_hamiltonian(**ORDINARY[name])
    monkeypatch.setattr(mham, "_follow_instabilities", lambda mf: None)
    free = build_molecular_hamiltonian(**ORDINARY[name])
    assert (ruled.qubit_hamiltonian.paulis.to_labels()
            == free.qubit_hamiltonian.paulis.to_labels()), name
    assert np.array_equal(ruled.qubit_hamiltonian.coeffs, free.qubit_hamiltonian.coeffs), name
    assert ruled.hf_energy == free.hf_energy, name


# --- G4: a persistent instability is an error, not a silent saddle ----------------------------------

def test_G4_an_instability_that_survives_the_steps_raises(monkeypatch):
    from pyscf.scf import stability

    mf = _free_rhf(1.0)
    mf.kernel()
    reruns = []
    monkeypatch.setattr(stability, "rhf_internal", lambda mf, **kw: (mf.mo_coeff, False))
    mf.kernel = lambda **kw: reruns.append(kw)
    with pytest.raises(RuntimeError, match="instability"):
        mham._follow_instabilities(mf, max_steps=3)
    assert len(reruns) == 3
