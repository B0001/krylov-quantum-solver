"""
Acceptance gates G1-G4 for specs/SPEC_odmd_spin.md (spin spectroscopy: the interlayer J).

Test-first: ``odmd_spin`` does not exist yet, so this file is RED until the spec is implemented.
Claim: kicking the dimer ground state with the staggered magnetization S1z - S2z (the operator
that cannot move charge) exposes exactly the state the polarization leaves dark -- the m=0
triplet at omega = J -- giving the interlayer exchange constants of the Nb3X8 family from the
paper's own cRPA parameters, pinned by the closed form sqrt(((U0-Us)/2)^2+4t^2) - (U0-Us)/2.
Falsifiable physics rider: the Heisenberg superexchange 4t^2/(U0-Us) fails progressively across
the family (46.5% for the iodide). Sz|psi0> is again an exact eigenstate, so G1 re-exercises the
SPEC_odmd_optical degenerate-reference fix on a second operator.

PySCF/qiskit, no block2; `make gates` runs it in its own process.
"""
import numpy as np

from hybrid_quantum_solver.model_hamiltonians import ModelIntegrals
from nb3x8_gaps import NB3X8_LT_BULK, dimer_cluster_integrals
from odmd_optical import absorption_lines, dimer_polarization
from odmd_spin import (
    dimer_exchange_analytic,
    dimer_exchange_heisenberg,
    dimer_staggered_moment,
    spin_excitation_lines,
)

ORDER = ("Nb3F8", "Nb3Cl8", "Nb3Br8", "Nb3I8")
_CACHE = {}


def _dimer(name):
    if name not in _CACHE:
        base = dimer_cluster_integrals(**NB3X8_LT_BULK[name])
        mh = ModelIntegrals(base.h1, base.eri, 0.0, (1, 1), 2).to_hamiltonian()
        w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
        psi_hf = np.asarray(mh.hf_state().data, dtype=complex)
        pops = np.abs(V.conj().T @ psi_hf) ** 2
        psi0 = V[:, int(np.flatnonzero(pops > 1e-8)[0])]
        _CACHE[name] = (mh, psi0)
    return _CACHE[name]


def test_G1_one_line_at_exactly_J():
    """Each dimer's Sz-kicked spectrum: exactly ONE line at the analytic J (< 1e-9 relative),
    weight ||Sz|psi0>||^2 -- terminating instantly via the eigenstate short-circuit."""
    Sz = dimer_staggered_moment()
    for name in ORDER:
        mh, psi0 = _dimer(name)
        om, wt = spin_excitation_lines(mh, reference=psi0)
        assert len(om) == 1, (name, om)
        J = dimer_exchange_analytic(**NB3X8_LT_BULK[name])
        assert abs(om[0] - J) < 1e-9 * max(J, 1.0), (name, om[0], J)
        nrm2 = float(np.real((Sz @ psi0).conj() @ (Sz @ psi0)))
        assert abs(wt[0] - nrm2) < 1e-9, (name, wt[0], nrm2)


def test_G2_channel_complementarity():
    """The spin and optical lines never coincide (> 100 meV apart) and the spin line lies BELOW
    the optical line on every material -- magnetism is the low-energy physics."""
    for name in ORDER:
        mh, psi0 = _dimer(name)
        om_s, _ = spin_excitation_lines(mh, reference=psi0)
        om_p, _ = absorption_lines(mh, dimer_polarization(), reference=psi0)
        dist = min(abs(a - b) for a in om_s for b in om_p)
        assert dist > 100.0, (name, dist)
        assert om_s[0] < min(om_p), (name, om_s[0], min(om_p))


def test_G3_J_table_and_heisenberg_failure():
    """DEFINITION OF DONE: J strictly increasing F -> I; the Heisenberg 4t^2/(U0-Us) estimate
    errs < 1% for F (perturbative anchor) and > 30% for I (measured 46.5%), with the relative
    error strictly increasing across the family -- the iodide is beyond the Heisenberg regime."""
    js, errs = [], []
    for name in ORDER:
        p = NB3X8_LT_BULK[name]
        J = dimer_exchange_analytic(**p)
        js.append(J)
        errs.append(abs(dimer_exchange_heisenberg(**p) / J - 1.0))
    assert all(a < b for a, b in zip(js, js[1:])), js
    assert errs[0] < 0.01, errs[0]
    assert errs[-1] > 0.30, errs[-1]
    assert all(a < b for a, b in zip(errs, errs[1:])), errs


def test_G4_local_moment_ladder():
    """||Sz|psi0>||^2 = 1 to < 1e-3 for F (pure-spin limit), strictly decreasing F -> I
    (measured 1.000, 0.944, 0.890, 0.759): charge fluctuations eat 24% of the iodide's moment."""
    Sz = dimer_staggered_moment()
    weights = []
    for name in ORDER:
        _, psi0 = _dimer(name)
        weights.append(float(np.real((Sz @ psi0).conj() @ (Sz @ psi0))))
    assert abs(weights[0] - 1.0) < 1e-3, weights[0]
    assert all(a > b for a, b in zip(weights, weights[1:])), weights
    assert weights[-1] < 0.8, weights[-1]
