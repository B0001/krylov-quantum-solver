"""
Acceptance gates G1-G4 for specs/SPEC_odmd_optical.md (optical absorption + exciton binding).

Test-first: ``odmd_optical`` does not exist yet, so this file is RED until the spec is
implemented -- and G2 additionally HANGS against the unfixed ``odmd_spectral.reference_signal``
(probing exposed that an operator kick which lands on an exact eigenstate -- what P|psi0> IS on
the inversion-symmetric dimer: P odd, psi0 even, exactly one odd singlet -- gives reachable
width 0, tau = pi*1e12, and ~1e17 expm_multiply substeps; two probe processes burned 35 CPU-min
before being killed).

Claims: DMD poles/weights of an operator-kicked ground state are the bright FCI excitations and
|<E_n|O|psi0>|^2 (HeH+ cross-pins SPEC_qksd_properties' bright |mu|~0.85 -> 0.7224); the Nb3X8
dimers' only bright state is the odd singlet at exactly U0 (analytic), so their optical gaps and
exciton bindings Delta_c - Delta_opt are pinned twice (dense ED + closed form) -- and the
binding collapses from ~Us (F, atomic limit) to 0.26 Us (I): the exciton unbinds with hopping.
PySCF/qiskit, no block2; `make gates` runs it in its own process.
"""
import numpy as np

from hybrid_quantum_solver.model_hamiltonians import ModelIntegrals
from hybrid_quantum_solver.molecular_hamiltonian import (
    build_dipole_operators,
    build_molecular_hamiltonian,
)
from nb3x8_gaps import NB3X8_LT_BULK, dimer_cluster_integrals
from odmd_optical import (
    absorption_lines,
    dimer_exciton_binding,
    dimer_optical_gap,
    dimer_polarization,
)

ORDER = ("Nb3F8", "Nb3Cl8", "Nb3Br8", "Nb3I8")
_CACHE = {}


def _ground(mh):
    w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    psi_hf = np.asarray(mh.hf_state().data, dtype=complex)
    pops = np.abs(V.conj().T @ psi_hf) ** 2
    idx = int(np.flatnonzero(pops > 1e-8)[0])
    return w, V, V[:, idx], float(w[idx])


def _dimer(name):
    if name not in _CACHE:
        base = dimer_cluster_integrals(**NB3X8_LT_BULK[name])
        mh = ModelIntegrals(base.h1, base.eri, 0.0, (1, 1), 2).to_hamiltonian()
        _CACHE[name] = (mh, *_ground(mh))
    return _CACHE[name]


def test_G1_absorption_lines_exact_on_hehplus():
    """Poles and weights of the mu_z-kicked exact ground state match dense FCI < 1e-10; the
    bright transition weight cross-pins SPEC_qksd_properties (0.85^2); elastic line = <mu>^2."""
    spec = dict(atom="He 0 0 0; H 0 0 0.772", charge=1)
    mh = build_molecular_hamiltonian(**spec)
    w, V, psi0, e0 = _ground(mh)
    mu_z = build_dipole_operators(**spec)[2].to_matrix(sparse=True).tocsc()
    om, wt = absorption_lines(mh, mu_z, reference=psi0)
    exact_w = np.abs(V.conj().T @ (mu_z @ psi0)) ** 2
    assert len(om) >= 3
    for o, x in zip(om, wt):
        j = int(np.argmin(np.abs(w - (e0 + o))))
        assert abs(w[j] - (e0 + o)) < 1e-10, (o, w[j] - (e0 + o))
        assert abs(exact_w[j] - x) < 1e-10, (o, exact_w[j] - x)
    elastic = wt[int(np.argmin(np.abs(om)))]
    perm = float(np.real(psi0.conj() @ (mu_z @ psi0)))
    assert abs(elastic - perm**2) < 1e-10
    bright = float(sorted(wt)[-2])
    assert abs(bright - 0.7224) < 1e-3, bright             # |mu|~0.85 from SPEC_qksd_properties


def test_G2_eigenstate_kick_terminates_and_is_exact():
    """THE FIX'S GATE: P|psi0> on the Nb3I8 dimer is an exact eigenstate -- pre-fix this hangs
    (~1e17 expm substeps). Fixed: returns with exactly ONE line at pole = U0 (< 1e-9 relative),
    weight = ||P|psi0>||^2."""
    mh, w, V, psi0, e0 = _dimer("Nb3I8")
    P = dimer_polarization()
    om, wt = absorption_lines(mh, P, reference=psi0)
    assert len(om) == 1, om
    U0 = NB3X8_LT_BULK["Nb3I8"]["U0"]
    assert abs((e0 + om[0]) - U0) < 1e-9 * U0, e0 + om[0]
    nrm2 = float(np.real((P @ psi0).conj() @ (P @ psi0)))
    assert abs(wt[0] - nrm2) < 1e-9, (wt[0], nrm2)


def test_G3_optical_gaps_and_exciton_binding():
    """DEFINITION OF DONE: ODMD optical gap == analytic odd-singlet formula (< 1e-6 meV) on all
    four materials; binding(F) within 2% of Us (the atomic limit); binding/Us strictly
    decreasing F -> Cl -> Br -> I -- the exciton unbinds with hopping."""
    ratios = []
    for name in ORDER:
        p = NB3X8_LT_BULK[name]
        mh, w, V, psi0, e0 = _dimer(name)
        om, wt = absorption_lines(mh, dimer_polarization(), reference=psi0)
        gap_odmd = float(min(o for o in om if o > 1e-6))
        assert abs(gap_odmd - dimer_optical_gap(**p)) < 1e-6, (name, gap_odmd)
        ratios.append(dimer_exciton_binding(**p) / p["Us"])
    assert abs(ratios[0] - 1.0) < 0.02, ratios[0]          # F: atomic limit, binding ~ Us
    assert all(a > b for a, b in zip(ratios, ratios[1:])), ratios


def test_G4_selection_rules_and_brightness_ladder():
    """Exactly one bright line per dimer while >= 4 levels exist (the odd-singlet selection
    rule); total oscillator weight ||P|psi0>||^2 strictly increases F -> I (the 4-orders
    polarizability ladder)."""
    weights = []
    for name in ORDER:
        mh, w, V, psi0, e0 = _dimer(name)
        P = dimer_polarization()
        om, wt = absorption_lines(mh, P, reference=psi0)
        assert len(om) == 1, (name, om)
        assert len(np.unique(np.round(w, 6))) >= 4, name
        weights.append(float(np.real((P @ psi0).conj() @ (P @ psi0))))
    assert all(a < b for a, b in zip(weights, weights[1:])), weights
    assert weights[0] < 1e-3 < weights[-1], weights        # near-dark F, bright I
