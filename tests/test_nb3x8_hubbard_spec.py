"""
Acceptance gates G1-G5 for specs/SPEC_nb3x8_hubbard.md -- the Nb3X8 / Hubbard model loader.

Maps a tight-binding hopping matrix + Hubbard/cRPA interaction into the universal
``(h1, eri, e_core, nelec, norb)`` tuple, then checks it against the *analytic* two-site Hubbard
dimer energy and against PySCF FCI in the correct particle-number sector. The handoff
(CLAUDE_CODE_HANDOFF.md) names the Nb3X8 Model Database as the cleanest validation target; the bulk
low-T manifold reduces to a Hubbard dimer with an exact ground-state energy.

Uses only pyscf/qiskit (no block2), but ``make gates`` runs every test_*_spec.py in its own process.
"""
import numpy as np

from hybrid_quantum_solver.model_hamiltonians import (
    NB3X8_BULK_DIMER_PARAMS,
    fixed_filling_energy,
    hubbard_dimer_energy,
    hubbard_dimer_gap,
    hubbard_integrals,
    load_from_nb3x8_database,
)
from hybrid_quantum_solver.molecular_hamiltonian import build_hamiltonian_from_integrals
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver

TIGHT = 1e-6   # Ha: solver-vs-analytic
EXACT = 1e-9   # Ha: mapping-vs-FCI (both exact, only rounding)


def _dimer(t, U):
    """Half-filled (2-electron) two-site Hubbard dimer as a ModelIntegrals."""
    return hubbard_integrals(np.array([[0.0, -t], [-t, 0.0]]), U, nelec=2)


def test_G1_analytic_dimer_all_Uovert():
    """Number-conserving Krylov matches (U - sqrt(U^2 + 16 t^2))/2 from weak to deep-Mott U/t."""
    t = 1.0
    for U in (0.0, 2.0, 4.0, 8.0, 20.0):
        e_kry = QuantumKrylovSolver(_dimer(t, U).to_hamiltonian()).solve(8).energy
        assert abs(e_kry - hubbard_dimer_energy(t, U)) < TIGHT, (U, e_kry)
    # named limits: U=0 -> -2t ; large U -> Heisenberg -4t^2/U
    assert abs(hubbard_dimer_energy(1.0, 0.0) - (-2.0)) < EXACT
    assert abs(hubbard_dimer_energy(1.0, 1e6) - (-4.0 / 1e6)) < 1e-9


def test_G2_variational_floor():
    """A number-conserving Rayleigh-Ritz estimate cannot dip below the in-sector ground state."""
    t = 1.0
    for U in (0.0, 2.0, 4.0, 8.0, 20.0):
        e_kry = QuantumKrylovSolver(_dimer(t, U).to_hamiltonian()).solve(8).energy
        assert e_kry >= hubbard_dimer_energy(t, U) - 1e-9, (U, e_kry)


def test_G3_rank4_crpa_mapping_vs_fci():
    """A full rank-4 interaction tensor + generic hopping reproduces PySCF FCI (faithful mapping)."""
    norb = 3
    # generic symmetric hopping (nearest neighbour + a small next-nearest term)
    hop = np.array([[0.1, -1.0, -0.3],
                    [-1.0, 0.0, -1.0],
                    [-0.3, -1.0, -0.2]])
    # density-density Coulomb tensor: on-site U_i on the diagonal + inter-site U_ij (chemist (ii|jj))
    U_site = np.array([3.0, 2.5, 4.0])
    U_inter = np.array([[0.0, 1.0, 0.4],
                        [1.0, 0.0, 0.8],
                        [0.4, 0.8, 0.0]])
    coulomb = np.zeros((norb,) * 4)
    for i in range(norb):
        coulomb[i, i, i, i] = U_site[i]
        for j in range(norb):
            if i != j:
                coulomb[i, i, j, j] = U_inter[i, j]   # (ii|jj) density-density
    model = hubbard_integrals(hop, coulomb, nelec=(2, 1))

    e_solver = QuantumKrylovSolver(model.to_hamiltonian()).solve(16).energy
    e_fci_direct = fixed_filling_energy(model)
    # independent FCI straight from the raw tuple (no ModelIntegrals helper in the loop)
    from hybrid_quantum_solver.dmrg_reference import fci_energy
    e_fci_tuple = fci_energy(hop, coulomb, (2, 1), 0.0)

    assert abs(e_fci_direct - e_fci_tuple) < EXACT
    assert abs(e_solver - e_fci_direct) < TIGHT, (e_solver, e_fci_direct)
    assert e_solver >= e_fci_direct - 1e-9


def test_G4_fixed_filling_pitfall():
    """FINDING: full-Fock-space ground_state_energy() gives the WRONG filling at large U/t.

    For a Hubbard dimer with U/t >> 1 the global Fock-space minimum is the one-electron bonding
    state (E = -t), not the half-filled singlet. The number-conserving path (FCI in sector / Krylov)
    is required -- the concrete content of SKQD checklist item 5 (U(1)/electron-number conservation).
    """
    t, U = 1.0, 8.0
    model = _dimer(t, U)
    mh = build_hamiltonian_from_integrals(model.h1, model.eri, model.nelec, model.e_core)

    e_fullspace = mh.ground_state_energy()        # diagonalizes ALL particle numbers
    e_insector = fixed_filling_energy(model)      # stays at half filling
    e_analytic = hubbard_dimer_energy(t, U)

    # the in-sector / analytic answer is the half-filled singlet ...
    assert abs(e_insector - e_analytic) < EXACT, (e_insector, e_analytic)
    # ... and the full-space minimum is strictly lower: it has leaked to the 1-electron state (-t).
    assert e_fullspace < e_insector - 1e-3, (e_fullspace, e_insector)
    assert abs(e_fullspace - (-t)) < TIGHT, e_fullspace


def test_G5_units_and_nb3i8_anchor():
    """meV inputs round-trip through eV->Ha scaling; Nb3I8 bulk-LT dimer gives a singlet ground."""
    p = NB3X8_BULK_DIMER_PARAMS["Nb3I8"]
    t, U = p["t"], p["U"]

    model = load_from_nb3x8_database("Nb3I8")      # meV inputs, scaled to Ha internally
    e_solver = QuantumKrylovSolver(model.to_hamiltonian()).solve(8).energy

    from hybrid_quantum_solver.model_hamiltonians import HARTREE_PER_EV
    e_analytic_ha = hubbard_dimer_energy(t, U) * HARTREE_PER_EV / 1000.0   # meV -> Ha
    assert abs(e_solver - e_analytic_ha) < EXACT, (e_solver, e_analytic_ha)

    # singlet below triplet: positive gap, ~194 meV from the published params
    gap = hubbard_dimer_gap(t, U)
    assert gap > 0.0
    assert abs(gap - 194.1) < 1.0, gap

    # database mode: parsed hopping + scalar Coulomb reproduce the same model
    model_db = load_from_nb3x8_database(
        hopping=np.array([[0.0, -t], [-t, 0.0]]), coulomb=U, nelec=(1, 1), units="meV"
    )
    assert np.allclose(model_db.h1, model.h1)
    assert np.allclose(model_db.eri, model.eri)
