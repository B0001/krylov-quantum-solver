"""Gates for specs/SPEC_eigenstate_reachability.md (chem-obf).

SPEC_reachability_tolerance showed no fixed |<HF|psi_k>|^2 threshold separates a physical overlap
from SCF-convergence residue on a symmetry-forbidden level. SPEC_symmetry_reachability only decided
WHETHER a symmetry filter is possible per system. This gate file pins the per-eigenstate decision
(majority weight in HF's exact symmetry sector) and checks that the threshold sites consume it.

Dense eigendecompositions here use the REAL matrix (these Hamiltonians are real): dsyevd is fast
and not hit by the macOS Accelerate ZHEEVD bug, so the gate exercises the decision, not LAPACK.
"""
import dataclasses
from functools import lru_cache

import numpy as np
import pytest
from qiskit.quantum_info import SparsePauliOp

from hybrid_quantum_solver.molecular_hamiltonian import (
    build_hamiltonian_from_integrals,
    build_molecular_hamiltonian,
)
from reachability import (
    REACHABLE_TOL_CERTIFIED,
    hf_symmetry_sector,
    orbital_parity_bounds,
    reachable_eigenpairs,
    reachable_mask,
)

ODMD_TOL = 1e-8                     # the ODMD/MSD/Trotter family's population cut


def _sq(a):                          # atom ORDER is load-bearing (SPEC_reachability_tolerance R2)
    return f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0"


LIN_H4 = "H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"
H2 = "H 0 0 0; H 0 0 0.74"
# (side, conv_tol, the site tol at which the residue level clears the population cut)
WITNESSES = [(1.10, 1e-6, REACHABLE_TOL_CERTIFIED), (1.19, 1e-9, ODMD_TOL)]


@lru_cache(maxsize=None)
def _system(atom, conv_tol=1e-9, **kw):
    mh = build_molecular_hamiltonian(atom=atom, conv_tol=conv_tol, **kw)
    h = mh.qubit_hamiltonian.to_matrix()
    assert np.abs(h.imag).max() < 1e-12
    w, V = np.linalg.eigh(h.real)
    pops = (V.T @ np.asarray(mh.hf_state().data).real) ** 2
    return mh, w, V, pops


@lru_cache(maxsize=None)
def _fci(atom, wfnsym, symmetry="D2h"):
    """PySCF symmetry-adapted FCI, S_z = 0, tight symmetric SCF -- independent of the qubit path."""
    from pyscf import fci, gto, scf

    mol = gto.M(atom=atom, basis="sto3g", symmetry=symmetry, verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-13
    mf.kernel()
    solver = fci.FCI(mf)
    solver.wfnsym, solver.nroots = wfnsym, 1000
    e, _ = solver.kernel()
    return np.sort(np.atleast_1d(e)) - mol.energy_nuc()


def _weights(mh, V):
    return (V[hf_symmetry_sector(mh)] ** 2).sum(axis=0)


# --- G1: the witnesses (DEFINITION OF DONE) -------------------------------------------------------

@pytest.mark.parametrize("a, conv_tol, tol", WITNESSES)
def test_G1_residue_level_rejected_and_the_Ag_ground_state_selected(a, conv_tol, tol):
    mh, w, V, pops = _system(_sq(a), conv_tol)
    keep, cut = reachable_mask(mh, w, V, pops, tol), pops > tol
    assert (cut & ~keep).any(), "population cut admits nothing the decision rejects"
    assert w[cut].min() == pytest.approx(_fci(_sq(a), "B1g")[0], abs=1e-8)   # the forbidden state
    assert w[keep].min() == pytest.approx(_fci(_sq(a), "Ag")[0], abs=1e-8)   # HF's irrep


# --- G2: brute force against an independent symmetry-adapted FCI ----------------------------------

G2_CASES = [(_sq(a), ct) for a in (1.05, 1.10, 1.19, 1.35) for ct in (1e-6, 1e-9)]
G2_CASES += [(LIN_H4, 1e-9), (H2, 1e-9)]


@pytest.mark.parametrize("atom, conv_tol", G2_CASES)
def test_G2_allowed_spectrum_equals_symmetry_adapted_fci(atom, conv_tol):
    mh, w, V, _ = _system(atom, conv_tol)
    allowed = np.sort(w[_weights(mh, V) > 0.5])          # eigenvectors IN the sector
    ref = _fci(atom, "Ag")
    assert len(allowed) == len(ref), (len(allowed), len(ref))
    assert np.abs(allowed - ref).max() < 1e-8


# --- G3: the bead's prediction -- no recorded number moves on the ordinary gated systems ----------

CONTROLS = {
    "H2": dict(atom=H2),
    "H2-stretched": dict(atom="H 0 0 0; H 0 0 2.0"),
    "H4-linear": dict(atom=LIN_H4),
    "HeH+": dict(atom="He 0 0 0; H 0 0 0.772", charge=1),
    "LiH": dict(atom="Li 0 0 0; H 0 0 1.6"),
    "LiH-CAS(2,5)": dict(atom="Li 0 0 0; H 0 0 1.6", active_electrons=2, active_orbitals=5),
    "N2-CAS(6,6)": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
}


@pytest.mark.parametrize("name", CONTROLS)
def test_G3_decision_equals_population_cut_on_ordinary_systems(name):
    mh, w, V, pops = _system(**CONTROLS[name])
    for tol in (ODMD_TOL, REACHABLE_TOL_CERTIFIED):
        assert np.array_equal(reachable_mask(mh, w, V, pops, tol), pops > tol), (name, tol)


# --- G4: the bound is load-bearing ----------------------------------------------------------------

def _naive_weights(mh, V):
    """Sector weights with EVERY parity bit, ignoring the bound -- what the decision must not do."""
    n, no = mh.num_qubits, mh.num_spatial_orbitals
    x = np.arange(2 ** n)
    occ = (x[:, None] >> np.arange(n)) & 1
    hf = int(np.argmax(np.abs(np.asarray(mh.hf_state().data))))
    sec = (occ[:, :no].sum(1) == occ[hf, :no].sum()) & (occ[:, no:].sum(1) == occ[hf, no:].sum())
    for odd, _ in orbital_parity_bounds(mh):
        par = (occ[:, :no] @ odd + occ[:, no:] @ odd) % 2
        sec &= par == par[hf]
    return (V[sec] ** 2).sum(axis=0)


@pytest.mark.parametrize("a", (1.20, 1.40))
def test_G4_broken_rhf_bit_is_excluded_and_nothing_is_vetoed(a):
    mh, w, V, pops = _system(_sq(a))
    assert max(beta for _, beta in orbital_parity_bounds(mh)) >= 1.0
    for tol in (ODMD_TOL, REACHABLE_TOL_CERTIFIED):
        assert np.array_equal(reachable_mask(mh, w, V, pops, tol), pops > tol), tol


def test_G4_naive_labels_would_veto_a_physical_level():
    mh, w, V, pops = _system(_sq(1.20))
    vetoed = (pops > 1e-3) & (_naive_weights(mh, V) < 0.5)
    assert vetoed.any(), "naive labels veto nothing -- the bound would be untested"


def test_G4_lih_pi_rotation_bit_is_excluded():
    """The free SCF rotates LiH's degenerate pi pair: the sigma_v bit is unusable, C2 is exact."""
    mh = _system(**CONTROLS["LiH"])[0]
    betas = sorted(beta for _, beta in orbital_parity_bounds(mh))
    assert betas[0] < 1e-10 and betas[-1] >= 1.0, betas


# --- G5: the fallback is the old behaviour --------------------------------------------------------

def test_G5_integrals_path_uses_only_particle_sectors():
    from pyscf import ao2mo, gto, scf

    mol = gto.M(atom=H2, basis="sto3g", verbose=0)
    mf = scf.RHF(mol).run()
    c = mf.mo_coeff
    eri = ao2mo.restore(1, ao2mo.full(mol, c), c.shape[1])
    mh = build_hamiltonian_from_integrals(c.T @ mf.get_hcore() @ c, eri, (1, 1), mol.energy_nuc())
    assert mh.build_args is None and orbital_parity_bounds(mh) is None
    assert int(hf_symmetry_sector(mh).sum()) == 4                       # C(2,1)^2: (N_a, N_b) only
    h = mh.qubit_hamiltonian.to_matrix().real
    w, V = np.linalg.eigh(h)
    pops = (V.T @ np.asarray(mh.hf_state().data).real) ** 2
    for tol in (ODMD_TOL, REACHABLE_TOL_CERTIFIED):
        assert np.array_equal(reachable_mask(mh, w, V, pops, tol), pops > tol)


def test_G5_stale_build_args_disable_spatial_bits_but_a_shift_does_not():
    mh = _system(_sq(1.19))[0]
    stale = dataclasses.replace(mh, qubit_hamiltonian=_system(_sq(1.10))[0].qubit_hamiltonian)
    assert orbital_parity_bounds(stale) is None
    shifted = dataclasses.replace(
        mh, qubit_hamiltonian=(mh.qubit_hamiltonian + SparsePauliOp("I" * 8, 0.37)).simplify())
    assert orbital_parity_bounds(shifted) is not None


# --- G6: the majority decision is never near its boundary -----------------------------------------

G6_SYSTEMS = [dict(atom=_sq(1.10), conv_tol=1e-6), dict(atom=_sq(1.19)), dict(atom=_sq(1.20)),
              dict(atom=LIN_H4), CONTROLS["LiH"], CONTROLS["N2-CAS(6,6)"]]


@pytest.mark.parametrize("kw", G6_SYSTEMS, ids=lambda kw: kw["atom"][:16])
def test_G6_populated_eigenvectors_are_sector_pure(kw):
    mh, w, V, pops = _system(**kw)
    wt = _weights(mh, V)[pops > REACHABLE_TOL_CERTIFIED]
    assert np.minimum(wt, 1 - wt).max() < 1e-3


# --- G7: the sites consume the decision -----------------------------------------------------------

@pytest.mark.parametrize("a, conv_tol, tol", WITNESSES)
def test_G7_threshold_sites_use_the_decision(a, conv_tol, tol):
    from device_odmd import centered_frame
    from hf_overlap_certificate import exact_reachable_overlap
    from hf_overlap_subspace import exact_hf_subspace_overlap
    from msd import build_msd_problem
    from odmd import build_odmd_problem
    from trotter_odmd import build_trotter_odmd_problem
    from trotter_resolution_floor import _centered

    mh, w, V, pops = _system(_sq(a), conv_tol)
    keep = reachable_mask(mh, w, V, pops, ODMD_TOL)
    reach = w[keep]
    mu, tau = 0.5 * (reach.max() + reach.min()), np.pi / (reach.max() - reach.min())
    old = w[pops > ODMD_TOL]
    assert abs(0.5 * (old.max() + old.min()) - mu) > 1e-3     # the witness really moves the frame

    _, tau_d, mu_d = centered_frame(mh)
    frames = {"device_odmd": (mu_d, tau_d), "trotter_resolution_floor": _centered(mh)[2:4]}
    for name, p in {"odmd": build_odmd_problem(mh, n=2), "msd": build_msd_problem(mh, n=2),
                    "trotter_odmd": build_trotter_odmd_problem(mh, n=2)}.items():
        frames[name] = (p.mu, p.tau)
    for name, (m, t) in frames.items():
        assert m == pytest.approx(mu, abs=1e-10) and t == pytest.approx(tau, rel=1e-10), name

    e_cert, v_cert = reachable_eigenpairs(mh)              # certified arc, tol 1e-10
    lowest = w[reachable_mask(mh, w, V, pops, REACHABLE_TOL_CERTIFIED)].min()
    assert e_cert[0] == pytest.approx(lowest, abs=1e-10)
    ov = float(np.sqrt(pops[np.flatnonzero(keep)[0]]))     # HF overlap of the Ag ground state
    assert exact_reachable_overlap(mh) == pytest.approx(ov, rel=1e-8)
    assert exact_hf_subspace_overlap(mh, 1) == pytest.approx(ov, rel=1e-8)
    assert ov > 0.5                                         # a physical overlap, not residue


# --- G8: open shells -- exactly degenerate M_s pairs (found in code review) -----------------------

OPEN_SHELL = {"H4-triplet": dict(atom=LIN_H4, spin=2),
              "OH-doublet": dict(atom="O 0 0 0; H 0 0 0.97", spin=1)}


@pytest.mark.parametrize("name", OPEN_SHELL)
def test_G8_open_shell_degenerate_pairs_are_decided_per_eigenspace(name):
    """eigh may mix the M_s = +S/-S copies (different (N_a, N_b) sectors, same energy). Deciding per
    cluster of degenerate eigenvalues keeps every populated level; the per-vector majority test
    (this module's first version) did not -- H4 triplet: 2 levels with p > 1e-3 dropped."""
    mh, w, V, pops = _system(**OPEN_SHELL[name])
    for tol in (ODMD_TOL, REACHABLE_TOL_CERTIFIED):
        assert np.array_equal(reachable_mask(mh, w, V, pops, tol), pops > tol), (name, tol)
    if name == "H4-triplet":   # non-vacuous: in a basis that mixes the M_s copies the per-vector rule
        # drops a populated level. Whether LAPACK's eigh mixes them depends on the build and on the
        # operator's low bits (canonical term order unmixes them on macOS), so mix deliberately:
        # a fixed orthogonal rotation of each populated degenerate cluster -- same span, so still
        # an orthonormal eigenbasis of H.
        Vm, rng, seen = V.copy(), np.random.default_rng(0), set()
        for p in np.flatnonzero(pops > 1e-3):
            cl = np.flatnonzero(np.abs(w - w[p]) < 1e-10)
            if len(cl) > 1 and cl[0] not in seen:
                seen.add(cl[0])
                Vm[:, cl] = V[:, cl] @ np.linalg.qr(rng.standard_normal((len(cl), len(cl))))[0]
        pm = (Vm.T @ np.asarray(mh.hf_state().data).real) ** 2
        for tol in (ODMD_TOL, REACHABLE_TOL_CERTIFIED):   # the decision is basis-invariant...
            assert np.array_equal(reachable_mask(mh, w, Vm, pm, tol), pm > tol), (name, tol)
        assert ((pm > 1e-3) & (_weights(mh, Vm) < 0.5)).any()   # ...the per-vector rule is not
