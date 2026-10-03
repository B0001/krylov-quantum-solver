"""Gates for specs/SPEC_eigenstate_reachability.md (chem-obf), revised after merge.

SPEC_reachability_tolerance showed no fixed |<HF|psi_k>|^2 threshold separates a physical overlap
from SCF-convergence residue on a symmetry-forbidden level. PR #62 vetoed such levels from every
reference; G9 kills that: QKSD from |HF> converges to the very level the veto removed, so references
stay on the population cut and the per-eigenstate symmetry decision only WARNS. This file pins the
diagnostic against an independent symmetric FCI (G1, G2), its silence on ordinary systems (G3-G6,
G8), the sites' population-cut frames (G7), the falsification (G9) and that it never raises (G10).

Dense eigendecompositions here use the REAL matrix (these Hamiltonians are real): dsyevd is fast
and not hit by the macOS Accelerate ZHEEVD bug, so the gate exercises the decision, not LAPACK.
"""
import dataclasses
import warnings
from functools import lru_cache

import numpy as np
import pytest
from qiskit.quantum_info import SparsePauliOp

import reachability
from hybrid_quantum_solver.molecular_hamiltonian import (
    build_hamiltonian_from_integrals,
    build_molecular_hamiltonian,
)
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from reachability import (
    REACHABLE_TOL_CERTIFIED,
    TIGHT_SCF_CONV_TOL,
    hf_symmetry_sector,
    orbital_parity_bounds,
    reachable_eigenpairs,
    reachable_mask,
    symmetry_allowed,
)

ODMD_TOL = 1e-8                     # the ODMD/MSD/Trotter family's population cut
TOLS = (ODMD_TOL, REACHABLE_TOL_CERTIFIED)
DEFINED = 1e-20                     # HF projection above roundoff: where the decision is defined


def _sq(a):                          # atom ORDER is load-bearing (SPEC_reachability_tolerance R2)
    return f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0"


LIN_H4 = "H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"
H2 = "H 0 0 0; H 0 0 0.74"
# (side, conv_tol). Both residue levels clear both site cuts (macOS: 2.66e-6 and 5.48e-8); a = 1.19
# left the default conv_tol, where its 1.40e-8 cleared the 1e-8 cut by only 1.4x.
WITNESSES = [(1.10, 1e-6), (1.19, 1e-6)]


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


def _residue_warnings(mh, w, V, pops, tol):
    """(reachable_mask, the residue warnings it emitted)."""
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        mask = reachable_mask(mh, w, V, pops, tol)
    return mask, [r for r in rec if issubclass(r.category, RuntimeWarning)
                  and "SCF residue" in str(r.message)]


# --- G1: the witnesses -- identified, warned, kept ------------------------------------------------

@pytest.mark.parametrize("tol", TOLS)
@pytest.mark.parametrize("a, conv_tol", WITNESSES)
def test_G1_residue_level_identified_warned_and_kept(a, conv_tol, tol):
    mh, w, V, pops = _system(_sq(a), conv_tol)
    mask, rec = _residue_warnings(mh, w, V, pops, tol)
    cut = pops > tol
    assert np.array_equal(mask, cut) and len(rec) == 1, len(rec)
    allowed = symmetry_allowed(mh, w, V)
    assert w[cut & ~allowed].min() == pytest.approx(_fci(_sq(a), "B1g")[0], abs=1e-8)  # forbidden
    assert w[cut & allowed].min() == pytest.approx(_fci(_sq(a), "Ag")[0], abs=1e-8)    # HF's irrep


# --- G2: brute force against an independent symmetry-adapted FCI, every populated level -----------

G2_CASES = [(_sq(a), ct) for a in (1.05, 1.10, 1.19, 1.35) for ct in (1e-6, 1e-9)]
G2_CASES += [(LIN_H4, 1e-9), (H2, 1e-9)]


@pytest.mark.parametrize("atom, conv_tol", G2_CASES)
def test_G2_decision_equals_symmetry_adapted_fci(atom, conv_tol):
    mh, w, V, pops = _system(atom, conv_tol)
    ref = _fci(atom, "Ag")
    in_sector = np.sort(w[_weights(mh, V) > 0.5])          # the sector itself, every eigenvector
    assert len(in_sector) == len(ref), (len(in_sector), len(ref))
    assert np.abs(in_sector - ref).max() < 1e-8
    defined = pops > DEFINED
    in_ag = np.array([np.abs(ref - e).min() < 1e-8 for e in w[defined]])
    assert in_ag.sum() >= 2                                 # above the minimum, not just the ground
    assert np.array_equal(symmetry_allowed(mh, w, V)[defined], in_ag)


# --- G3: no false alarm on the ordinary gated systems ---------------------------------------------

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
def test_G3_no_false_alarm_on_ordinary_systems(name):
    mh, w, V, pops = _system(**CONTROLS[name])
    allowed = symmetry_allowed(mh, w, V)
    for tol in TOLS:
        mask, rec = _residue_warnings(mh, w, V, pops, tol)
        assert np.array_equal(mask, pops > tol) and not rec, (name, tol)
        assert allowed[pops > tol].all(), (name, tol)


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
def test_G4_broken_rhf_bit_is_excluded_and_nothing_is_flagged(a):
    mh, w, V, pops = _system(_sq(a))
    assert max(beta for _, beta in orbital_parity_bounds(mh)) >= 1.0
    assert symmetry_allowed(mh, w, V)[pops > REACHABLE_TOL_CERTIFIED].all()


def test_G4_naive_labels_would_flag_a_physical_level():
    mh, w, V, pops = _system(_sq(1.20))
    flagged = (pops > 1e-3) & (_naive_weights(mh, V) < 0.5)
    assert flagged.any(), "naive labels flag nothing -- the bound would be untested"


def test_G4_lih_c2_bit_is_exact():
    """PR #62 also asserted LiH's sigma_v bit has beta >= 1. Dropped: that beta measures the
    arbitrary rotation SCF/LAPACK leave inside LiH's degenerate pi pair, not LiH (3.2 here, 0.646 at
    conv_tol=1e-7, 0.606 with the atom order swapped -- SPEC section 5, G4)."""
    mh = _system(**CONTROLS["LiH"])[0]
    assert min(beta for _, beta in orbital_parity_bounds(mh)) < 1e-10


# --- G5: the fallback ----------------------------------------------------------------------------

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
    for tol in TOLS:
        mask, rec = _residue_warnings(mh, w, V, pops, tol)
        assert np.array_equal(mask, pops > tol) and not rec


def test_G5_stale_build_args_disable_spatial_bits_but_a_shift_does_not():
    mh = _system(_sq(1.19))[0]
    stale = dataclasses.replace(mh, qubit_hamiltonian=_system(_sq(1.10))[0].qubit_hamiltonian)
    assert orbital_parity_bounds(stale) is None
    shifted = dataclasses.replace(
        mh, qubit_hamiltonian=(mh.qubit_hamiltonian + SparsePauliOp("I" * 8, 0.37)).simplify())
    assert orbital_parity_bounds(shifted) is not None


# --- G6: the decision is the rounding of a quantity never near its boundary -----------------------

G6_SYSTEMS = [dict(atom=_sq(1.10), conv_tol=1e-6), dict(atom=_sq(1.19), conv_tol=1e-6),
              dict(atom=_sq(1.20)), dict(atom=LIN_H4), CONTROLS["LiH"], CONTROLS["N2-CAS(6,6)"]]


@pytest.mark.parametrize("kw", G6_SYSTEMS, ids=lambda kw: kw["atom"][:16])
def test_G6_populated_eigenvectors_are_sector_pure_and_decided_by_rounding(kw):
    mh, w, V, pops = _system(**kw)
    populated = pops > REACHABLE_TOL_CERTIFIED
    wt = _weights(mh, V)[populated]
    assert np.minimum(wt, 1 - wt).max() < 1e-3
    assert np.array_equal(symmetry_allowed(mh, w, V)[populated], wt > 0.5)


# --- G7: the sites use the population cut ---------------------------------------------------------

@pytest.mark.parametrize("a, conv_tol", WITNESSES)
def test_G7_threshold_sites_use_the_population_cut(a, conv_tol):
    from device_odmd import centered_frame
    from hf_overlap_certificate import exact_reachable_overlap
    from hf_overlap_subspace import exact_hf_subspace_overlap
    from msd import build_msd_problem
    from odmd import build_odmd_problem
    from trotter_odmd import build_trotter_odmd_problem
    from trotter_resolution_floor import _centered

    mh, w, V, pops = _system(_sq(a), conv_tol)
    reach = w[pops > ODMD_TOL]
    mu, tau = 0.5 * (reach.max() + reach.min()), np.pi / (reach.max() - reach.min())
    _, tau_d, mu_d = centered_frame(mh)
    frames = {"device_odmd": (mu_d, tau_d), "trotter_resolution_floor": _centered(mh)[2:4]}
    for name, p in {"odmd": build_odmd_problem(mh, n=2), "msd": build_msd_problem(mh, n=2),
                    "trotter_odmd": build_trotter_odmd_problem(mh, n=2)}.items():
        frames[name] = (p.mu, p.tau)
    for name, (m, t) in frames.items():
        assert m == pytest.approx(mu, abs=1e-10) and t == pytest.approx(tau, rel=1e-10), name

    e_b1g = _fci(_sq(a), "B1g")[0]
    assert reachable_eigenpairs(mh)[0][0] == pytest.approx(e_b1g, abs=1e-8)   # certified arc
    ov = float(np.sqrt(pops[np.argmin(np.abs(w - e_b1g))]))   # the residue level's HF overlap
    assert exact_reachable_overlap(mh) == pytest.approx(ov, rel=1e-8)
    assert exact_hf_subspace_overlap(mh, 1) == pytest.approx(ov, rel=1e-8)


# --- G8: open shells -- exactly degenerate M_s pairs (found in the first code review) -------------

OPEN_SHELL = {"H4-triplet": dict(atom=LIN_H4, spin=2),
              "OH-doublet": dict(atom="O 0 0 0; H 0 0 0.97", spin=1)}


@pytest.mark.parametrize("name", OPEN_SHELL)
def test_G8_open_shell_degenerate_pairs_are_decided_per_eigenspace(name):
    """eigh may mix the M_s = +S/-S copies (different (N_a, N_b) sectors, same energy). Deciding per
    cluster of degenerate eigenvalues allows every populated level; the per-vector majority test
    (this module's first version) did not -- H4 triplet: 2 levels with p > 1e-3 dropped."""
    mh, w, V, pops = _system(**OPEN_SHELL[name])
    for tol in TOLS:
        mask, rec = _residue_warnings(mh, w, V, pops, tol)
        assert np.array_equal(mask, pops > tol) and not rec, (name, tol)
    assert symmetry_allowed(mh, w, V)[pops > REACHABLE_TOL_CERTIFIED].all()
    if name == "H4-triplet":   # non-vacuous: in a basis that mixes the M_s copies the per-vector rule
        # flags a populated level. Whether LAPACK's eigh mixes them depends on the build and on the
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
        assert symmetry_allowed(mh, w, Vm)[pm > REACHABLE_TOL_CERTIFIED].all()  # basis-invariant...
        assert ((pm > 1e-3) & (_weights(mh, Vm) < 0.5)).any()   # ...the per-vector rule is not


# --- G9: the falsification -- QKSD from |HF> converges to the residue level -----------------------

G9_M = range(28, 33)                 # every M in [28, 32], pre-registered (SPEC section 5)


def _qksd_electronic(mh):
    steps = QuantumKrylovSolver(mh).convergence(max(G9_M))
    return np.array([s.energy for s in steps]) - mh.energy_offset       # index M - 1


@pytest.mark.parametrize("a, conv_tol", WITNESSES + [(1.19, None)])    # None: builder default
def test_G9_qksd_from_hf_converges_to_the_residue_level(a, conv_tol):
    mh = build_molecular_hamiltonian(atom=_sq(a), **({} if conv_tol is None
                                                     else {"conv_tol": conv_tol}))
    e_b1g, e_ag = _fci(_sq(a), "B1g")[0], _fci(_sq(a), "Ag")[0]
    e = _qksd_electronic(mh)[[m - 1 for m in G9_M]]
    assert (e - e_b1g).min() >= -1e-9, e - e_b1g            # the variational floor
    assert (e - e_b1g).max() < 1e-2, e - e_b1g              # ...and the solver sits on B1g
    assert (e < e_ag - 0.1).all()                           # far below the level the veto kept
    assert reachable_eigenpairs(mh)[0][0] == pytest.approx(e_b1g, abs=1e-8)   # the reference


@pytest.mark.parametrize("a", (1.10, 1.19))
def test_G9_control_tight_scf_stays_on_the_ag_level(a):
    mh = build_molecular_hamiltonian(atom=_sq(a), conv_tol=TIGHT_SCF_CONV_TOL)
    e_ag = _fci(_sq(a), "Ag")[0]
    e = _qksd_electronic(mh)                                # every M <= 32
    assert e.min() >= e_ag - 1e-8, e.min() - e_ag
    assert reachable_eigenpairs(mh)[0][0] == pytest.approx(e_ag, abs=1e-8)


# --- G10: the diagnostic never raises -------------------------------------------------------------

HALF_DETECTED = ["H 0 0 0; H 1.100015 0 0; H 1.1 1.1 0; H 0 1.1 0",   # IndexError before (PySCF)
                 "H 0 0 0; H 1.10002 0 0; H 1.1 1.1 0; H 0 1.1 0"]    # PointGroupSymmetryError


@pytest.mark.parametrize("atom", HALF_DETECTED)
def test_G10_half_detected_symmetry_does_not_raise(atom):
    mh, w, V, pops = _system(atom)
    for tol in TOLS:
        assert np.array_equal(reachable_mask(mh, w, V, pops, tol), pops > tol)


def test_G10_any_diagnostic_failure_means_no_diagnostic(monkeypatch):
    mh, w, V, pops = _system(_sq(1.10), 1e-6)              # a witness: the diagnostic would warn

    def boom(*_):
        raise RuntimeError("diagnostic failure")

    monkeypatch.setattr(reachability, "symmetry_allowed", boom)
    mask, rec = _residue_warnings(mh, w, V, pops, ODMD_TOL)
    assert np.array_equal(mask, pops > ODMD_TOL) and not rec


def test_G10_odd_electron_active_space_is_verified():
    mh, w, V, pops = _system("O 0 0 0; H 0 0 0.97", spin=1, active_electrons=(3, 2),
                             active_orbitals=4)
    assert orbital_parity_bounds(mh) is not None             # the rebuild reproduced the operator
    assert symmetry_allowed(mh, w, V)[pops > REACHABLE_TOL_CERTIFIED].all()
