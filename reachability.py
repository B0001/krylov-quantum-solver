#!/usr/bin/env python3
"""
reachability.py -- deciding which eigenstates a Hartree-Fock reference can actually reach.

The certified arc and the ODMD/MSD family both define the "HF-reachable" sector by thresholding
the HF population, ``|<HF|psi_k>|^2 > tol``. specs/SPEC_reachability_tolerance.md showed that this
is broken: on square H4 the level admitted at tol=1e-10 is symmetry-FORBIDDEN to the HF determinant
(FCI coefficient exactly zero) and its apparent amplitude is pure SCF convergence residue -- it
moves 19 orders of magnitude with ``PySCFDriver(conv_tol)``. No fixed constant fixes this: at
a = 1.190 A the same residue reaches 1.4e-8, above the looser threshold too.

The principled replacement is to ask which SPATIAL IRREP each eigenstate belongs to and keep only
those matching the HF determinant -- no amplitude threshold at all. The obvious objection is that
RHF breaks symmetry on exactly the strongly-correlated systems this matters for, and a
symmetry-broken determinant has no irrep to match.

That objection turns out to be self-consistent rather than fatal, which is the finding this module
encodes (specs/SPEC_symmetry_reachability.md, gated 9/9 on the square-H4 family):

    symmetric SCF   ->  the artifact is present (a forbidden level carries SCF residue)
                        AND the irrep labelling succeeds, so the filter can remove it.
    broken-sym SCF  ->  the reference is genuinely lower in energy and genuinely overlaps that
                        level (population ~0.45, not ~1e-10), so there is no artifact to remove
                        AND the irrep labelling correctly refuses.

The filter is therefore available exactly when it is needed. This module does NOT change the
DECISION PROCEDURE at the call sites -- swapping the threshold for the irrep filter is a separate
change with its own blast radius. What it does own is the single implementation: the certified arc
extracts its reachable sector through ``reachable_eigenpairs`` here, so when the filter does land
it is one edit rather than eight.
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

# Below this the two RHF solutions are the same one; above it the unsymmetrized solve found a
# genuinely lower, symmetry-broken determinant. 1e-6 Ha is far below the ~0.08 Ha gaps actually
# observed on square H4 and far above SCF convergence noise at conv_tol=1e-13.
SYMMETRY_BREAK_TOL = 1e-6

# Tight enough that the forbidden-state residue collapses to the machine-zero floor (~1e-28); at
# PySCFDriver's default 1e-9 the same residue sits at 5e-10 and is admitted as "reachable".
TIGHT_SCF_CONV_TOL = 1e-13

# The HF-reachable sector is defined by |<HF|psi_k>|^2 > tol. This is the CERTIFIED ARC's value,
# shared by certified_gaps.py, certified_dipole.py, certified_noise.py and
# hf_overlap_certificate.py -- but NOT by hf_overlap_subspace.py, which uses 1e-8. That divergence
# is not cosmetic: at square H4 a = 1.1 A the two thresholds select DIFFERENT ground states (HF
# overlap 2.25e-5 vs 0.667), so the d=1 and d=2 certificates whose head-to-head is
# SPEC_hf_overlap_subspace's headline are certifying different targets there. Named here so the
# divergence is visible and gated (tests/test_reachability_tolerance_spec.py); unifying it decides
# which of two recorded findings is right and is deliberately left to a follow-up.
# See specs/SPEC_reachability_tolerance.md.
REACHABLE_TOL_CERTIFIED = 1e-10


def scf_symmetry_status(atom: str, basis: str = "sto-3g",
                        conv_tol: float = TIGHT_SCF_CONV_TOL) -> Tuple[bool, float]:
    """(is_broken, dE) for the RHF reference: does the unsymmetrized solve find a LOWER solution?

    ``dE = E_unsymmetrized - E_symmetry_enforced``; negative beyond SYMMETRY_BREAK_TOL means RHF
    broke symmetry. On square H4 this is textbook behaviour for a strongly-correlated square, not a
    defect -- the broken solution is variationally better by ~0.08 Ha.
    """
    from pyscf import gto, scf

    free = scf.RHF(gto.M(atom=atom, basis=basis, symmetry=False, verbose=0))
    free.conv_tol = conv_tol
    free.kernel()
    enforced = scf.RHF(gto.M(atom=atom, basis=basis, symmetry=True, verbose=0))
    enforced.conv_tol = conv_tol
    enforced.kernel()
    d_e = float(free.e_tot - enforced.e_tot)
    return (d_e < -SYMMETRY_BREAK_TOL), d_e


def hf_orbital_irreps(atom: str, basis: str = "sto-3g",
                      conv_tol: float = TIGHT_SCF_CONV_TOL) -> Optional[list]:
    """Spatial irrep labels for the RHF orbitals, or None if the reference cannot be labelled.

    Returns None exactly when RHF has broken spatial symmetry, so the determinant is not a symmetry
    eigenfunction and no irrep assignment exists. That is a correct refusal, not a failure: see the
    module docstring for why those are also the cases with no artifact to remove.

    TWO DISTINCT FAILURES, which an earlier version of this conflated (LiH forced the distinction):

      * genuine symmetry BREAKING -- the unsymmetrized solve finds a strictly LOWER determinant
        (square H4 a=1.20: 0.076 Ha lower). No irrep exists. Refuse.
      * arbitrary rotation within a DEGENERATE irrep block -- same solution, same energy, but the
        free solve returns an arbitrary mixture of e.g. LiH's E1x/E1y pair, which cannot be labelled
        column-by-column. Nothing is broken; relabel using the symmetry-enforced solve, which is the
        same determinant in an irrep-adapted basis.

    Only the first is a refusal. Distinguishing them by ENERGY (not by whether labelling threw) is
    what makes the refusal meaningful.
    """
    from pyscf import gto, scf, symm

    broken, _ = scf_symmetry_status(atom, basis=basis, conv_tol=conv_tol)
    if broken:
        return None            # symmetry-broken reference -- no irrep to match against

    mol_sym = gto.M(atom=atom, basis=basis, symmetry=True, verbose=0)
    mf_sym = scf.RHF(mol_sym)
    mf_sym.conv_tol = conv_tol
    mf_sym.kernel()
    try:
        return list(symm.label_orb_symm(mol_sym, mol_sym.irrep_name, mol_sym.symm_orb,
                                        mf_sym.mo_coeff))
    except Exception:
        return None            # unlabelable even when symmetry-adapted -- refuse rather than guess


def symmetry_filter_available(atom: str, basis: str = "sto-3g",
                              conv_tol: float = TIGHT_SCF_CONV_TOL) -> bool:
    """Can a symmetry-aware reachability test be applied to this system at all?"""
    return hf_orbital_irreps(atom, basis=basis, conv_tol=conv_tol) is not None


def _dense_hf_projection(mh):
    """(eigenvalues ascending, eigenvectors, HF populations) -- dense, O(2^n)."""
    w, vecs = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    u = np.asarray(mh.hf_state().data, dtype=complex)
    return w, vecs, np.abs(vecs.conj().T @ u) ** 2


def hf_population_spectrum(mh) -> np.ndarray:
    """|<HF|psi_k>|^2 over the full eigenbasis (dense, O(2^n)) -- validation scale only."""
    return _dense_hf_projection(mh)[2]


def reachable_eigenpairs(mh, tol: float = REACHABLE_TOL_CERTIFIED):
    """(energies, eigenvectors) of the HF-reachable sector, ascending -- dense, O(2^n).

    THE one implementation of the ``|<HF|psi_k>|^2 > tol`` sector cut the certified arc runs on;
    energies are in the ELECTRONIC frame (add ``mh.energy_offset`` for totals). REFERENCE ONLY --
    it diagonalizes H exactly, so it is a validation oracle, never a live path.

    The threshold is the known-broken criterion this module's docstring dissects; it lives here so
    the eventual symmetry-filter replacement is one edit. Note the sector cut is what keeps a
    CHARGED species honest (HeH+): the global lowest eigenvector sits in a different
    particle-number sector, and only the reachable cut selects the state QKSD actually converges to.
    """
    w, vecs, pops = _dense_hf_projection(mh)
    keep = pops > tol
    return w[keep], vecs[:, keep]


if __name__ == "__main__":
    print(f"{'a':>5s} {'p0':>11s} {'SCF':>10s} {'irrep filter':>13s}")
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    for a in (1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40):
        geom = f"H 0 0 0; H {a} 0 0; H {a} {a} 0; H 0 {a} 0"
        p0 = float(hf_population_spectrum(build_molecular_hamiltonian(atom=geom))[0])
        broken, _ = scf_symmetry_status(geom)
        avail = symmetry_filter_available(geom)
        print(f"{a:5.2f} {p0:11.3e} {'broken' if broken else 'symmetric':>10s} "
              f"{'available' if avail else 'unavailable':>13s}")
