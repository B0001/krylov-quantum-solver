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

The filter is therefore available exactly when it is needed.

PER-EIGENSTATE DECISION (specs/SPEC_eigenstate_reachability.md, chem-obf). ``reachable_mask`` keeps
eigenvector k iff HF's component in k's (degenerate) eigenspace lies mostly in HF's exact symmetry
sector -- (N_alpha, N_beta) plus every point-group parity whose Z-string provably stays within 1/2
of the exact operator in ``mh``'s own MO basis -- AND its HF population clears the site's ``tol``.
The sector needs geometry:
``MolecularHamiltonian.build_args``, verified by rebuilding the operator. Without it the decision is
the old population cut. ``reachable_eigenpairs`` and the ODMD/MSD/Trotter centered frames run on it.
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
import scipy.linalg

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
# See specs/SPEC_reachability_tolerance.md. Since chem-obf the sites cut through ``reachable_mask``,
# which vetoes that symmetry-forbidden level at either tol; what the split still decides is only
# which faint symmetry-ALLOWED levels count (LiH: 4 in (1e-10, 1e-8]). SPEC_eigenstate_reachability.
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
    h = mh.qubit_hamiltonian.to_matrix()
    try:
        w, vecs = np.linalg.eigh(h)
    except np.linalg.LinAlgError:
        # ponytail: macOS Accelerate's complex ZHEEVD rejects its own LRWORK for N >~ 2900
        # (12+ qubits: "parameter number 10 had an illegal value", chem-a0y). ZHEEVR is
        # unaffected; used only on failure so Linux/OpenBLAS results stay bit-identical.
        w, vecs = scipy.linalg.eigh(h, driver="evr")
    u = np.asarray(mh.hf_state().data, dtype=complex)
    return w, vecs, np.abs(vecs.conj().T @ u) ** 2


def hf_population_spectrum(mh) -> np.ndarray:
    """|<HF|psi_k>|^2 over the full eigenbasis (dense, O(2^n)) -- validation scale only."""
    return _dense_hf_projection(mh)[2]


# A rebuild from ``build_args`` must reproduce the operator to this. Summation-order roundoff is
# ~1e-16; a different MO basis (another degenerate-pair rotation, another geometry) moves the
# coefficients by orders of magnitude more. A reproducibility check, not a physics threshold.
_SAME_OPERATOR_ATOL = 1e-12


def _driver_orbitals(mh):
    """(PySCF mol, active MO coefficients) of the run that built ``mh`` -- or None when ``mh`` has
    no ``build_args`` or rebuilding from them does not reproduce ``mh.qubit_hamiltonian`` (stale or
    edited operator). The identity term is ignored: a centered frame is the same operator."""
    args = getattr(mh, "build_args", None)
    if not args:
        return None
    from qiskit_nature.second_q.drivers import PySCFDriver
    from qiskit_nature.second_q.mappers import JordanWignerMapper
    from qiskit_nature.second_q.transformers import ActiveSpaceTransformer

    drv = PySCFDriver(atom=args["atom"], basis=args["basis"], charge=args["charge"],
                      spin=args["spin"], conv_tol=args["conv_tol"])
    problem = drv.run()
    mo = drv._calc.mo_coeff                 # private, but the driver has no public MO accessor
    n_e, n_o = args["active_electrons"], args["active_orbitals"]
    if n_e is not None:
        problem = ActiveSpaceTransformer(n_e, n_o).transform(problem)
        first = (drv._mol.nelectron - n_e) // 2          # the transformer's own active window
        mo = mo[:, first:first + n_o]
    if mh.qubit_hamiltonian.num_qubits != 2 * mo.shape[1]:
        return None
    rebuilt = JordanWignerMapper().map(problem.hamiltonian.second_q_op())
    diff = (rebuilt - mh.qubit_hamiltonian).simplify(atol=0)
    non_identity = diff.paulis.x.any(axis=1) | diff.paulis.z.any(axis=1)
    if np.abs(diff.coeffs[non_identity]).max(initial=0.0) > _SAME_OPERATOR_ATOL:
        return None
    return drv._mol, mo


def orbital_parity_bounds(mh):
    """[(odd, beta)] per parity bit of the point group's abelian irrep ids, or None if ``mh``'s MO
    basis cannot be verified (no geometry, or a stale/edited operator).

    ``odd`` marks the active MOs the bit's Z-string treats as odd -- the sign of the diagonal of the
    exact symmetry operation g_b written in ``mh``'s own MO basis. ``beta`` bounds
    ||Z-string - exact operator|| on N-electron states: the sum of the N largest spin-orbital
    |1 - lambda| over the eigenvalues of D_b R_b (specs/SPEC_eigenstate_reachability.md section 3).
    Rigorous for the full orbital space, where D_b R_b is orthogonal; heuristic for active spaces.
    """
    found = _driver_orbitals(mh)
    if found is None:
        return None
    from pyscf import gto, symm

    mol, mo = found
    ref = gto.M(atom=mol.atom, unit=mol.unit, basis=mol.basis, charge=mol.charge, spin=mol.spin,
                symmetry=True, verbose=0)
    # symmetry-adapted AOs in the DRIVER's frame (the symmetric mol is re-oriented)
    so, ids = symm.symm_adapted_basis(mol, ref.groupname, ref._symm_orig, ref._symm_axes)
    ids = np.asarray(ids) % 10              # linear groups (Dooh/Coov): the D2h/C2v parent irrep
    s = mol.intor_symmetric("int1e_ovlp")
    proj = [mo.T @ s @ c @ np.linalg.solve(c.T @ s @ c, c.T @ s @ mo) for c in so]
    n_elec = sum(mh.num_particles)
    bounds = []
    for b in range(int(ids.max()).bit_length()):
        r = sum((-1.0) ** ((i >> b) & 1) * p for i, p in zip(ids, proj))   # g_b, MO basis
        odd = np.diag(r) < 0
        lam = np.linalg.eigvals(np.where(odd, -1.0, 1.0)[:, None] * r)
        dev = np.sort(np.repeat(np.abs(1.0 - lam), 2))[::-1]               # alpha and beta
        bounds.append((odd, float(dev[:n_elec].sum())))
    return bounds


def hf_symmetry_sector(mh) -> Optional[np.ndarray]:
    """Boolean mask over the 2^n computational basis: HF's exact symmetry sector.

    Same (N_alpha, N_beta) as the HF bitstring -- always exact -- plus the same parity under every
    point-group bit of the greedy smallest-beta set with sum(beta) < 1, which keeps the product
    projector within 1/2 of the exact one so the majority test in ``symmetry_allowed`` cannot flip
    for a non-degenerate symmetry eigenstate. None when |HF> is not a Jordan-Wigner HF bitstring
    (other mapper, tapered operator): no sector is formed and callers keep the population cut.
    """
    cached = mh.__dict__.get("_hf_sector")
    if cached is not None and cached[0] is mh.qubit_hamiltonian:
        return cached[1]        # ponytail: per-instance cache; the rebuild is an SCF + a JW map
    n, n_orb = mh.num_qubits, mh.num_spatial_orbitals
    n_a, n_b = mh.num_particles
    hf = int(np.argmax(np.abs(np.asarray(mh.hf_state().data))))
    sector = None
    if n == 2 * n_orb and hf == (1 << n_a) - 1 + (((1 << n_b) - 1) << n_orb):
        occ = ((np.arange(2 ** n)[:, None] >> np.arange(n)) & 1).astype(np.int8)
        sector = (occ[:, :n_orb].sum(1) == n_a) & (occ[:, n_orb:].sum(1) == n_b)
        budget = 0.0
        for odd, beta in sorted(orbital_parity_bounds(mh) or [], key=lambda t: t[1]):
            if budget + beta >= 1.0:
                break
            budget += beta
            parity = (occ[:, :n_orb] @ odd + occ[:, n_orb:] @ odd) % 2
            sector &= parity == parity[hf]
    mh.__dict__["_hf_sector"] = (mh.qubit_hamiltonian, sector)
    return sector


# Eigenvalues closer than this are one eigenspace. eigh splits exactly degenerate levels (open-shell
# M_s pairs, Pi/E pairs) by ~1e-14 here; the witnesses' symmetry gaps are ~0.1 Ha. Numerical, and
# fail-safe: a wider cluster only falls back toward the population cut, it never drops a level.
_DEGENERATE_ATOL = 1e-9


def symmetry_allowed(mh, w, vecs) -> np.ndarray:
    """Per eigenvector: does HF's component in its eigenspace lie mostly in HF's symmetry sector?

    Decided per cluster C of numerically degenerate eigenvalues (``w`` ascending) on the
    basis-independent projection phi_C = V_C V_C^dagger |HF>: allowed iff ||P_sector phi_C||^2 >
    ||phi_C||^2 / 2 -- exactly 1 or 0 under exact symmetry, so 1/2 is its rounding, not a tuned
    constant. Per vector this is the majority test; per cluster it survives ``eigh`` mixing exactly
    degenerate levels across sectors (open-shell M_s pairs), which a per-vector test does not.
    All True when no sector can be formed (the old behaviour)."""
    sector = hf_symmetry_sector(mh)
    allowed = np.ones(vecs.shape[1], dtype=bool)
    if sector is None:
        return allowed
    hf = int(np.argmax(np.abs(np.asarray(mh.hf_state().data))))
    for idx in np.split(np.arange(len(w)), np.flatnonzero(np.diff(w) > _DEGENERATE_ATOL) + 1):
        c = vecs[hf, idx].conj()                         # <psi_k|HF>, k in the cluster
        inside = vecs[np.ix_(sector, idx)] @ c           # P_sector phi_C
        allowed[idx] = np.vdot(inside, inside).real > 0.5 * np.vdot(c, c).real
    return allowed


def reachable_mask(mh, w, vecs, pops, tol: float) -> np.ndarray:
    """THE HF-reachability decision: symmetry-allowed AND populated above the site's ``tol``."""
    return symmetry_allowed(mh, w, vecs) & (np.asarray(pops) > tol)


def reachable_eigenpairs(mh, tol: float = REACHABLE_TOL_CERTIFIED):
    """(energies, eigenvectors) of the HF-reachable sector, ascending -- dense, O(2^n).

    THE one implementation of the reachable sector the certified arc runs on: ``reachable_mask``
    (symmetry-allowed AND ``|<HF|psi_k>|^2 > tol``; just the population cut when ``mh`` carries no
    verifiable geometry). Energies are in the ELECTRONIC frame (add ``mh.energy_offset`` for
    totals). REFERENCE ONLY -- it diagonalizes H exactly, so it is a validation oracle, never a live
    path. Note the sector cut is what keeps a CHARGED species honest (HeH+): the global lowest
    eigenvector sits in a different particle-number sector, and only the reachable cut selects the
    state QKSD actually converges to.
    """
    w, vecs, pops = _dense_hf_projection(mh)
    keep = reachable_mask(mh, w, vecs, pops, tol)
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
