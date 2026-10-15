#!/usr/bin/env python3
"""
reachability.py -- which eigenstates a Hartree-Fock reference reaches, plus a symmetry diagnostic.

The certified arc and the ODMD/MSD family define the "HF-reachable" sector by thresholding the HF
population, ``|<HF|psi_k>|^2 > tol``. specs/SPEC_reachability_tolerance.md showed what that admits
on square H4: a level symmetry-FORBIDDEN to the exact HF determinant (FCI coefficient exactly
zero), whose apparent amplitude is pure SCF convergence residue -- it moves 19 orders of magnitude
with ``PySCFDriver(conv_tol)``, and at a = 1.190 A it reaches 1.4e-8, above the looser cut too.

Vetoing that level is nevertheless WRONG for a reference (specs/SPEC_eigenstate_reachability.md,
G9, which killed chem-obf's first version): the qubit Hamiltonian is built from the loose SCF's
orbitals, in which |HF> really does carry the residue, and QKSD from |HF> converges to that level,
~150 mHa below the symmetric ground state. A reference must describe what the data targets, so
``reachable_mask`` IS the population cut. The residue's remedy is the SCF: rebuild with
``conv_tol=TIGHT_SCF_CONV_TOL``.

What survives is a best-effort DIAGNOSTIC. ``symmetry_allowed`` asks per eigenstate whether HF's
component in its eigenspace lies mostly in HF's exact symmetry sector -- (N_alpha, N_beta) plus
every point-group parity whose Z-string provably stays within 1/2 of the exact operator in ``mh``'s
own MO basis, which needs geometry (``MolecularHamiltonian.build_args``, verified by rebuilding the
operator). ``reachable_mask`` warns when the cut admits a level it judges forbidden; the diagnostic
never changes the mask and never raises. Its known blind spots: the spec's section 8.

The per-SYSTEM question (specs/SPEC_symmetry_reachability.md, ``symmetry_filter_available``): with
a symmetric SCF the residue is present AND the irrep labelling succeeds; with a broken-symmetry SCF
the reference genuinely overlaps that level (population ~0.45, not ~1e-10) AND labelling refuses.
"""
from __future__ import annotations

import warnings
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
# See specs/SPEC_reachability_tolerance.md. The sites cut through ``reachable_mask`` (this cut plus
# a residue warning); SPEC_eigenstate_reachability G9 is why a symmetry veto cannot settle it.
REACHABLE_TOL_CERTIFIED = 1e-10


def scf_symmetry_status(atom: str, basis: str = "sto-3g",
                        conv_tol: float = TIGHT_SCF_CONV_TOL) -> Tuple[bool, float]:
    """(is_broken, dE) for the RHF reference: does the unsymmetrized solve find a LOWER solution?

    ``dE = E_unsymmetrized - E_symmetry_enforced``; negative beyond SYMMETRY_BREAK_TOL means RHF
    broke symmetry. On square H4 this is textbook behaviour for a strongly-correlated square, not a
    defect -- the broken solution is variationally better by ~0.08 Ha.

    Both RHF solves run with ``pyscf.lib.num_threads(1)`` pinned (ambient count restored in a
    ``finally``), the same fix chem-7cx applied to ``build_molecular_hamiltonian``: PySCF's
    multi-threaded integral evaluation is not bit-associative, which on a near-degenerate geometry
    can tip the free (``symmetry=False``) solve onto a different locally-stable branch between
    otherwise-identical calls (chem-s6c, specs/SPEC_reachability_scf_determinism.md).
    """
    from pyscf import gto, lib, scf

    ambient_threads = lib.num_threads()
    lib.num_threads(1)
    try:
        free = scf.RHF(gto.M(atom=atom, basis=basis, symmetry=False, verbose=0))
        free.conv_tol = conv_tol
        free.kernel()
        enforced = scf.RHF(gto.M(atom=atom, basis=basis, symmetry=True, verbose=0))
        enforced.conv_tol = conv_tol
        enforced.kernel()
    finally:
        lib.num_threads(ambient_threads)
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
    from pyscf import gto, lib, scf, symm

    broken, _ = scf_symmetry_status(atom, basis=basis, conv_tol=conv_tol)
    if broken:
        return None            # symmetry-broken reference -- no irrep to match against

    mol_sym = gto.M(atom=atom, basis=basis, symmetry=True, verbose=0)
    mf_sym = scf.RHF(mol_sym)
    mf_sym.conv_tol = conv_tol
    # Same thread-pin as scf_symmetry_status (chem-s6c): this is a second, independent RHF solve
    # (symmetry-enforced, but still subject to the same non-bit-associative integral evaluation).
    ambient_threads = lib.num_threads()
    lib.num_threads(1)
    try:
        mf_sym.kernel()
    finally:
        lib.num_threads(ambient_threads)
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


# A rebuild from ``build_args``, AFTER its MO gauge is reconciled with ``mh`` (see
# _ENERGY_MATCH_ATOL/_MAX_SIGN_SEARCH_ORBITALS below), must reproduce the operator to this.
# Summation-order roundoff alone is ~1e-16 (measured: two gauge-matched rebuilds at
# conv_tol=TIGHT_SCF_CONV_TOL=1e-13 differ by ~1e-16 to 3.5e-15, chem-g14). But at the conv_tol
# values this module's gate file actually exercises (1e-6, 1e-9, the builder default), independent
# SCF solves land on genuinely different last-few-digits of the SAME MO basis, not just roundoff:
# measured max residual after gauge reconciliation is 8.1e-5 at conv_tol=1e-6 and 2.4e-6 at 1e-9
# (square H4, a in {1.05, 1.10, 1.19, 1.35}). 1e-4 sits above that with >10x margin and more than
# three orders of magnitude below the smallest genuinely-different-solution residual observed
# (~0.067, see _ENERGY_MATCH_ATOL) -- a reproducibility check at the SCF's own convergence floor,
# not a physics threshold.
_SAME_OPERATOR_ATOL = 1e-4

# PySCFDriver hardcodes symmetry=False, so a rebuild's free RHF is not pinned to the original run's
# branch (chem-g14). Two failure modes, told apart by energy: a genuinely different, WORSE local
# minimum (gap >= ~0.067 Ha on square H4) -- reject and retry -- or an exactly (or near-exactly)
# energy-degenerate alternate MO gauge (gap as small as ~1e-11 at conv_tol=1e-9 -- a real
# degeneracy, not SCF noise, consistent with PySCF's symmetry detector resolving square H4's true
# D4h only as far as D2h). Measured (chem-g14): the unflipped rebuild of a gauge-degenerate branch
# still differs from ``mh`` by ~0.068 (the SAME size as a genuinely different solution!) in most
# such cases, but searching the discrete per-orbital sign gauge (A^T S B, the AO metric, not a plain
# dot product) always finds a pattern bringing it down to the ordinary noise floor above --
# confirming the degeneracy is a discrete sign choice, never a continuous rotation, once searched
# correctly. 1e-6 Ha sits well above the conv_tol=1e-9 energy-noise floor (~1e-11) and four orders
# of magnitude below the smallest observed wrong-branch energy gap.
_ENERGY_MATCH_ATOL = 1e-6

# 2**(k-1) sign patterns per retry (one column's sign is a global gauge redundancy for both the one-
# and two-electron integrals, so it is fixed rather than searched). No system in this module's own
# gate file needs more than 6 active orbitals; above the bound, fall through to None as before --
# an honest refusal rather than an unbounded search.
_MAX_SIGN_SEARCH_ORBITALS = 10
_MAX_REBUILD_RETRIES = 5


def _driver_orbitals(mh):
    """(PySCF mol, active MO coefficients) of the run that built ``mh`` -- or None when ``mh`` has
    no ``build_args`` or rebuilding from them does not reproduce ``mh.qubit_hamiltonian`` (stale or
    edited operator). The identity term is ignored: a centered frame is the same operator.

    Retries the SCF (bounded) past a different, worse local minimum, and -- once the energy
    matches -- searches the discrete per-orbital sign gauge (bounded) past an exactly degenerate
    alternate MO basis. See ``_ENERGY_MATCH_ATOL`` / ``_MAX_SIGN_SEARCH_ORBITALS`` above and
    chem-g14 / specs/SPEC_rebuild_gauge_search.md. A rebuild whose plain, unflipped comparison
    already matches costs exactly one PySCFDriver call, same as before this fix.

    Each retry's ``drv.run_pyscf()`` runs with ``pyscf.lib.num_threads(1)`` pinned (ambient count
    restored in a ``finally``) -- the same chem-7cx/chem-s6c fix for the same mechanism (PySCF's
    multi-threaded integral evaluation is not bit-associative), applied here to this function's own
    independent RHF solve (chem-50w, specs/SPEC_driver_orbitals_scf_determinism.md). The sign search
    itself makes no new SCF call (just a fresh AO2MO transform on already-computed ``mo_coeff``), so
    nothing else here touches thread count.
    """
    args = getattr(mh, "build_args", None)
    if not args:
        return None
    from pyscf import lib
    from qiskit_nature.second_q.drivers import PySCFDriver
    from qiskit_nature.second_q.mappers import JordanWignerMapper
    from qiskit_nature.second_q.transformers import ActiveSpaceTransformer

    mapper = JordanWignerMapper()
    n_e, n_o = args["active_electrons"], args["active_orbitals"]
    hf_energy = getattr(mh, "hf_energy", None)

    def rebuild_from(drv):
        """(operator or None, active MO slice, active window's start column) for drv's CURRENT
        drv._calc.mo_coeff -- re-reads it fresh, so callers may flip signs in place first."""
        problem = drv.to_problem()
        mo = drv._calc.mo_coeff
        window_start = 0
        if n_e is not None:
            problem = ActiveSpaceTransformer(n_e, n_o).transform(problem)
            # the transformer's own active window; a tuple (n_alpha, n_beta) is summed, as it does
            window_start = (drv._mol.nelectron - int(np.sum(n_e))) // 2
            mo = mo[:, window_start:window_start + n_o]
        if mh.qubit_hamiltonian.num_qubits != 2 * mo.shape[1]:
            return None, mo, window_start
        return mapper.map(problem.hamiltonian.second_q_op()), mo, window_start

    def same_operator(op):
        if op is None:
            return False
        diff = (op - mh.qubit_hamiltonian).simplify(atol=0)
        non_identity = diff.paulis.x.any(axis=1) | diff.paulis.z.any(axis=1)
        return np.abs(diff.coeffs[non_identity]).max(initial=0.0) <= _SAME_OPERATOR_ATOL

    for _ in range(_MAX_REBUILD_RETRIES):
        drv = PySCFDriver(atom=args["atom"], basis=args["basis"], charge=args["charge"],
                          spin=args["spin"], conv_tol=args["conv_tol"])
        # PySCFDriver hardcodes symmetry=False; its RHF solve shares chem-7cx/chem-s6c's
        # non-bit-associative multi-threaded integral evaluation (chem-50w). Pin/restore around
        # just this call -- the sign search below makes no further SCF call.
        ambient_threads = lib.num_threads()
        lib.num_threads(1)
        try:
            drv.run_pyscf()
        finally:
            lib.num_threads(ambient_threads)
        if hf_energy is not None and abs(drv._calc.e_tot - hf_energy) > _ENERGY_MATCH_ATOL:
            continue                            # a different, worse SCF local minimum -- retry

        base_mo = drv._calc.mo_coeff.copy()
        op, mo, window_start = rebuild_from(drv)
        if same_operator(op):
            return drv._mol, mo
        if op is None:
            continue                            # qubit-count mismatch -- not a gauge question

        k = mo.shape[1]
        if k - 1 > _MAX_SIGN_SEARCH_ORBITALS:
            continue                            # too large to search; stay with the plain refusal
        for pattern in range(1, 1 << (k - 1)):
            signs = np.ones(k)
            signs[1:] = [1.0 if (pattern >> b) & 1 == 0 else -1.0 for b in range(k - 1)]
            flipped = base_mo.copy()
            flipped[:, window_start:window_start + k] *= signs
            drv._calc.mo_coeff = flipped
            op, mo, _ = rebuild_from(drv)
            if same_operator(op):
                return drv._mol, mo
    return None


# D_b R_b bounds the Z-string error only where R_b (the point-group operation g_b restricted to the
# active MO window, oblique-projected) is itself orthogonal. An active window that splits a
# point-group-degenerate shell compresses R_b into something far from orthogonal, and the beta bound
# built from it can then be small (and so get greedily ACCEPTED by hf_symmetry_sector) while being
# physically wrong -- chem-648 L1, CH4 Td CAS(2,2): a populated level solve_excited reports as theta_1
# got flagged forbidden this way. Measured separation on that witness: ||R_b R_b^T - I|| ~ 1e-14 for
# the full orbital space or a window that does not split a degenerate shell, >= 5.7e-4 once it does.
_R_ORTHOGONALITY_ATOL = 1e-6

# Bit -> the Cartesian operation it is the sign of, derived from pyscf.symm.param.CHARACTER_TABLE
# (column order matches OPERATOR_TABLE; bit b's column is the one whose +-1 pattern across irrep ids
# equals (-1)**((id >> b) & 1)). Only the finite abelian point groups symm_adapted_basis actually
# returns for a non-linear molecule; Dooh/Coov (linear) are out of scope for the geometric check
# below -- their ids are folded into a D2h/C2v parent by ``% 10`` and the bit layout differs.
_BIT_OPERATION = {
    "D2h": ("C2y", "C2z", "i"),
    "D2": ("C2y", "C2z"),
    "C2h": ("C2z", "i"),
    "C2v": ("sy", "C2z"),
    "Cs": ("sz",),
    "Ci": ("i",),
    "C2": ("C2z",),
    "C1": (),
}

# A real symmetry operation applied to the ACTUAL (unsymmetrized) geometry, independent of whatever
# geometric tolerance PySCF's own detection used. Exact-to-double-precision geometries land at
# ~1e-15 Bohr; chem-648 L2's witness (square H4, one atom moved 1e-6 Angstrom -- still inside PySCF's
# natm-scaled detection tolerance, so ``ref.groupname`` still reports D2h) lands at ~1e-6 Bohr under
# the operations that atom's shift actually breaks. Fixed and size-independent, unlike
# ``symm.geom.TOLERANCE`` (which both PySCF's detection AND symm_adapted_basis's own atom-equivalence
# grouping use, scaled by ``1/sqrt(1+natm)``) -- that is the point: this check does not inherit it.
_GEOMETRY_RESIDUAL_ATOL = 1e-8


def _geometry_residuals(mol, ref):
    """{op_name: max atom-mapping residual (Bohr)} for each bit's claimed operation, checked against
    the REAL geometry directly -- not PySCF's own (size-scaled) detection tolerance."""
    from pyscf import symm

    ops = _BIT_OPERATION.get(ref.groupname, ())
    if not ops:
        return {}
    prop_atoms = mol.format_atom(mol._atom, ref._symm_orig, ref._symm_axes, "Bohr")
    coords = np.array([a[1] for a in prop_atoms])
    elems = [a[0] for a in prop_atoms]
    opdic = symm.geom.symm_ops(ref.groupname)
    out = {}
    for name in ops:
        op = opdic[name]
        newc = coords * op if np.ndim(op) == 0 else coords @ op
        residual = 0.0
        for i, elem in enumerate(elems):
            same = [j for j, e in enumerate(elems) if e == elem]
            residual = max(residual, float(np.min(np.linalg.norm(coords[same] - newc[i], axis=1))))
        out[name] = residual
    return out


def orbital_parity_bounds(mh):
    """[(odd, beta)] per parity bit of the point group's abelian irrep ids, or None if ``mh``'s MO
    basis cannot be verified (no geometry, or a stale/edited operator).

    ``odd`` marks the active MOs the bit's Z-string treats as odd -- the sign of the diagonal of the
    exact symmetry operation g_b written in ``mh``'s own MO basis. ``beta`` bounds
    ||Z-string - exact operator|| on N-electron states: the sum of the N largest spin-orbital
    |1 - lambda| over the eigenvalues of D_b R_b (specs/SPEC_eigenstate_reachability.md section 3).
    Rigorous for the full orbital space, where D_b R_b is orthogonal; heuristic for active spaces
    (and wrong when the window splits a degenerate shell) -- guarded by an explicit orthogonality
    check on R_b, below which a bit is excluded (beta = inf) rather than trusted (chem-648 L1). A bit
    whose claimed operation does not actually map the real geometry onto itself to
    ``_GEOMETRY_RESIDUAL_ATOL``, independent of PySCF's own detection tolerance, is excluded the same
    way (chem-648 L2). May raise where PySCF only half-detects the symmetry; ``reachable_mask`` treats
    any failure here as "no diagnostic".
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
    geom_residuals = _geometry_residuals(mol, ref)
    op_names = _BIT_OPERATION.get(ref.groupname, ())
    bounds = []
    for b in range(int(ids.max()).bit_length()):
        r = sum((-1.0) ** ((i >> b) & 1) * p for i, p in zip(ids, proj))   # g_b, MO basis
        odd = np.diag(r) < 0
        untrustworthy = np.linalg.norm(r @ r.T - np.eye(r.shape[0])) > _R_ORTHOGONALITY_ATOL
        if b < len(op_names):
            untrustworthy |= geom_residuals.get(op_names[b], 0.0) > _GEOMETRY_RESIDUAL_ATOL
        if untrustworthy:
            bounds.append((odd, float("inf")))   # compressed/not-really-symmetric: not trustworthy
            continue
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
    (other mapper, tapered operator): no sector, so no diagnostic. The layout check sees only the HF
    index, which an interleaved-mapper operator CAN share with a Jordan-Wigner one for some
    (n_alpha, n_beta, n_orb) -- chem-648 L3 -- so it additionally requires ``build_args`` (when
    present) to record ``mapper="JordanWignerMapper"``; absent ``build_args`` (no provenance, e.g.
    ``build_hamiltonian_from_integrals``, which only ever uses Jordan-Wigner) still trusts the index.
    """
    cached = mh.__dict__.get("_hf_sector")
    if cached is not None and cached[0] is mh.qubit_hamiltonian:
        return cached[1]        # ponytail: per-instance cache; the rebuild is an SCF + a JW map
    build_args = getattr(mh, "build_args", None)
    mapper_name = build_args.get("mapper", "JordanWignerMapper") if build_args else "JordanWignerMapper"
    n, n_orb = mh.num_qubits, mh.num_spatial_orbitals
    n_a, n_b = mh.num_particles
    hf = int(np.argmax(np.abs(np.asarray(mh.hf_state().data))))
    sector = None
    if (mapper_name == "JordanWignerMapper" and n == 2 * n_orb
            and hf == (1 << n_a) - 1 + (((1 << n_b) - 1) << n_orb)):
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
# M_s pairs, Pi/E pairs) by ~1e-14 here; the witnesses' symmetry gaps are ~0.1 Ha. Numerical. A
# cluster shares one majority verdict, so levels of different symmetry closer than this can get a
# wrong one -- a missed or spurious diagnostic warning; the reachable mask never depends on it.
_DEGENERATE_ATOL = 1e-9

# A cluster's sector weight is 0 or 1 to within this on every witnessed gated system (also G6's own
# bound) -- real symmetry (even a true degeneracy spanning two sectors, e.g. an open-shell M_s pair)
# concentrates HF's overlap on one side. A weight meaningfully away from 0 or 1 is the signature of
# chem-648 L4: an ACCIDENTAL <1e-9 Ha degeneracy merging two different-symmetry, both-populated
# levels into one cluster, where the majority rule's 1/2 cut is no longer backed by physics. Abstain
# (do not flag) rather than guess which side the majority verdict actually belongs to.
_AMBIGUOUS_WEIGHT_MARGIN = 1e-3


def symmetry_allowed(mh, w, vecs) -> np.ndarray:
    """Per eigenvector: does HF's component in its eigenspace lie mostly in HF's symmetry sector?

    Decided per cluster C of numerically degenerate eigenvalues (``w`` ascending) on the
    basis-independent projection phi_C = V_C V_C^dagger |HF>: allowed iff ||P_sector phi_C||^2 >
    ||phi_C||^2 / 2 -- exactly 1 or 0 under exact symmetry, so 1/2 is its rounding, not a tuned
    constant. Per vector this is the majority test; per cluster it survives ``eigh`` mixing exactly
    degenerate levels across sectors (open-shell M_s pairs), which a per-vector test does not. A
    cluster whose weight falls in the ambiguous band around 1/2 (chem-648 L4: an accidental,
    non-symmetry degeneracy merging two real sectors) abstains -- treated as allowed -- rather than
    applying the majority cut, since neither side of 1/2 is trustworthy there.
    Meaningful only where HF's projection is above roundoff (populated levels). All True when no
    sector can be formed. A DIAGNOSTIC: no reference depends on it (see ``reachable_mask``)."""
    sector = hf_symmetry_sector(mh)
    allowed = np.ones(vecs.shape[1], dtype=bool)
    if sector is None:
        return allowed
    hf = int(np.argmax(np.abs(np.asarray(mh.hf_state().data))))
    for idx in np.split(np.arange(len(w)), np.flatnonzero(np.diff(w) > _DEGENERATE_ATOL) + 1):
        c = vecs[hf, idx].conj()                         # <psi_k|HF>, k in the cluster
        inside = vecs[np.ix_(sector, idx)] @ c           # P_sector phi_C
        c_norm = np.vdot(c, c).real
        weight = np.vdot(inside, inside).real / c_norm if c_norm > 0 else 0.0
        if _AMBIGUOUS_WEIGHT_MARGIN < weight < 1.0 - _AMBIGUOUS_WEIGHT_MARGIN:
            continue                                    # accidental cross-sector degeneracy: abstain
        allowed[idx] = weight > 0.5
    return allowed


def reachable_mask(mh, w, vecs, pops, tol: float) -> np.ndarray:
    """THE HF-reachability decision: ``|<HF|psi_k>|^2 > tol`` -- residue levels included, because
    QKSD from |HF> can converge to them (and does to the lowest: SPEC_eigenstate_reachability G9).
    The symmetry decision runs only as a diagnostic: it warns when the cut admits a level it judges
    forbidden, and never changes the mask."""
    pops = np.asarray(pops)
    keep = pops > tol
    try:
        flagged = np.flatnonzero(keep & ~symmetry_allowed(mh, w, vecs))
    except Exception:
        # ponytail: broad catch is fine -- diagnostic only, the reference path never depends on it
        # (PySCF raises on geometries it half-detects as symmetric; that just means no diagnostic).
        flagged = ()
    if len(flagged):
        levels = ", ".join(f"E = {float(np.real(w[k])):.6f} Ha (electronic) at "
                           f"|<HF|psi>|^2 = {float(pops[k]):.2e}" for k in flagged)
        warnings.warn(f"HF-reachable level(s) {levels} look symmetry-forbidden to |HF>: that "
                      f"population is SCF residue. The reference keeps them because QKSD from |HF> "
                      f"can converge to them; rebuild with conv_tol=TIGHT_SCF_CONV_TOL "
                      f"({TIGHT_SCF_CONV_TOL:g}) to remove the residue.",
                      RuntimeWarning, stacklevel=2)
    return keep


def reachable_eigenpairs(mh, tol: float = REACHABLE_TOL_CERTIFIED):
    """(energies, eigenvectors) of the HF-reachable sector, ascending -- dense, O(2^n).

    THE one implementation of the reachable sector the certified arc runs on: ``reachable_mask``,
    i.e. ``|<HF|psi_k>|^2 > tol`` with a warning on SCF-residue levels (kept: QKSD from |HF> can
    converge to them). Energies are in the ELECTRONIC frame (add ``mh.energy_offset`` for
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
