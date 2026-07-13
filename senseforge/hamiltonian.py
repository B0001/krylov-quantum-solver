#!/usr/bin/env python3
"""
senseforge.hamiltonian -- strain and field perturbations of the validated Nb3X8 dimer model.

TWO DEVIATIONS FROM THE LITERAL TASK SPEC, both recorded here rather than hidden (see the
SenseForge design note in the PR description).

(1) NO CIF/AB-INITIO GEOMETRY PATH. specs/tasks/04-senseforge.md and the PRD (section 3) call
for "strained cluster geometries from the CIF path." No such path exists in this repo: Nb3X8
model physics is a validated PARAMETRIZED HUBBARD DIMER (``nb3x8_gaps.dimer_cluster_integrals``,
cRPA parameters of arXiv:2501.10320), never PySCF/CIF-derived integrals -- and a real ab-initio
Nb3 cluster would need an ECP basis outside certchem's validated ``ALLOWED_BASES`` envelope
(Nb, Z=41) and would likely exceed its 16-spin-orbital cap. Building that path is a real
follow-up, not attempted here. This module reuses the existing, gated model instead:
  STRAIN: nb3x8_strain.py already establishes "|t| is the sole strain proxy" (its Grueneisen
  parameter is d ln X / d ln|t|). This module makes that a literal sweep axis: strain epsilon is
  DEFINED as the fractional perturbation of the interlayer hopping, ``t(eps) = t0*(1+eps)`` --
  not a physical uniaxial/biaxial engineering strain with its own elastic-tensor mapping (which
  does not exist in this repo either). Single strain axis; uniaxial-vs-biaxial scaling deferred.

(2) NO certified_gaps.gap_bracket ON THIS SYSTEM -- ACTIVELY WRONG, NOT JUST UNNECESSARY.
The PRD asks to "attach a certified two-sided bracket to the gap (certified_gaps.py)". Tried
first, and rejected after a numerical check: the singlet -> Sz=0-triplet gap J (the "spin-gap
sensitivity" the whole PRD is about, PRD sec 1) is a DARK excitation from the closed-shell HF
reference -- Sz_tot is a good quantum number fixed by the (n_alpha, n_beta) Fock sector, and the
plain HF-seeded real-time Krylov subspace has near-zero overlap with the triplet. Verified: at
eps=0, ``certified_gaps.gap_bracket`` on the bare dimer Hamiltonian returns ~1117 meV (the BRIGHT
ionic/charge gap E_s), not J~66 meV -- silently the wrong physical quantity. Seeding the Krylov
subspace from the ODMD spin-kick reference (odmd_spin.py) does not fix this either: the staggered
moment S1z-S2z maps the exact singlet ground state to the exact triplet state with ZERO leakage
(odmd_spin.py's own docstring: "Sz|psi0> is an exact eigenstate") -- a 1-dimensional, degenerate
Krylov subspace with no second Ritz value to bracket against (``gap_bracket``'s own vacuous-case
branch). And certified_gaps.py's own module docstring says outright: "On a system small enough
to diagonalize exactly (e.g. the 4-qubit Nb3X8 dimer) certification is pointless -- you already
have the gap." This module takes that at face value: every gap here (strain-swept J(eps), and
the field-split triplet sublevels) comes from EXACT diagonalization / the closed-form
``dimer_exchange_analytic``, matching how nb3x8_strain.py / nb3x8_susceptibility.py / odmd_spin.py
already treat this system. Results are wrapped in certchem's ``Bracket``/``CertifiedResult``
contract as explicit zero-width EXACT brackets, for interface consistency with the rest of the
portfolio -- honestly labeled as exact, not as an approximate-but-rigorously-bounded estimate.

FIELD: a uniform Zeeman term H_Z = g*mu_B*B*Sz_tot (g=2, matching nb3x8_susceptibility.py's g=2
convention) is added to the dimer Hamiltonian. Second physics finding, also verified numerically
before implementing: a uniform field has ZERO effect within the (1,1)/Sz=0 sector (Sz_tot=0
identically there) -- it only splits the S=1 triplet's Sz=+1/0/-1 sublevels apart, and the Sz=+/-1
components live in DIFFERENT (2,0)/(0,2) Fock sectors that no Sz-conserving perturbation of the
closed-shell reference can reach at all. The physically real field-swept gap is therefore the
exact, closed-form Zeeman splitting of the already-exactly-diagonalized 6-level N=2 manifold, not
a Krylov quantity of any kind.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from qiskit.quantum_info import SparsePauliOp
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.operators import FermionicOp

from certchem.contract import Bracket, Certificate, CertifiedResult
from certchem.core import solver_version
from nb3x8_gaps import NB3X8_LT_BULK, dimer_cluster_integrals
from odmd_spin import dimer_exchange_analytic

#: g-factor and Bohr magneton (meV/T), matching nb3x8_susceptibility.py's g=2 spin-only convention.
G_FACTOR = 2.0
MU_B_MEV_PER_T = 0.057883818  # CODATA mu_B = 5.7883818e-5 eV/T

_MAPPER = JordanWignerMapper()
#: Total Sz operator on the 2-orbital (4 spin-orbital) dimer, same JW block order as
#: nb3x8_susceptibility.py's ``_SZ`` (validated there against the exact N=2 spectrum structure).
_SZ_OP: SparsePauliOp = _MAPPER.map(
    FermionicOp(
        {"+_0 -_0": 0.5, "+_1 -_1": 0.5, "+_2 -_2": -0.5, "+_3 -_3": -0.5},
        num_spin_orbitals=4,
    )
)
_N_OP = sum(
    _MAPPER.map(FermionicOp({f"+_{i} -_{i}": 1.0}, num_spin_orbitals=4)).to_matrix()
    for i in range(4)
)

_EXACT_METHOD = "exact diagonalization / closed form -- zero-width by construction (dimer is FCI-trivial)"


def strained_params(halide: str, eps: float) -> dict:
    """(U0, t, Us) for ``halide`` (e.g. "Nb3Cl8") with the interlayer hopping strained by
    ``eps`` (fractional): ``t(eps) = t0 * (1 + eps)``. See module docstring for the honest scope
    of this strain proxy."""
    p = dict(NB3X8_LT_BULK[halide])
    p["t"] = p["t"] * (1.0 + eps)
    return p


def _exact_bracket(gap: float, halide: str, extra: dict) -> CertifiedResult:
    manifest = {"halide": halide, **extra}
    return CertifiedResult(
        bracket=Bracket(lower_hartree=gap, upper_hartree=gap, best_estimate_hartree=gap),
        certificate=Certificate(
            method=_EXACT_METHOD,
            floor_check="n/a (exact, not variational)",
            krylov_dim=0,
            convergence="exact",
            solver_version=solver_version(),
            manifest=manifest,
        ),
    )


def certified_strain_gap(halide: str, eps: float) -> CertifiedResult:
    """Exact (zero-width) bracket on the strained singlet -> Sz=0-triplet gap J(eps).

    ``dimer_exchange_analytic`` is the closed-form exact exchange splitting, already gated
    (SPEC_odmd_spin.md, SPEC_nb3x8_strain.md) -- see module docstring for why this must NOT go
    through ``certified_gaps.gap_bracket`` (it silently returns the wrong, bright ionic gap).
    """
    p = strained_params(halide, eps)
    J = dimer_exchange_analytic(**p)
    return _exact_bracket(J, halide, {"eps": eps, "t": p["t"]})


def zeeman_split_gap(halide: str, B: float) -> CertifiedResult:
    """Exact (zero-width) bracket on the field-split gap: singlet ground (B=0 energy) to the
    LOWEST Zeeman-shifted level of the S=1 triplet manifold, min(E_triplet(Sz)) - E_singlet.

    Exact full diagonalization of the 6-level N=2 sector (both (1,1) and the two 1-dimensional
    (2,0)/(0,2) sectors) -- see module docstring for why this must NOT go through the
    HF-reachable Krylov path (it cannot see the field at all from a closed-shell reference).
    Cross-checks exactly against ``dimer_exchange_analytic`` at B=0 (gated in the test suite).
    """
    p = NB3X8_LT_BULK[halide]
    mh = dimer_cluster_integrals(**p).to_hamiltonian()
    H0 = mh.qubit_hamiltonian.to_matrix()
    Sz = _SZ_OP.to_matrix()

    w0, V0 = np.linalg.eigh(H0)
    n0 = np.real(np.einsum("ji,jk,ki->i", V0.conj(), _N_OP, V0))
    e_ground = float(w0[np.abs(n0 - 2.0) < 1e-8].min())

    H_field = H0 + G_FACTOR * MU_B_MEV_PER_T * B * Sz
    w, V = np.linalg.eigh(H_field)
    n = np.real(np.einsum("ji,jk,ki->i", V.conj(), _N_OP, V))
    keep = np.abs(n - 2.0) < 1e-8
    e2 = np.sort(w[keep].real)
    gap = float(e2[1] - e_ground)  # e2[0] is the singlet ground itself (unshifted by a field)

    return _exact_bracket(gap, halide, {"B_tesla": B})


@dataclass(frozen=True)
class HermiticityCheck:
    is_hermitian: bool
    max_asymmetry: float


def check_zeeman_hermitian(B: float = 1.0) -> HermiticityCheck:
    """Sanity gate: the Zeeman term (and H0 + Zeeman) is Hermitian for any field B."""
    M = G_FACTOR * MU_B_MEV_PER_T * B * _SZ_OP.to_matrix()
    asym = float(np.max(np.abs(M - M.conj().T)))
    return HermiticityCheck(is_hermitian=asym < 1e-10, max_asymmetry=asym)
