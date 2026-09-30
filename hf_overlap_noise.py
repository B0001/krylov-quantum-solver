#!/usr/bin/env python3
"""
The overlap certificate's noise sibling -- and its noise is structurally different.

`hf_overlap_certificate.certify_hf_overlap` certifies gamma_min <= |<HF|psi_0>| from a certified
E_1 floor beta, via Davis-Kahan: gamma_min = sqrt(1 - r^2/delta^2), delta = beta - lambda_u. Unlike
every other certificate in the noise arc (`certified_noise`, `gap_selfcheck_noise`,
`certified_thermochem_noise`), for a DETERMINANT guiding state u = HF, lambda_u = <HF|H|HF> and
r = ||(H - lambda_u)HF|| are exact classical computations on a known computational-basis state --
zero shot noise. ALL the noise in this certificate enters through the single gap input beta (the
self-mode Weinstein floor theta_1 - sigma_1, which DOES need <H>, <H^2> measurements on the Krylov
first-excited Ritz state).

THE HYPOTHESIS (specs/SPEC_hf_overlap_noise.md): gamma_min = sqrt(1 - r^2/delta^2) has an UNBOUNDED
derivative in delta as delta -> r. So unlike the other certificates in this arc -- whose coverage
degrades roughly N-independently but *smoothly* across systems -- this one should show a sharper
break: near-perfect coverage at wide margin (delta_exact - r large), collapsing fast as the margin
thins, because a small downward fluctuation in beta that would cost little at wide margin instead
pushes r >= delta and votes the certificate vacuous (or, on the other side, a small upward
fluctuation of beta lets gamma_min_noisy exceed the true overlap -- an invalid, over-claiming
certificate). The two failure modes are tracked separately (vacuous vs invalid) because they are
not equally bad: vacuous is a conservative "I don't know", invalid is a silent overclaim.

NOISE MODEL (reuses `certified_noise`'s i.i.d.-Gaussian, lambda-1-norm machinery and
`gap_selfcheck_noise._pad`'s post-hoc padding convention -- see module docstrings there):
  * lambda_u, r: EXACT (no noise) -- the determinant-only shortcut this module exists to test.
  * beta = theta_1 - sigma_1 (self mode): theta_1 and <H^2>_1 get independent noisy realizations at
    N shots (se ~ lambda_H/sqrt(N), lambda_{H^2}/sqrt(N)), exactly `certified_gaps.gap_bracket`'s
    formula rebuilt from the noisy pair -- then padded ONE-SIDED by z*lambda_H/sqrt(N) (subtracted,
    never added: beta is a lower bound and the pad must only make it more conservative, mirroring
    `gap_selfcheck_noise._pad`'s "-half_width" on gap_lower).

SECTOR/SCOPE BOUNDARY (state this up front, not discovered later): the shot-free-residual argument
is SPECIFIC to determinant guiding states. A chained-Ritz state v (`krylov_refine.py`,
specs/SPEC_chained_overlap.md) is a Krylov Ritz vector, not a computational basis state -- its own
lambda_v, r_v need <v|H|v> and <v|H^2|v>, which DO cost shots. Chaining puts r back on the noise
budget; nothing here generalizes to it.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from certified_gaps import gap_bracket
from certified_noise import certified_half_width, hamiltonian_one_norms
from hf_overlap_certificate import exact_reachable_overlap
from hybrid_quantum_solver.certified_overlap import rayleigh_quotient, residual_norm
from hybrid_quantum_solver.molecular_hamiltonian import MolecularHamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from temple_bounds import mean_and_variance

_VALIDITY_TOL = 1e-9


def exact_margin(mh: MolecularHamiltonian, m: int,
                 solver: Optional[QuantumKrylovSolver] = None) -> dict:
    """Zero-noise reference point for one system: lambda_u, r (exact, classical), the self-mode
    delta = eps1 - lambda_u, the margin delta - r, and the exact reachable overlap. This is the
    x-axis of the coverage-vs-margin curve; it never touches the noise model below."""
    solver = solver if solver is not None else QuantumKrylovSolver(mh)
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    u = np.asarray(mh.hf_state().data, dtype=complex)
    lam = rayleigh_quotient(H, u)
    r = residual_norm(H, u, lam)
    gb = gap_bracket(mh, m, solver=solver)
    delta = gb.eps1 - lam
    exact = exact_reachable_overlap(mh)
    return dict(lam=lam, r=r, delta=delta, margin=delta - r,
               r_over_delta=(r / delta if delta > 0 else np.inf), exact_overlap=exact)


def hf_overlap_noise_coverage(mh: MolecularHamiltonian, m: int, shots: float, z: float = 2.0,
                              trials: int = 6000, seed: int = 0,
                              solver: Optional[QuantumKrylovSolver] = None) -> dict:
    """Monte-Carlo coverage of the exact reachable overlap by the noisy HF overlap certificate.

    Only beta = theta_1 - sigma_1 is noisy (N shots on the first-excited Ritz state); lambda_u and
    r (u = HF, a determinant) are exact in every trial -- the structural difference this module
    measures. ``z`` pads beta one-sided (more conservative), the `gap_selfcheck_noise._pad`
    convention applied to a single lower-bound input instead of a two-sided bracket.

    Returns per-trial fates partitioned into three disjoint buckets:
      * ``coverage``    -- non-vacuous AND valid (gamma_min_noisy <= exact overlap): the useful,
                           correct case the certificate is FOR.
      * ``frac_vacuous``-- r >= delta_noisy: a conservative "can't certify", never wrong.
      * ``frac_invalid``-- non-vacuous but gamma_min_noisy > exact overlap: a silent overclaim,
                           the failure mode a certificate must never produce.
    plus the exact (zero-noise) margin diagnostics from `exact_margin` for the coverage-vs-margin
    curve.
    """
    solver = solver if solver is not None else QuantumKrylovSolver(mh)
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    u = np.asarray(mh.hf_state().data, dtype=complex)
    lam = rayleigh_quotient(H, u)
    r = residual_norm(H, u, lam)
    exact_overlap = exact_reachable_overlap(mh)

    _, states = solver.eigenstates(m, n_states=2)
    if len(states) < 2:
        raise ValueError(f"Krylov subspace at M={m} is rank-1: no first-excited Ritz state to "
                         "floor E1 from. Increase m.")
    th1, var1 = mean_and_variance(H, states[1])
    h2_1 = var1 + th1 * th1
    lam_h, lam_h2 = hamiltonian_one_norms(mh)
    se_h, se_h2 = lam_h / np.sqrt(shots), lam_h2 / np.sqrt(shots)
    half_width = certified_half_width(lam_h, shots, z)

    rng = np.random.default_rng(seed)
    n_th1 = th1 + rng.normal(0.0, se_h, trials)
    n_h2_1 = h2_1 + rng.normal(0.0, se_h2, trials)
    n_var1 = np.clip(n_h2_1 - n_th1 * n_th1, 0.0, None)
    eps1_noisy = n_th1 - np.sqrt(n_var1)
    beta = eps1_noisy - half_width                      # one-sided pad: more conservative only
    delta = beta - lam                                  # lam is exact -- the only noisy term above

    vacuous = delta <= r                                # covers delta <= 0 too (r >= 0 always)
    safe_delta = np.where(vacuous, 1.0, delta)           # avoid div-by-zero on vacuous entries
    ratio = r / safe_delta
    gamma = np.where(vacuous, 0.0, np.sqrt(np.clip(1.0 - ratio * ratio, 0.0, 1.0)))

    invalid = (~vacuous) & (gamma > exact_overlap + _VALIDITY_TOL)
    covered = (~vacuous) & ~invalid

    out = dict(
        coverage=float(np.mean(covered)),
        frac_vacuous=float(np.mean(vacuous)),
        frac_invalid=float(np.mean(invalid)),
        lam_h=lam_h, lam_h2=lam_h2, half_width=half_width,
    )
    out.update(exact_margin(mh, m, solver=solver))
    return out


if __name__ == "__main__":
    from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian

    cases = {
        "H2 eq (0.74A, wide)": dict(atom="H 0 0 0; H 0 0 0.74"),
        "H4 chain (0.9/1.8/2.7A)": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
        "H4 chain (1.0/2.0/3.0A, thin)": dict(atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0"),
        "square H4 a=1.2A (vacuous)": dict(atom="H 0 0 0; H 1.2 0 0; H 1.2 1.2 0; H 0 1.2 0"),
    }
    M = 8
    print("=" * 100)
    print("HF overlap certificate under shot noise -- noise enters ONLY through the gap floor beta")
    print("  system                          | margin  | shots  | z | coverage | vacuous | invalid")
    for name, spec in cases.items():
        mh = build_molecular_hamiltonian(**spec)
        solver = QuantumKrylovSolver(mh)
        for shots in (1e4, 1e5, 1e6):
            for z in (0, 1, 2, 3):
                r = hf_overlap_noise_coverage(mh, M, shots, z=float(z), solver=solver)
                print(f"  {name:32s}| {r['margin']:7.4f} | {shots:.0e} | {z} | "
                      f"{r['coverage']:.3f}    | {r['frac_vacuous']:.3f}   | {r['frac_invalid']:.3f}")
    print("=" * 100)
