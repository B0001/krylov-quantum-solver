#!/usr/bin/env python3
"""
Shift both sides -- the BLISS/SCDF symmetry shift helps the NEAR-TERM arc too, moving the crossover.

`precision_cost` / `cost_advisor` compared the near-term arc using the RAW measurement 1-norm
lambda_meas against fault-tolerant qubitization using the SHIFTED lambda_DF. That is an
inconsistency: the number-operator shift H -> H + (b1 N + b2)(N - n_e) is a spectrum-preserving
constant on the target sector (`df_factorization.symmetry_shift`), so its effect on ANY 1-norm is
fair game for EITHER method. Shift both sides and the FT advantage shrinks.

THE HONEST METRIC (and the finding). The shot count is N = (z lambda / eps)^2 where lambda must sum
only the NON-IDENTITY Pauli coefficients: the identity term is a constant with ZERO variance, so it
costs zero shots. The repo's `precision_cost.measurement_lambda` summed ALL coefficients, identity
included -- overstating near-term cost (fixed since specs/SPEC_lambda_meas_identity.md; this module
found it). That matters here because a large part of what the shift
does to the qubit Hamiltonian is dump weight INTO the identity term (N2: identity 8.55 -> 0.10 Ha).
Scoring the shift with the identity-inclusive 1-norm therefore FLATTERS it:

    reduction in lambda_meas      H2      H2O(4,3)   N2(6,6)
      identity-INCLUDED         53.9%      51.4%      73.2%   <- overstated
      identity-EXCLUDED (shot)  42.9%      39.4%      57.9%   <- the honest number

Both tell the same qualitative story -- the shift really does cut the near-term 1-norm -- but the
crossover movement is 2.7-5.7x, not the 4-14x the inclusive metric advertises. `include_identity`
is exposed so both can be computed; it defaults to False (the honest, shot-relevant one).

WHAT WE CLAIM: the shift lowers lambda_meas materially (>= 35%); the certified shot cost drops by
(lambda_raw/lambda_shift)^2 (2.7-5.7x); the fair both-sided flip-rho is that same factor below the
one-sided bridge value -- so the one-sided bridge OVERSTATES FT's advantage. WHAT WE DO NOT CLAIM:
that the shift changes the EXPONENTS (still 1/eps^2 vs 1/eps -- it moves constants only), nor that
it is free on hardware (it adds number-operator terms that must themselves be measured).

See specs/SPEC_shift_both_sides.md. Reference: spectrum preservation is exact and is checked
against FCI (G4); the flip-rho ratio is the algebra (lambda_raw/lambda_shift)^2.
"""
from __future__ import annotations

import numpy as np

from df_factorization import symmetry_shift
from hybrid_quantum_solver.molecular_hamiltonian import build_hamiltonian_from_integrals
from precision_cost import ft_queries, near_term_shots, qubitization_lambda


def shot_lambda(mh, *, include_identity: bool = False) -> float:
    """The measurement 1-norm that actually costs shots: sum |c_P| over NON-identity Paulis P.

    The identity term is a constant of zero variance -- it consumes no shots. Set
    ``include_identity=True`` to reproduce the identity-inclusive 1-norm (what
    `precision_cost.measurement_lambda` computed before SPEC_lambda_meas_identity fixed it).
    """
    op = mh.qubit_hamiltonian
    total = 0.0
    for coeff, pauli in zip(op.coeffs, op.paulis):
        label = pauli.to_label()
        if not include_identity and set(label) == {"I"}:
            continue
        total += abs(complex(coeff))
    return float(total)


def measurement_lambda_from_integrals(h1, eri, nelec, e_core: float = 0.0, *,
                                      include_identity: bool = False) -> float:
    """lambda_meas of the qubit Hamiltonian built from these (possibly shifted) integrals."""
    mh = build_hamiltonian_from_integrals(np.asarray(h1), np.asarray(eri), nelec, float(e_core))
    return shot_lambda(mh, include_identity=include_identity)


def optimize_meas_shift(h1, eri, norb: int, nelec, *, maxiter: int = 60):
    """Re-optimize (b1, b2) to minimise lambda_meas DIRECTLY (rather than lambda_DF).

    Seeded from the lambda_DF-optimal (SCDF) shift, so it can only match or improve on it.
    Returns ``(b1, b2)``.
    """
    from scipy.optimize import minimize

    _, _, _, (b1_df, b2_df) = symmetry_shift(h1, eri, norb, nelec)

    def objective(params):
        h1_s, eri_s, e_shift, _ = symmetry_shift(h1, eri, norb, nelec, b1=params[0], b2=params[1])
        return measurement_lambda_from_integrals(h1_s, eri_s, nelec, e_shift)

    res = minimize(objective, x0=[b1_df, b2_df], method="Nelder-Mead",
                   options={"xatol": 1e-5, "fatol": 1e-7, "maxiter": maxiter})
    return float(res.x[0]), float(res.x[1])


def shifted_measurement_lambda(h1, eri, norb: int, nelec, *, target: str = "df",
                               include_identity: bool = False) -> float:
    """lambda_meas after the symmetry shift.

    ``target="df"``   -- the SCDF shift (b1, b2 optimised for lambda_DF), i.e. the very same shift
                         the FT side already takes credit for. This is the fair comparison.
    ``target="meas"`` -- (b1, b2) re-optimised for lambda_meas itself (an upper bound on what the
                         shift can do for the near-term arc).
    ``target="raw"``  -- no shift (the baseline).
    """
    if target == "raw":
        return measurement_lambda_from_integrals(h1, eri, nelec, 0.0,
                                                 include_identity=include_identity)
    if target == "df":
        h1_s, eri_s, e_shift, _ = symmetry_shift(h1, eri, norb, nelec)
    elif target == "meas":
        b1, b2 = optimize_meas_shift(h1, eri, norb, nelec)
        h1_s, eri_s, e_shift, _ = symmetry_shift(h1, eri, norb, nelec, b1=b1, b2=b2)
    else:
        raise ValueError(f"target must be 'raw', 'df' or 'meas' (got {target!r})")
    return measurement_lambda_from_integrals(h1_s, eri_s, nelec, e_shift,
                                             include_identity=include_identity)


def flip_rho(lam_meas: float, lam_df: float, eps: float, z: float = 2.0) -> float:
    """The per-query/per-shot cost ratio rho at which the FT-vs-near-term verdict flips: N/Q.

    Same quantity as `cost_advisor.Verdict.rho_flip`. FT is cheaper only if the true rho is BELOW
    this; a smaller flip-rho is a WEAKER case for FT.
    """
    return near_term_shots(lam_meas, eps, z) / ft_queries(lam_df, eps)


def fair_flip_rho(h1, eri, norb: int, nelec, eps: float, z: float = 2.0, *,
                  include_identity: bool = False):
    """(one_sided, both_sided) flip-rho at precision ``eps``.

    ``one_sided``  -- the bridge as it stood: RAW lambda_meas vs SHIFTED lambda_DF (inconsistent).
    ``both_sided`` -- the fair comparison: SHIFTED lambda_meas vs SHIFTED lambda_DF.

    both_sided == one_sided / (lam_raw/lam_shift)^2 -- the shot cost is quadratic in the 1-norm, so
    the same shift that buys the FT side its constant buys the near-term side a bigger one.
    """
    lam_df_shift = qubitization_lambda(h1, eri, norb, nelec=nelec, shift=True)
    lam_meas_raw = shifted_measurement_lambda(h1, eri, norb, nelec, target="raw",
                                              include_identity=include_identity)
    lam_meas_shift = shifted_measurement_lambda(h1, eri, norb, nelec, target="df",
                                                include_identity=include_identity)
    return (flip_rho(lam_meas_raw, lam_df_shift, eps, z),
            flip_rho(lam_meas_shift, lam_df_shift, eps, z))


def spectrum_preserved(h1, eri, norb: int, nelec, tol: float = 1e-8) -> bool:
    """FCI(shifted) + e_shift == FCI(raw)? The shift is exact, not an approximation (G4)."""
    from pyscf import fci

    h1_s, eri_s, e_shift, _ = symmetry_shift(h1, eri, norb, nelec)
    e_raw = fci.direct_spin1.kernel(np.asarray(h1), np.asarray(eri), norb, nelec)[0]
    e_shifted = fci.direct_spin1.kernel(np.asarray(h1_s), np.asarray(eri_s), norb, nelec)[0]
    return bool(abs((e_shifted + e_shift) - e_raw) < tol)


if __name__ == "__main__":
    from pyscf import ao2mo, gto, mcscf, scf

    def cas(atom, norb, nelec):
        mf = scf.RHF(gto.M(atom=atom, basis="sto-3g", verbose=0)).run()
        c = mcscf.CASCI(mf, norb, nelec)
        c.kernel()
        h1, ec = c.get_h1eff()
        na = (nelec + nelec % 2) // 2
        return (np.asarray(h1), ao2mo.restore(1, c.get_h2eff(), norb), float(ec),
                (na, nelec - na), norb)

    mols = {"H2": ("H 0 0 0; H 0 0 0.74", 2, 2),
            "H2O(4,3)": ("O 0 0 0.117; H 0 0.757 -0.467; H 0 -0.757 -0.467", 3, 4),
            "N2(6,6)": ("N 0 0 1.10", 6, 6)}
    mols["N2(6,6)"] = ("N 0 0 0; N 0 0 1.10", 6, 6)

    eps = 1.6e-3
    print("Shift BOTH sides: the symmetry shift cuts the near-term 1-norm too, moving the crossover")
    print(f"(eps = {eps:.1e} Ha, z = 2; lambda_meas EXCLUDES the identity term -- it costs no shots)\n")
    print(f"{'mol':9s} | {'lam_raw':>7} | {'lam_shift':>9} | {'red%':>5} | {'shot gain':>9} | "
          f"{'flip-rho 1-sided':>16} | {'flip-rho fair':>13} | {'overstated by':>13}")
    for name, (atom, norb, nel) in mols.items():
        h1, eri, ec, ne, no = cas(atom, norb, nel)
        lam_raw = shifted_measurement_lambda(h1, eri, no, ne, target="raw")
        lam_shift = shifted_measurement_lambda(h1, eri, no, ne, target="df")
        one, both = fair_flip_rho(h1, eri, no, ne, eps)
        gain = (lam_raw / lam_shift) ** 2
        print(f"{name:9s} | {lam_raw:7.3f} | {lam_shift:9.3f} | {100*(1-lam_shift/lam_raw):5.1f} | "
              f"{gain:8.2f}x | {one:16.2e} | {both:13.2e} | {one/both:12.2f}x")
    print("\nThe one-sided bridge (raw lambda_meas vs shifted lambda_DF) OVERSTATES FT's advantage by")
    print("the shot gain. Honest caveat: with the identity term wrongly included in lambda_meas the")
    print("gain looks like 4-14x instead of 2.7-5.7x -- see the module docstring.")
