#!/usr/bin/env python3
"""Emit certkit certificates from real QKSD Krylov energy bounds (certkit-jn1.1, jn1.3).

The solver produces (eigenvector estimate, gap parameter, bracket). This wraps them in
the certificate format and writes them to disk, then runs certkit's checker as a separate
process over those files. Per certkit INTEGRATION.md the checker is consumed as a protocol,
not imported: nothing here may call check() in-process.

Nothing in this file is trusted: it is a producer. The only question it answers is
whether the solver's bound survives an independent re-derivation.

Run:  uv run python certkit_bridge.py [H2|H4|H4-stretched|N2]

Since portfolio-0pc the producer itself lives in ``hybrid_quantum_solver.certification``,
because every energy the solver returns is certified through it; the names this file used to
define are re-exported below. What stays here is the fixed regression cases and their CLI.

The gate over what this emits is tests/test_certkit_regression_gate_spec.py, which pins, per
case, the solver's own return-path certificate too (run_case's solve/ entry). Every return
path certifies through certification.certify_energies; all of them are gated on H2 by
tests/test_certkit_return_paths_spec.py and tests/test_certificate_or_abstention_spec.py.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from certkit.operators import DENSE_LIMIT
from certkit.producer import pad_claim     # public producer-side padding convention

from hybrid_quantum_solver.certification import (  # noqa: F401  -- re-exported producer names
    Verdict,
    build_certificate,
    certify_ground_state,
    check_certificate,
    emit,
    encode_hamiltonian,
    full_space_beta,
    gershgorin_lower,
    real_witness,
)
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from hybrid_quantum_solver.quantum_krylov_solver import QuantumKrylovSolver
from temple_bounds import mean_and_variance

CASES = {
    "H2": dict(atom="H 0 0 0; H 0 0 0.74"),
    "H4": dict(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7"),
    "H4-stretched": dict(atom="H 0 0 0; H 0 0 2.0; H 0 0 4.0; H 0 0 6.0"),
    "N2": dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6),
}
KRYLOV_DIM = 8


def krylov_witness(H, mh, m: int):
    """The solver's own numbers: trial vector and its mean/variance.

    Electronic frame -- the certificate is about the matrix, and mh.energy_offset is
    a trusted addition that happens outside it.
    """
    solver = QuantumKrylovSolver(mh)
    # The Ritz decomposition only: the energies are diagnostics here, and the certificates
    # run_case builds from these states are checked one by one below.
    spectrum, states = solver.uncertified_eigenstates(m, n_states=2)
    energies = spectrum.unverified_energies
    psi0 = states[0]

    # The real-time Krylov basis is complex; H is real symmetric, so its ground
    # eigenvector can be chosen real. Take the dominant real direction.
    weight = float(psi0.imag @ psi0.imag)
    x = psi0.real if weight < 0.5 else psi0.imag
    x = x / np.linalg.norm(x)

    # Mean and variance must describe the vector the certificate carries, not psi0:
    # the Rayleigh quotient of the real part is a different, higher number.
    th0, var0 = mean_and_variance(H, x)
    th1, var1 = mean_and_variance(H, states[1])
    sector_eps = th1 - np.sqrt(var1)                # self mode: unverifiable premise
    return x, th0, var0, sector_eps, dict(imag_weight=weight, energies=energies)


def run_case(case: str, out: Path, conv_tol: float = 1e-9) -> tuple[float, list[Verdict]]:
    """Emit every certificate this case supports, and check each one independently.

    Returns the exact electronic-frame lambda_min -- for comparison only, no certificate
    depends on it -- and one Verdict per certificate actually emitted.
    """
    mh = build_molecular_hamiltonian(**CASES[case], conv_tol=conv_tol)
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    lam = mh.ground_state_energy() - mh.energy_offset     # comparison only

    enc = encode_hamiltonian(mh)
    x, th0, var0, sector_eps, info = krylov_witness(H, mh, KRYLOV_DIM)
    n = 1 << mh.qubit_hamiltonian.num_qubits

    print(f"=== {case} ===")
    print(f"  dim {n}   pauli terms {len(enc['terms'])}   krylov_dim {KRYLOV_DIM}")
    print(f"  discarded imaginary weight  {info['imag_weight']:.3e}")
    print(f"  witness variance            {var0:.3e}")
    print(f"  exact lambda_min            {lam!r}")

    out.mkdir(parents=True, exist_ok=True)
    (out / "operator.json").write_text(json.dumps(enc, indent=2))

    def emit_and_report(*args, **kwargs) -> Verdict:
        """emit, plus the checker's line on stdout: this is the report, not the library."""
        v = emit(*args, **kwargs)
        print(f"  [{v.name:23s}] {v.line or '(no checker output)'}")
        return v

    verdicts: list[Verdict] = []
    if n <= DENSE_LIMIT:
        beta, vals = full_space_beta(H)
        print(f"  sector eps {sector_eps:.6f} vs matrix lambda_2 {vals[1]:.6f}"
              f"  -> premise {'holds' if sector_eps <= vals[1] else 'FAILS'}")
        # (1) the solver's bound transcribed as-is, with the sector's own separator.
        # Temple says nothing finite unless eps > theta -- a different comparison from the
        # premise printed above, and the one that decides whether a certificate exists at
        # all. When it fails, the honest output is no certificate: -inf is "valid but
        # vacuous" in temple_bounds' vocabulary, and certkit refuses non-finite floats.
        if sector_eps > th0:
            lo = th0 - var0 / (sector_eps - th0)
            verdicts.append(emit_and_report(enc, "temple_inertia", x, lo, th0, out,
                                            "certificate_sector", beta=sector_eps))
        else:
            print(f"  [{'certificate_sector':23s}] not emitted: sector eps {sector_eps:.6f}"
                  f" <= theta {th0:.6f}, so Temple gives no finite bound")
        # (2) same witness, separator taken over the whole matrix
        lo = th0 - var0 / (beta - th0)
        pad = pad_claim(th0, 1e-9, len(x), th0 - lo)
        verdicts.append(emit_and_report(enc, "temple_inertia", x, lo - pad, th0 + pad, out,
                                        "certificate_temple", beta=beta))
    else:
        print(f"  n > DENSE_LIMIT ({DENSE_LIMIT}): no inertia route, so no gap discharge")

    # The loose route needs no gap and no factorisation, and always applies.
    lo = gershgorin_lower(enc)
    pad = pad_claim(th0, 1e-9, len(x), th0 - lo)
    verdicts.append(emit_and_report(enc, "gershgorin_rayleigh", x, lo - pad, th0 + pad, out,
                                    "certificate_gershgorin"))

    # What the solver's own return path emits for the energy it returns (portfolio-cir; on
    # every call since portfolio-0pc): certify_ground_state, whose ceiling is the returned
    # number rather than a witness's Rayleigh quotient. One verdict, from the route that won.
    # It writes its own operator.json and reuses the route filenames, hence solve/.
    v = QuantumKrylovSolver(mh).solve(KRYLOV_DIM, certificate_dir=out / "solve").certificate
    print(f"  [{'solve/' + v.name:23s}] {v.line or '(no checker output)'}")
    verdicts.append(v._replace(name=f"solve/{v.name}"))

    for v in verdicts:
        if v.ok:
            print(f"  {v.name:23s}   width {v.hi - v.lo:.3e}"
                  f"   exact inside: {v.lo <= lam <= v.hi}")
    return lam, verdicts


def main() -> int:
    case = sys.argv[1] if len(sys.argv) > 1 else "H2"
    _, verdicts = run_case(case, Path("certkit_out") / case)
    # "Some route verified" is weaker than "this energy is certified". The Gershgorin route
    # needs no gap, always applies, and certifies [gershgorin_lower(H), <x|H|x>] for ANY unit
    # vector x -- noise included. What the certificates are actually worth is decided per
    # route and against a reference, by tests/test_certkit_regression_gate_spec.py, not here.
    #
    # This disjunction is now over every emitted certificate, including the sector one, whose
    # verdict the pre-refactor code dropped on the floor. No observable change -- the sector
    # route has never verified -- but it is a wider disjunction than what it replaced, so say
    # so rather than let a reader assume the old exclusion still holds. It also covers the
    # solver's own return-path certificate (solve/), since chem-cgb.
    return 0 if any(v.ok for v in verdicts) else 1


if __name__ == "__main__":
    sys.exit(main())
