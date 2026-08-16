#!/usr/bin/env python3
"""Emit certkit certificates from real QKSD Krylov energy bounds (certkit-jn1.1, jn1.3).

The solver produces (eigenvector estimate, gap parameter, bracket). This wraps them in
the certificate format and hands them to certkit's checker, which re-derives the bracket
in interval arithmetic without importing anything from here.

Nothing in this file is trusted: it is a producer. The only question it answers is
whether the solver's bound survives an independent re-derivation.

Run:  PYTHONPATH=../certkit python certkit_bridge.py [H2|H4|H4-stretched|N2]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from certkit.checker import check
from certkit.interval import Iv
from certkit.operators import DENSE_LIMIT, decode_operator, encode_pauli, operator_ref
from certkit.producer import _pad          # producer-side padding convention
from certkit.schema import f2h, seal

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


def encode_hamiltonian(mh) -> dict:
    """The qubit Hamiltonian as a certkit Pauli-sum operator.

    NOT as CSR. Materialising the matrix sums each Pauli term in floating point, and
    the (i,j) and (j,i) entries accumulate in different orders -- H4 comes out
    asymmetric by 8.3e-17, which certkit refuses (correctly: it is a different
    operator than the symmetric one you meant). The Pauli encoding carries the terms
    themselves, so symmetry holds by construction rather than by luck.

    certkit's Pauli strings are LSB-first; qiskit's labels are MSB-first.
    """
    op = mh.qubit_hamiltonian
    terms = []
    for label, coeff in op.to_list():
        c = complex(coeff)
        if c.imag != 0.0:
            raise ValueError(f"complex coefficient on {label}: {c}")
        if label.count("Y") % 2:
            raise ValueError(f"odd number of Y factors in {label} -- operator is not real")
        terms.append((c.real, label[::-1]))
    return encode_pauli(op.num_qubits, terms)


def krylov_witness(H, mh, m: int):
    """The solver's own numbers: trial vector and its mean/variance.

    Electronic frame -- the certificate is about the matrix, and mh.energy_offset is
    a trusted addition that happens outside it.
    """
    solver = QuantumKrylovSolver(mh)
    energies, states = solver.eigenstates(m, n_states=2)
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


def full_space_beta(H):
    """A separator for the WHOLE matrix, not the reachable Krylov sector.

    Untrusted like every producer number: the checker discharges beta by inertia
    count, so a wrong guess costs coverage and never soundness.
    """
    from scipy.sparse.linalg import eigsh

    vals = np.sort(eigsh(H.real.astype(float), k=2, which="SA", return_eigenvectors=False))
    return 0.5 * (vals[0] + vals[1]), vals


def gershgorin_lower(enc: dict) -> float:
    """A floor on every eigenvalue, read straight off the rows. No gap needed."""
    op = decode_operator(enc)
    lower = float("inf")
    for i in range(op.n):
        entries = op.row(i)
        diag = entries.get(i, Iv.exact(0.0)).lo
        radius = sum(v.mag_ub for j, v in entries.items() if j != i)
        lower = min(lower, diag - radius)
    return lower


def build_certificate(enc: dict, rule: str, x, lo: float, hi: float,
                      beta: float | None = None) -> dict:
    """The solver's claim in the certificate format. No certkit producer involved."""
    witness = {"rule": rule, "vector": [f2h(float(v)) for v in x]}
    if beta is not None:
        witness["beta"] = f2h(float(beta))
    cert = {
        "schema": "certkit/1",
        "claim": {
            "operator_ref": operator_ref(enc),
            "kind": "lambda_min_enclosure",
            "enclosure": {"lo": f2h(float(lo)), "hi": f2h(float(hi))},
        },
        "witness": witness,
        "producer": {"name": "hybrid_quantum_solver.QuantumKrylovSolver", "backend": "numpy"},
    }
    return seal(cert)


def report(tag: str, enc: dict, cert: dict, fci: float, out: Path, name: str) -> bool:
    (out / f"{name}.json").write_text(json.dumps(cert, indent=2))
    v = check(cert, enc)
    print(f"  [{tag:19s}] {'VERIFIED' if v.ok else 'ABSTAIN '}  {v.reason}")
    if v.ok:
        lo, hi = v.enclosure
        print(f"  {'':19s}   width {hi - lo:.3e}   exact inside: {lo <= fci <= hi}")
    return v.ok


def main() -> int:
    case = sys.argv[1] if len(sys.argv) > 1 else "H2"
    mh = build_molecular_hamiltonian(**CASES[case])
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    fci = mh.ground_state_energy() - mh.energy_offset     # comparison only

    enc = encode_hamiltonian(mh)
    x, th0, var0, sector_eps, info = krylov_witness(H, mh, KRYLOV_DIM)
    n = 1 << mh.qubit_hamiltonian.num_qubits

    print(f"=== {case} ===")
    print(f"  dim {n}   pauli terms {len(enc['terms'])}   krylov_dim {KRYLOV_DIM}")
    print(f"  discarded imaginary weight  {info['imag_weight']:.3e}")
    print(f"  witness variance            {var0:.3e}")
    print(f"  exact lambda_min            {fci!r}")

    out = Path("certkit_out") / case
    out.mkdir(parents=True, exist_ok=True)
    (out / "operator.json").write_text(json.dumps(enc, indent=2))

    ok = False
    if n <= DENSE_LIMIT:
        beta, vals = full_space_beta(H)
        print(f"  sector eps {sector_eps:.6f} vs matrix lambda_2 {vals[1]:.6f}"
              f"  -> premise {'holds' if sector_eps <= vals[1] else 'FAILS'}")
        # (1) the solver's bound transcribed as-is, with the sector's own separator
        lo = th0 - var0 / (sector_eps - th0) if sector_eps > th0 else -np.inf
        report("temple/sector beta", enc,
               build_certificate(enc, "temple_inertia", x, lo, th0, sector_eps),
               fci, out, "certificate_sector")
        # (2) same witness, separator taken over the whole matrix
        lo = th0 - var0 / (beta - th0)
        pad = _pad(th0, 1e-9, len(x), th0 - lo)
        ok = report("temple/matrix beta", enc,
                    build_certificate(enc, "temple_inertia", x, lo - pad, th0 + pad, beta),
                    fci, out, "certificate")
    else:
        print(f"  n > DENSE_LIMIT ({DENSE_LIMIT}): no inertia route, so no gap discharge")

    # The loose route needs no gap and no factorisation, and always applies.
    lo = gershgorin_lower(enc)
    pad = _pad(th0, 1e-9, len(x), th0 - lo)
    name = "certificate" if not ok else "certificate_gershgorin"
    ok |= report("gershgorin_rayleigh", enc,
                 build_certificate(enc, "gershgorin_rayleigh", x, lo - pad, th0 + pad),
                 fci, out, name)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
