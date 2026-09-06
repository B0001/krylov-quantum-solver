#!/usr/bin/env python3
"""Emit certkit certificates from real QKSD Krylov energy bounds (certkit-jn1.1, jn1.3).

The solver produces (eigenvector estimate, gap parameter, bracket). This wraps them in
the certificate format and writes them to disk, then runs certkit's checker as a separate
process over those files. Per certkit INTEGRATION.md the checker is consumed as a protocol,
not imported: nothing here may call check() in-process.

Nothing in this file is trusted: it is a producer. The only question it answers is
whether the solver's bound survives an independent re-derivation.

Run:  uv pip install -e ".[certkit]"  then  uv run python certkit_bridge.py [H2|H4|H4-stretched|N2]

The gate over what this emits is tests/test_certkit_regression_gate_spec.py.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import NamedTuple

import numpy as np

from certkit.interval import Iv
from certkit.operators import DENSE_LIMIT, decode_operator, encode_pauli, operator_ref
from certkit.producer import pad_claim     # public producer-side padding convention
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


class Verdict(NamedTuple):
    """One certificate, and what the independent checker said about it.

    ``line`` is the checker's first stdout line, verbatim, or a ``CRASH``-prefixed
    diagnostic if it printed no verdict at all. Both an ABSTAIN and a crash exit 1, so
    the line is the only thing that tells them apart.
    """

    name: str
    rule: str
    ok: bool
    line: str
    lo: float
    hi: float


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


def check_certificate(cert_path: Path, operator_path: Path) -> tuple[bool, str]:
    """Run the checker as a separate process over the files, and report what it said.

    The verdict is the checker's exit status. This file is a producer and is
    untrusted; it must not adjudicate its own claims in its own process.

    Exit 1 covers both ABSTAIN and a crash, and confusing the two would let a broken
    checker pass for an expected abstention. The checker prints its verdict to stdout
    and nothing else there, so an empty stdout means it died: that is returned with a
    CRASH prefix, never as a bare traceback line, which a caller matching on the text
    would otherwise read as a verdict.
    """
    proc = subprocess.run(
        [sys.executable, "-m", "certkit.cli", "check", str(cert_path), str(operator_path)],
        capture_output=True, text=True,
    )
    verdict = proc.stdout.strip().splitlines()
    if verdict:
        return proc.returncode == 0, verdict[0]
    err = proc.stderr.strip().splitlines()
    return False, f"CRASH  {err[-1] if err else '(no output at all)'}"


def emit(enc: dict, rule: str, x, lo: float, hi: float, out: Path, name: str,
         beta: float | None = None) -> Verdict:
    """Write one certificate, then have it checked independently over the files."""
    cert = build_certificate(enc, rule, x, lo, hi, beta)
    cert_path = out / f"{name}.json"
    cert_path.write_text(json.dumps(cert, indent=2))
    ok, line = check_certificate(cert_path, out / "operator.json")
    print(f"  [{name:23s}] {line or '(no checker output)'}")
    return Verdict(name, rule, ok, line, float(lo), float(hi))


def run_case(case: str, out: Path) -> tuple[float, list[Verdict]]:
    """Emit every certificate this case supports, and check each one independently.

    Returns the exact electronic-frame lambda_min -- for comparison only, no certificate
    depends on it -- and one Verdict per certificate actually emitted.
    """
    mh = build_molecular_hamiltonian(**CASES[case])
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
            verdicts.append(emit(enc, "temple_inertia", x, lo, th0, out,
                                 "certificate_sector", beta=sector_eps))
        else:
            print(f"  [{'certificate_sector':23s}] not emitted: sector eps {sector_eps:.6f}"
                  f" <= theta {th0:.6f}, so Temple gives no finite bound")
        # (2) same witness, separator taken over the whole matrix
        lo = th0 - var0 / (beta - th0)
        pad = pad_claim(th0, 1e-9, len(x), th0 - lo)
        verdicts.append(emit(enc, "temple_inertia", x, lo - pad, th0 + pad, out,
                             "certificate_temple", beta=beta))
    else:
        print(f"  n > DENSE_LIMIT ({DENSE_LIMIT}): no inertia route, so no gap discharge")

    # The loose route needs no gap and no factorisation, and always applies.
    lo = gershgorin_lower(enc)
    pad = pad_claim(th0, 1e-9, len(x), th0 - lo)
    verdicts.append(emit(enc, "gershgorin_rayleigh", x, lo - pad, th0 + pad, out,
                         "certificate_gershgorin"))

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
    # so rather than let a reader assume the old exclusion still holds.
    return 0 if any(v.ok for v in verdicts) else 1


if __name__ == "__main__":
    sys.exit(main())
