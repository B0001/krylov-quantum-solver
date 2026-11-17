#!/usr/bin/env python3
"""
certification.py -- certificate-or-abstention for every ground-state energy a solver returns.

specs/SPEC_certificate_or_abstention.md (portfolio-0pc). Two things live here, and every
energy-returning path of this package goes through both:

  * the PRODUCER: ``certify_energies`` wraps a returned energy and the state it is the Ritz
    value of in a certkit ``lambda_min_enclosure`` certificate, writes it to disk, and has
    certkit's checker judge the files in a separate process;
  * the GATE on the number: ``certified`` hands an energy out only when that checker said
    VERIFIED, and raises ``EnergyAbstained`` otherwise. The result types read their
    ``energy`` through it, so an unverified number can only be reached under a name that
    says ``unverified``.

Nothing in this file is trusted, and nothing in it adjudicates: per certkit's INTEGRATION.md
the checker is consumed as a protocol, never imported. certkit is a core dependency, imported
here at import time, so a missing checker is a broken install and not a silent skip.

The fixed regression cases (H2, H4, H4-stretched, N2) and their CLI stay in the top-level
``certkit_bridge.py``, which re-exports the producer names it used to define.
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Callable, NamedTuple, Sequence

import numpy as np
from certkit.interval import Iv
from certkit.operators import DENSE_LIMIT, decode_operator, encode_pauli, operator_ref
from certkit.producer import pad_claim     # public producer-side padding convention
from certkit.schema import f2h, seal

# ponytail: one fixed ceiling for every system size. The slowest check measured is 37 min
# (16 qubits, under load; specs/SPEC_certificate_or_abstention.md section 3), so 2 h fails
# only a hung checker. A constant, not a setting: G1 pins that nothing outside the call can
# change the return path. Make it scale with the operator if larger systems need longer
# honest checks.
CHECK_TIMEOUT_S = 7200


class Verdict(NamedTuple):
    """One certificate, and what the independent checker said about it.

    ``line`` is the checker's first stdout line, verbatim, or a ``CRASH``-prefixed
    diagnostic if it gave no verdict at all (see ``check_certificate``). A crash can exit 1
    just as an ABSTAIN does, so the line is the only thing that tells them apart.
    """

    name: str
    rule: str
    ok: bool
    line: str
    lo: float
    hi: float


class EnergyAbstained(RuntimeError):
    """Raised on reading an energy the independent checker did not VERIFY.

    Loud by design (the house style of specs/SPEC_certified_overlap_bounds.md section 4):
    an abstention is never returned as a number. The solver's estimate is still on the
    result, as ``unverified_energy``, for a caller that says so.
    """


def certified(value, verdict: Verdict):
    """``value`` if the checker VERIFIED its certificate; raise ``EnergyAbstained`` if not.

    The one place a number becomes a certified energy. ``ok`` and the checker's own word
    must agree: a crash can share an exit code with a verdict, never a line.
    """
    if verdict.ok and verdict.line.startswith("VERIFIED"):
        return value
    raise EnergyAbstained(
        f"no certified energy: {verdict.line or '(no checker output)'}. The solver's "
        "estimate is on the result as unverified_energy; nothing vouches for it.")


def require_verdict(result) -> None:
    """Refuse to build a result that has no checker verdict (called from ``__post_init__``)."""
    if not isinstance(result.certificate, Verdict):
        raise TypeError(
            f"{type(result).__name__} needs the checker's Verdict as its certificate, got "
            f"{result.certificate!r}: no energy without certificate-or-abstention.")


def not_submitted(reason: str) -> Verdict:
    """The abstention of the uncertified door: no certificate was written or checked."""
    nan = float("nan")
    return Verdict("not_submitted", "", False,
                   f"ABSTAIN   not submitted to the checker: {reason}", nan, nan)


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


def full_space_beta(H):
    """A separator for the WHOLE matrix, not the reachable Krylov sector.

    Untrusted like every producer number: the checker discharges beta by inertia
    count, so a wrong guess costs coverage and never soundness. Dense on purpose: it is
    only called when certkit will materialise the matrix (n <= DENSE_LIMIT), and unlike a
    Lanczos solve it has no random start vector and cannot fail to converge -- this runs on
    every return path, so it must not be more fragile than the solve it certifies.
    """
    vals = np.linalg.eigvalsh(H.toarray().real)[:2]
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


def check_certificate(cert_path: Path, operator_path: Path, op_ref: str) -> tuple[bool, str]:
    """Run the checker as a separate process over the files, and report what it said.

    The verdict is the checker's exit status, on the question certkit's INTEGRATION.md says
    to ask: is this a lambda_min_enclosure (--require-kind) of the operator the producer
    encoded, ``op_ref`` (--operator-ref)? A true claim of another kind, or a swapped operator
    file, then abstains instead of verifying. This file is a producer and is untrusted; it
    must not adjudicate its own claims in its own process.

    A crash must never pass for an expected abstention, and the checker can crash three
    ways: exit 2 (no verdict, certkit-3zr); dying before it prints, which still exits 1;
    and, also since certkit-3zr, an exception inside a rule, which it prints as an ABSTAIN
    whose reason is "internal error: ...". All three come back CRASH-prefixed, never as a
    bare line that a caller matching on "ABSTAIN" would read as the expected refusal.

    A checker still running after ``CHECK_TIMEOUT_S`` is killed and reported the same way,
    so a hung checker abstains instead of hanging the solver.
    """
    try:
        proc = subprocess.run(
            [sys.executable, "-m", "certkit.cli", "check", "--require-kind",
             "lambda_min_enclosure", "--operator-ref", op_ref, str(cert_path),
             str(operator_path)],
            capture_output=True, text=True, timeout=CHECK_TIMEOUT_S,
        )
    except subprocess.TimeoutExpired:
        return False, f"CRASH  timeout: no verdict within {CHECK_TIMEOUT_S} s"
    line = proc.stdout.strip().split("\n", 1)[0]
    internal = line.removeprefix("ABSTAIN").lstrip().startswith("internal error:")
    if line and proc.returncode in (0, 1) and not internal:
        return proc.returncode == 0, line
    err = proc.stderr.strip().splitlines()
    why = line or (err[-1] if err else "(no output at all)")
    return False, f"CRASH  exit {proc.returncode}: {why}"


def emit(enc: dict, rule: str, x, lo: float, hi: float, out: Path, name: str,
         beta: float | None = None) -> Verdict:
    """Write one certificate, then have it checked independently over the files."""
    cert = build_certificate(enc, rule, x, lo, hi, beta)
    cert_path = out / f"{name}.json"
    cert_path.write_text(json.dumps(cert, indent=2))
    ok, line = check_certificate(cert_path, out / "operator.json", operator_ref(enc))
    return Verdict(name, rule, ok, line, float(lo), float(hi))


def real_witness(H, psi):
    """The real unit vector a certificate can carry, from a complex solver state.

    H is real symmetric, so <psi|H|psi> is the norm-weighted mean of the Rayleigh quotients
    of Re(psi) and Im(psi): the lower of the two is never above psi's own. Picking by weight
    instead (as certkit_bridge.krylov_witness does) can land above the energy the solver
    returned.
    """
    parts = [p / np.linalg.norm(p) for p in (psi.real, psi.imag) if np.linalg.norm(p) > 0]
    return min(parts, key=lambda p: float(p @ (H @ p)))


def certify_ground_state(mh, psi, energy: float, out: Path) -> Verdict:
    """Certificate for one ground-state energy a solver returned, checked out of process.

    ``energy`` is the returned value in the ELECTRONIC frame and ``psi`` the state it is the
    Ritz value of (portfolio-cir). The claim is lambda_min in [floor, energy]: the certified
    ceiling is the returned number itself, padded, so VERIFIED means that number is a proven
    upper bound on the ground state -- not merely that some vector has a low Rayleigh
    quotient (G5 of the regression gate shows how little that is worth). A noisy energy
    below its own witness's Rayleigh quotient therefore ABSTAINs, as does one below the
    floor (an inverted enclosure).

    Routes, tightest first: temple_inertia when certkit will materialise the matrix and
    Temple's bound is finite, then gershgorin_rayleigh. The first VERIFIED wins; otherwise
    the last abstention is returned. Writes ``out/operator.json`` plus one file per route
    tried.
    """
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    enc = encode_hamiltonian(mh)
    (out / "operator.json").write_text(json.dumps(enc, indent=2))
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).real.tocsc()
    x = real_witness(H, np.asarray(psi))
    # Mean and variance of the witness, from one matvec (temple_bounds.mean_and_variance,
    # which lives outside the package and imports the solver).
    hx = H @ x
    th = float(x @ hx)
    var = max(float(hx @ hx) - th ** 2, 0.0)

    routes = []
    if len(x) <= DENSE_LIMIT:
        beta, _ = full_space_beta(H)
        if th < beta:
            routes.append(("temple_inertia", th - var / (beta - th), beta, "certificate_temple"))
    routes.append(("gershgorin_rayleigh", gershgorin_lower(enc), None, "certificate_gershgorin"))
    for rule, lo, beta, name in routes:
        pad = pad_claim(energy, 1e-9, len(x), energy - lo)
        v = emit(enc, rule, x, lo - pad, energy + pad, out, name, beta=beta)
        if v.ok:
            break
    return v


def certify_energies(mh, candidates: Sequence[tuple[float, Callable[[], np.ndarray]]],
                     out: str | Path | None = None,
                     names: Sequence[str] | None = None) -> list[Verdict]:
    """One checker ``Verdict`` per returned energy, from as few checker runs as the proof needs.

    ``candidates`` are ``(electronic-frame energy, zero-argument callable returning the state
    it is the Ritz value of)`` for the energies ONE call returns for ONE Hamiltonian; the
    verdicts come back in the same order. Every energy-returning path of this package ends
    here, and nothing about a call can make it skip the checker.

    The lowest energy is submitted first. A VERIFIED certificate for energy ``E`` proves
    ``lambda_min <= E + pad``, which makes every energy ``E' >= E`` of the same operator a
    proven upper bound too, so those share that verdict instead of paying for a checker run
    each: a noiseless sweep costs one run, not one per step. Their certificate's ceiling is
    then at or below their own energy, not equal to it. An energy below the first verified
    one keeps its own abstention, with no enclosure attached.

    ``out`` is where the certificate files are kept (``names[i]`` a subdirectory per
    candidate); without it they go to a temporary directory that is removed afterwards.
    """
    order = sorted(range(len(candidates)), key=lambda i: candidates[i][0])
    verdicts: list = [None] * len(candidates)
    nan = float("nan")
    with tempfile.TemporaryDirectory(prefix="certkit-") as scratch:
        base = Path(scratch if out is None else out)
        for k, i in enumerate(order):
            energy, witness = candidates[i]
            v = certify_ground_state(mh, witness(), float(energy),
                                     base / names[i] if names else base)
            if v.ok:
                for j in order[k:]:
                    verdicts[j] = v
                break
            verdicts[i] = v._replace(lo=nan, hi=nan)   # a refused claim encloses nothing
    return verdicts
