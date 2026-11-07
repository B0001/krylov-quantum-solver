"""
Acceptance gates G1-G5 for specs/SPEC_certkit_regression_gate.md (the regression-certificate
CI gate).

Claim: every certificate this repo's producer emits still verifies under the exact certkit
commit pinned in pyproject.toml, and the energy each verified certificate certifies is still
within chemical accuracy of the exact lowest eigenvalue of the qubit matrix -- the GLOBAL one
over every symmetry sector, not the HF-reachable minimum SPEC_temple_bracket.md brackets; a
lambda_min_enclosure is a claim about the whole matrix. The checker runs out of process over
the files, per certkit's INTEGRATION.md -- nothing here imports check(). "Every certificate"
includes, per case, the one the solver's own return path emits (solve/<route>, portfolio-cir).

Why the checker is not enough on its own, and why G4 exists: the gershgorin_rayleigh route
certifies [gershgorin_lower(H), <x|H|x>], a true enclosure for ANY unit vector x. A witness of
pure random noise, 0.84 Ha above the ground state, is VERIFIED. G5 demonstrates that rather
than assuming it. Soundness comes from the checker; that the energy is still correct comes
from comparing the certified upper bound against a reference. Both are needed.

Electronic frame throughout: certificates are about the qubit matrix, and mh.energy_offset is
a trusted addition made outside it. The references below are therefore NOT the total energies
the rest of the suite asserts.

Needs the pinned checker (`uv sync --inexact --extra certkit`). PySCF/qiskit, no block2;
`make gates` runs it in its own process.
"""
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pytest

if importlib.util.find_spec("certkit") is None:
    # Locally the extra is optional, and skipping says so. In CI its absence is the gate not
    # running at all, which must never be reported as a pass. GitHub Actions sets CI=true.
    if os.environ.get("CI"):
        raise RuntimeError(
            "certkit is not installed: the regression-certificate gate cannot run. "
            'Install the pinned checker with `uv sync --inexact --extra certkit`.'
        )
    pytest.skip('certkit extra not installed (uv sync --inexact --extra certkit)',
                allow_module_level=True)

from certkit.producer import pad_claim  # noqa: E402  -- only once the checker is present

import certkit_bridge  # noqa: E402
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian  # noqa: E402
from temple_bounds import mean_and_variance  # noqa: E402

# Ha. The threshold this repo measures energy claims against; ~1.6 mHa.
CHEMICAL_ACCURACY = 1.6e-3

# Where the producer writes. Resolved from this file, not the CWD, so the gate behaves the
# same under `make gates`, `make test`, and a bare pytest run.
OUT = Path(__file__).resolve().parent.parent / "certkit_out"

# case -> (exact electronic-frame lambda_min, {certificate name: expected checker status}).
#
# Provenance: dense diagonalization of the qubit Hamiltonian,
# mh.ground_state_energy() - mh.energy_offset, which tests/test_reference_energies.py
# cross-checks against an independent PySCF FCI solve in the total frame. G3 re-derives these
# at test time as well as comparing to the literal, so neither can drift alone.
#
# The certificate set is pinned, not just the verdicts: a route that quietly stops emitting
# would otherwise read as a pass. The sector certificates are expected to ABSTAIN -- their
# eps comes from the solver's own second Ritz pair and the checker cannot discharge that
# premise by inertia count, which is the honest outcome and worth pinning in both directions.
# Only the status is pinned, never the reason string: the reason varies with Pauli term
# ordering, which is not stable across processes (PYTHONHASHSEED).
# Each certificate maps to (expected checker status, how far below exact its certified floor
# may sit). The floor needs its own pin because the checker accepts ANY lower bound below the
# one it re-derives -- a certificate claiming lambda_min > -1e6 Ha is perfectly VERIFIED and
# perfectly useless. Values are ~2x what the producer emits today.
#
# solve/<route> is what the solver's own return path emits for the energy it returns
# (QuantumKrylovSolver.solve(certify=), portfolio-cir): certify_ground_state, whose ceiling is
# the returned energy itself, so G4's upper end gates the solver's answer and not a witness.
# It returns the first route that verifies, so pinning the name pins which route won. Its
# floors are the same numbers as run_case's for that route, and so are the pins.
EXPECTED = {
    "H2": (-1.8523881735695822, {
        "certificate_sector": ("ABSTAIN", None),
        "certificate_temple": ("VERIFIED", 1e-8),
        "certificate_gershgorin": ("VERIFIED", 0.2),
        "solve/certificate_temple": ("VERIFIED", 1e-8),
    }),
    "H4": (-4.728206889123854, {
        "certificate_sector": ("ABSTAIN", None),
        "certificate_temple": ("VERIFIED", 1e-4),
        "certificate_gershgorin": ("VERIFIED", 1.0),
        "solve/certificate_temple": ("VERIFIED", 1e-4),
    }),
    # Self mode gives no finite Temple bound here (eps <= theta), so no sector certificate is
    # emitted at all -- certificate-or-abstention, applied by the producer.
    "H4-stretched": (-3.04433126964987, {
        "certificate_temple": ("VERIFIED", 2e-2),
        "certificate_gershgorin": ("VERIFIED", 1.5),
        "solve/certificate_temple": ("VERIFIED", 2e-2),
    }),
    # 12 qubits: n = 4096 exceeds certkit's DENSE_LIMIT, so there is no inertia route and the
    # loose Gershgorin enclosure (2.6 Ha wide) is all this case can certify. G4 is the only
    # thing standing between that certificate and vacuity, at both ends.
    "N2": (-11.392972038179622, {
        "certificate_gershgorin": ("VERIFIED", 3.0),
        "solve/certificate_gershgorin": ("VERIFIED", 3.0),
    }),
}

_CACHE = {}


def _case(name):
    """(exact electronic lambda_min, [Verdict]) for one case. ~30 s for all four, so once."""
    if name not in _CACHE:
        _CACHE[name] = certkit_bridge.run_case(name, OUT / name)
    return _CACHE[name]


def test_G1_every_certificate_gets_the_verdict_it_is_pinned_to():
    """The independent checker's status on every emitted certificate matches the pin."""
    for case, (_, expected) in EXPECTED.items():
        _, verdicts = _case(case)
        for v in verdicts:
            got = "VERIFIED" if v.ok else "ABSTAIN"
            # Exit 1 means ABSTAIN *or* a crash, and a crash must never be accepted as an
            # expected abstention. The checker's own word has to agree with its exit code;
            # check_certificate makes a crash's word CRASH (test_G1_crash_guard below).
            assert v.line.startswith(got), (
                f"{case}/{v.name} ({v.rule}): exit code says {got}, the checker's output "
                f"does not -- {v.line}")
            assert got == expected[v.name][0], (
                f"{case}/{v.name} ({v.rule}): checker said {got}, pinned "
                f"{expected[v.name][0]} -- {v.line}")


def test_G1_crash_guard_a_crashing_checker_is_never_read_as_an_abstention(tmp_path, monkeypatch):
    """G1's crash guard holds for each way the pinned checker can crash.

    Since certkit-3zr an exception inside a rule comes out as "ABSTAIN   internal error: ...",
    exit 1 -- exactly what the sector certificates are pinned to -- and the CLI's own errors
    exit 2. Both are injected into the real checker process here and must come back CRASH.
    Without the first, a checker crashing on every sector certificate passes G1.
    """
    _case("H2")
    cert = OUT / "H2" / "certificate_sector.json"
    op_ref = json.loads(cert.read_text())["claim"]["operator_ref"]

    (tmp_path / "sitecustomize.py").write_text(
        "import certkit.checker as c\n"
        "def _boom(*a, **k):\n"
        "    raise RuntimeError('injected by the gate')\n"
        "kind, _, needs_vector = c.RULES['temple_inertia']\n"
        "c.RULES['temple_inertia'] = (kind, _boom, needs_vector)\n")
    monkeypatch.setenv("PYTHONPATH", os.pathsep.join(
        filter(None, [str(tmp_path), os.environ.get("PYTHONPATH")])))
    ok, line = certkit_bridge.check_certificate(cert, OUT / "H2" / "operator.json", op_ref)
    assert not ok and line.startswith("CRASH") and "internal error" in line, line

    monkeypatch.undo()
    ok, line = certkit_bridge.check_certificate(cert, tmp_path / "missing.json", op_ref)
    assert not ok and line.startswith("CRASH  exit 2"), line


def test_G2_the_emitted_certificate_set_is_exactly_the_pinned_set():
    """A route that silently stops emitting is a loss of coverage, not a pass."""
    for case, (_, expected) in EXPECTED.items():
        _, verdicts = _case(case)
        assert {v.name for v in verdicts} == set(expected), (
            f"{case}: emitted {sorted(v.name for v in verdicts)}, pinned {sorted(expected)}")


def test_G3_every_verified_enclosure_contains_the_exact_energy():
    """DEFINITION OF DONE: soundness. One containment escape kills the claim."""
    for case, (pinned, _) in EXPECTED.items():
        lam, verdicts = _case(case)
        assert abs(lam - pinned) < 1e-6, (
            f"{case}: reference drifted -- re-derived {lam!r}, pinned {pinned!r}")
        for v in verdicts:
            if v.ok:
                assert v.lo <= lam <= v.hi, (
                    f"{case}/{v.name} ({v.rule}): exact {lam!r} escaped the certified "
                    f"enclosure [{v.lo!r}, {v.hi!r}]")


def test_G4_every_verified_enclosure_is_still_tight_at_both_ends():
    """What makes this a regression gate and not only a soundness check.

    The certified upper bound is the solver's energy, so it is gated against chemical
    accuracy. The certified floor is gated too, and separately: the checker accepts any
    lower bound below the one it re-derives, so an enclosure can stay sound while going
    slack enough to say nothing. Both ends can be VERIFIED and worthless -- see G5.
    """
    for case, (_, expected) in EXPECTED.items():
        lam, verdicts = _case(case)
        for v in verdicts:
            if not v.ok:
                continue
            assert v.hi - lam <= CHEMICAL_ACCURACY, (
                f"{case}/{v.name} ({v.rule}): certified upper bound {v.hi!r} is "
                f"{v.hi - lam:.3e} Ha above exact {lam!r}")
            max_slack = expected[v.name][1]
            assert lam - v.lo <= max_slack, (
                f"{case}/{v.name} ({v.rule}): certified floor {v.lo!r} sits "
                f"{lam - v.lo:.3e} Ha below exact, pinned at {max_slack:.3e}")


def test_G5_a_noise_witness_is_verified_by_the_checker_and_caught_by_G4():
    """Non-vacuity, demonstrated rather than asserted.

    A Gershgorin certificate built on random noise is a true enclosure, so the checker
    verifies it. If certkit ever grows strong enough to reject it, this fails -- and that is
    a finding about the checker, not a regression here.
    """
    lam, _ = _case("H2")
    mh = build_molecular_hamiltonian(**certkit_bridge.CASES["H2"])
    H = mh.qubit_hamiltonian.to_matrix(sparse=True).tocsc()
    enc = certkit_bridge.encode_hamiltonian(mh)

    out = OUT / "H2-noise"
    out.mkdir(parents=True, exist_ok=True)
    (out / "operator.json").write_text(json.dumps(enc, indent=2))

    x = np.random.default_rng(0).standard_normal(1 << mh.qubit_hamiltonian.num_qubits)
    x /= np.linalg.norm(x)
    theta, _var = mean_and_variance(H, x)
    lo = certkit_bridge.gershgorin_lower(enc)
    # Padded exactly as the producer pads. Without it the claim comes out tighter than the
    # checker's re-derivation and abstains on rounding alone, which would prove nothing.
    pad = pad_claim(theta, 1e-9, len(x), theta - lo)
    v = certkit_bridge.emit(enc, "gershgorin_rayleigh", x, lo - pad, theta + pad,
                            out, "certificate_noise")

    assert v.ok, f"the checker no longer verifies a noise witness -- {v.line}"
    assert v.hi - lam > CHEMICAL_ACCURACY, (
        f"G4 would not have caught a noise witness: its certified upper bound {v.hi!r} is "
        f"only {v.hi - lam:.3e} Ha above exact {lam!r}")
