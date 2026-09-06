"""
Acceptance gates G1-G5 for specs/SPEC_certkit_regression_gate.md (the regression-certificate
CI gate).

Claim: every certificate this repo's producer emits still verifies under the exact certkit
release pinned in pyproject.toml, and the energy each verified certificate certifies is still
within chemical accuracy of the exact reachable minimum. The checker runs out of process over
the files, per certkit's INTEGRATION.md -- nothing here imports check().

Why the checker is not enough on its own, and why G4 exists: the gershgorin_rayleigh route
certifies [gershgorin_lower(H), <x|H|x>], a true enclosure for ANY unit vector x. A witness of
pure random noise, 0.84 Ha above the ground state, is VERIFIED. G5 demonstrates that rather
than assuming it. Soundness comes from the checker; that the energy is still correct comes
from comparing the certified upper bound against a reference. Both are needed.

Electronic frame throughout: certificates are about the qubit matrix, and mh.energy_offset is
a trusted addition made outside it. The references below are therefore NOT the total energies
the rest of the suite asserts.

Needs the pinned checker (`uv pip install -e ".[certkit]"`). PySCF/qiskit, no block2;
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
            'Install the pinned checker with `uv pip install -e ".[certkit]"`.'
        )
    pytest.skip('certkit extra not installed (uv pip install -e ".[certkit]")',
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
EXPECTED = {
    "H2": (-1.8523881735695822, {
        "certificate_sector": "ABSTAIN",
        "certificate_temple": "VERIFIED",
        "certificate_gershgorin": "VERIFIED",
    }),
    "H4": (-4.728206889123854, {
        "certificate_sector": "ABSTAIN",
        "certificate_temple": "VERIFIED",
        "certificate_gershgorin": "VERIFIED",
    }),
    # Self mode gives no finite Temple bound here (eps <= theta), so no sector certificate is
    # emitted at all -- certificate-or-abstention, applied by the producer.
    "H4-stretched": (-3.04433126964987, {
        "certificate_temple": "VERIFIED",
        "certificate_gershgorin": "VERIFIED",
    }),
    # 12 qubits: n = 4096 exceeds certkit's DENSE_LIMIT, so there is no inertia route and the
    # loose Gershgorin enclosure (2.6 Ha wide) is all this case can certify. G4 is the only
    # thing standing between that certificate and vacuity.
    "N2": (-11.392972038179622, {"certificate_gershgorin": "VERIFIED"}),
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
            # Exit 1 means ABSTAIN or a crash; only the output line tells them apart.
            assert v.line, f"{case}/{v.name}: checker produced no verdict line -- it crashed"
            got = "VERIFIED" if v.ok else "ABSTAIN"
            assert got == expected.get(v.name), (
                f"{case}/{v.name} ({v.rule}): checker said {got}, pinned "
                f"{expected.get(v.name)} -- {v.line}")


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


def test_G4_every_verified_upper_bound_is_within_chemical_accuracy():
    """What makes this a regression gate and not only a soundness check.

    The certified upper bound is the solver's energy. The checker will happily verify a
    sound enclosure around a badly wrong one -- see G5.
    """
    for case, (_, _expected) in EXPECTED.items():
        lam, verdicts = _case(case)
        for v in verdicts:
            if v.ok:
                assert v.hi - lam <= CHEMICAL_ACCURACY, (
                    f"{case}/{v.name} ({v.rule}): certified upper bound {v.hi!r} is "
                    f"{v.hi - lam:.3e} Ha above exact {lam!r}")


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
