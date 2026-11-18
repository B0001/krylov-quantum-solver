"""
Acceptance gates G1-G5 for specs/SPEC_certificate_or_abstention.md (portfolio-0pc).

Claim: every public path of hybrid_quantum_solver that returns a ground-state energy hands
the caller either the energy with a certkit certificate the pinned checker VERIFIED out of
process, or an explicit abstention -- never a bare number, and no argument, default or
environment variable that skips the checker while still returning something that reads as a
certified energy.

The whole gate turns on one detector, ``_assert_certificate_or_abstention``: G2 runs it on
certified results, G3 on abstained ones, and G5 shows it rejects four mutants, including a
path patched to return the bare Ritz value.

Electronic frame for certificates, as in the sibling gates: lambda_min is the GLOBAL lowest
eigenvalue of the qubit matrix and mh.energy_offset is added outside it.

PySCF/qiskit and the pinned checker (a core dependency since this spec). No block2.
"""
import dataclasses
import inspect
import json
import math
import shutil
import subprocess
import sys
import types
from pathlib import Path

import numpy as np
import pytest
from certkit.operators import operator_ref
from certkit.schema import h2f

from hybrid_quantum_solver import (
    HardwareKrylovSolver,
    QuantumKrylovSolver,
    TrotterKrylovSolver,
    build_hamiltonian_from_integrals,
    build_molecular_hamiltonian,
    certification,
    hardware_krylov,
    hubbard_chain_integrals,
    pipeline,
    quantum_krylov_solver,
    run_from_integrals,
    run_geometry,
    skqd,
    trotter_krylov,
)
from hybrid_quantum_solver.certification import EnergyAbstained, Verdict
from hybrid_quantum_solver.skqd import SampleKrylovSolver

H2 = "H 0 0 0; H 0 0 0.74"
M = 3
# pad_claim's float padding on the certified ceiling is 1e-9 relative; past 1e-7 Ha the
# certificate is no longer about the returned number. Same constant as the sibling gate.
PAD_MAX = 1e-7
# Written where the sibling gates write, so a CI failure can upload the certificates.
OUT = Path(__file__).resolve().parent.parent / "certkit_out" / "certificate_or_abstention"

_CACHE = {}


def _h2():
    if "h2" not in _CACHE:
        _CACHE["h2"] = build_molecular_hamiltonian(atom=H2)
    return _CACHE["h2"]


def _dimer():
    """Hubbard dimer (U=4) for run_from_integrals: the global lambda_min lies in another
    particle sector, so the returned energy is a true but loose upper bound (Gershgorin)."""
    if "dimer" not in _CACHE:
        mi = hubbard_chain_integrals(2, 4.0)
        mh = build_hamiltonian_from_integrals(mi.h1, mi.eri, num_particles=mi.nelec,
                                              energy_offset=mi.e_core)
        _CACHE["dimer"] = (mi, mh)
    return _CACHE["dimer"]


def _from_integrals(d, **kw):
    mi, mh = _dimer()
    return mh, run_from_integrals(mi.h1, mi.eri, mi.nelec, mi.e_core, krylov_dim=M,
                                  certificate_dir=d, **kw)


# Every public ENERGY path: name -> (its parameters, how the gate calls it). ``call(d)``
# returns (mh, result) with certificates kept in directory ``d`` (None: a temporary one).
# The parameter names are pinned on purpose: a new argument on an energy path is exactly
# where a switch would appear, so it has to be looked at here before the gate goes green.
_PIPELINE_COMMON = ("krylov_dim", "dt", "shots", "seed")
PATHS = {
    "QuantumKrylovSolver.solve": (
        ("krylov_dim", "certificate_dir"),
        lambda d: (_h2(), QuantumKrylovSolver(_h2()).solve(M, d))),
    "QuantumKrylovSolver.convergence": (
        ("max_dim", "certificate_dir"),
        lambda d: (_h2(), QuantumKrylovSolver(_h2()).convergence(M, d))),
    "QuantumKrylovSolver.solve_excited": (
        ("krylov_dim", "n_states", "certificate_dir"),
        lambda d: (_h2(), QuantumKrylovSolver(_h2()).solve_excited(M, 2, d))),
    "QuantumKrylovSolver.eigenstates": (
        ("krylov_dim", "n_states", "certificate_dir"),
        lambda d: (_h2(), QuantumKrylovSolver(_h2()).eigenstates(M, 2, d))),
    "TrotterKrylovSolver.solve": (
        ("krylov_dim", "certificate_dir"),
        lambda d: (_h2(), TrotterKrylovSolver(_h2(), trotter_reps=2).solve(M, d))),
    "TrotterKrylovSolver.convergence": (
        ("max_dim", "certificate_dir"),
        lambda d: (_h2(), TrotterKrylovSolver(_h2(), trotter_reps=2).convergence(M, d))),
    "HardwareKrylovSolver.solve": (
        ("krylov_dim", "certificate_dir"),
        lambda d: (_h2(), HardwareKrylovSolver(_h2()).solve(M, d))),
    "HardwareKrylovSolver.convergence": (
        ("max_dim", "certificate_dir"),
        lambda d: (_h2(), HardwareKrylovSolver(_h2()).convergence(M, d))),
    "SampleKrylovSolver.solve": (
        ("certificate_dir",),
        lambda d: (_h2(), SampleKrylovSolver(_h2(), depth=4, n_shots=50, seed=0).solve(d))),
    "SampleKrylovSolver.convergence": (
        ("shot_schedule", "certificate_dir"),
        lambda d: (_h2(), SampleKrylovSolver(_h2(), depth=4, seed=0).convergence((5, 50), d))),
    "pipeline.run_geometry": (
        ("atom", "basis", "charge", "spin", "active_electrons", "active_orbitals",
         *_PIPELINE_COMMON, "reference", "track_convergence", "conv_tol", "certificate_dir",
         "uncertified_estimates_only"),
        lambda d: (_h2(), run_geometry(H2, krylov_dim=M, reference="none", certificate_dir=d))),
    "pipeline.run_from_integrals": (
        ("h1", "eri", "num_particles", "e_core", *_PIPELINE_COMMON, "reference_energy",
         "track_convergence", "certificate_dir", "uncertified_estimates_only"),
        _from_integrals),
}

# Every other public callable of the same modules, and why it is not an energy path.
NOT_ENERGY = {
    "QuantumKrylovSolver.uncertified_eigenstates":
        "the one uncertified door: returns abstentions, covered by G3",
    "HardwareKrylovSolver.measure_signal": "overlaps <phi0|U^k|phi0>; no Hamiltonian enters",
    "HardwareKrylovSolver.resource_report": "circuit and shot counts",
    "quantum_krylov_solver.expm_multiply_taylor": "a propagator: returns a state",
    "quantum_krylov_solver.ritz_pairs": "kernel: eigenpairs of matrices the CALLER supplies",
    "quantum_krylov_solver.ritz_spectrum": "kernel: eigenvalues of matrices the CALLER supplies",
    "quantum_krylov_solver.solve_generalized_eig":
        "kernel: lowest eigenvalue of matrices the CALLER supplies; it sees no operator and "
        "no state, so nothing could certify it (spec section 2, item 2)",
    "trotter_krylov.canonical_term_order": "reorders Pauli terms",
    "trotter_krylov.build_trotter_step": "returns a circuit",
    "trotter_krylov.estimate_energy_aer":
        "expectation value of a circuit the CALLER supplies, not a ground-state energy",
}

# Constructor arguments are pinned for the same reason as the path arguments.
CONSTRUCTORS = {
    QuantumKrylovSolver: ("molecular_hamiltonian", "dt", "threshold", "noise_sigma", "seed",
                          "device"),
    TrotterKrylovSolver: ("molecular_hamiltonian", "dt", "trotter_order", "trotter_reps",
                          "threshold"),
    HardwareKrylovSolver: ("molecular_hamiltonian", "dt", "trotter_order", "trotter_reps",
                           "threshold", "shots", "noise_model", "seed", "zne_scale_factors",
                           "zne_order"),
    SampleKrylovSolver: ("molecular_hamiltonian", "dt", "n_shots", "depth", "seed"),
}
RETURN_PATH_MODULES = (certification, quantum_krylov_solver, trotter_krylov, hardware_krylov,
                       skqd, pipeline)


# -- the detector ----------------------------------------------------------------------------
def _records(result):
    """The objects in a path's result that each carry one ground-state energy."""
    if isinstance(result, tuple):                 # eigenstates -> (spectrum step, states)
        result = result[0]
    if isinstance(result, list):
        return [r for item in result for r in _records(item)]
    if hasattr(result, "steps"):                  # PipelineResult
        return _records(result.steps)
    return [result]


def _raw(rec):
    """The solver's number, which must live under a name that says it is unverified."""
    try:
        if hasattr(rec, "unverified_energies"):
            return float(rec.unverified_energies[0])
        return float(rec.unverified_energy)
    except AttributeError as exc:
        raise AssertionError(f"{type(rec).__name__} does not name its number 'unverified': "
                             f"{exc}") from None


def _gated(rec):
    """The certified reading of the same number: returns it, or raises EnergyAbstained."""
    return rec.energies[0] if hasattr(rec, "unverified_energies") else rec.energy


# Classical references a PipelineResult also carries. They come from the SCF and from the
# caller, not from the solver, and are out of scope (spec section 2, item 2). ``hf_energy``
# has to be named here because it IS the M=1 Krylov energy: a one-vector Krylov space is the
# Hartree-Fock determinant.
REFERENCES = (".hf_energy", ".reference_energy")


def _leaks(obj, values, path="", depth=0):
    """Public attribute paths under which one of ``values`` can be read without the path
    saying ``unverified``. A gated accessor that raises EnergyAbstained is not a leak, and
    integers (dimensions, ranks, qubit counts) are not energies. "Read" means to within the
    certificate's padding: a refused claim's ceiling is the energy plus that padding, so an
    abstention that kept it would be handing the number back (found by mutation)."""
    found = []
    if isinstance(obj, (bool, int, np.integer, str, bytes, type(None))) or depth > 4:
        return found
    if path in REFERENCES:
        return found
    if isinstance(obj, (float, np.floating)):
        if "unverified" not in path and any(
                math.isfinite(v) and abs(float(obj) - v) <= PAD_MAX for v in values):
            found.append(path)
        return found
    if isinstance(obj, np.ndarray):
        return [] if obj.size > 16 else _leaks(obj.tolist(), values, path, depth + 1)
    if isinstance(obj, (list, tuple)):
        for i, item in enumerate(obj):
            found += _leaks(item, values, f"{path}[{i}]", depth + 1)
        return found
    for name in dir(obj):
        if name.startswith("_"):
            continue
        try:
            item = getattr(obj, name)
        except EnergyAbstained:
            continue
        if callable(item):
            continue
        found += _leaks(item, values, f"{path}.{name}", depth + 1)
    return found


def _recheck(verdict, mh, directory):
    """The certificate file behind a VERIFIED verdict, verified again by the gate itself.

    A verdict the package made up would have no file, or one the checker refuses. The
    operator file must be this Hamiltonian's, not merely some operator the claim is true of.
    """
    assert directory is not None, "a VERIFIED verdict was returned with no directory to audit"
    want = operator_ref(certification.encode_hamiltonian(mh))
    for cert_path in sorted(Path(directory).rglob(f"{verdict.name}.json")):
        enclosure = json.loads(cert_path.read_text())["claim"]["enclosure"]
        if (h2f(enclosure["lo"]), h2f(enclosure["hi"])) != (verdict.lo, verdict.hi):
            continue
        op_path = cert_path.parent / "operator.json"
        assert operator_ref(json.loads(op_path.read_text())) == want, (
            f"{cert_path}: certificate is not about this Hamiltonian")
        proc = subprocess.run(
            [sys.executable, "-m", "certkit.cli", "check", "--require-kind",
             "lambda_min_enclosure", "--operator-ref", want, str(cert_path), str(op_path)],
            capture_output=True, text=True)
        assert proc.returncode == 0 and proc.stdout.startswith("VERIFIED"), (
            f"{cert_path}: the gate's own check says {proc.stdout.strip() or proc.stderr!r}")
        return
    raise AssertionError(f"VERIFIED verdict {verdict.name} [{verdict.lo!r}, {verdict.hi!r}] "
                         f"has no certificate file under {directory}")


def _assert_certificate_or_abstention(result, mh, directory=None):
    """Every energy in ``result`` is either certified or an explicit abstention.

    Returns (raw total energies, verdicts). Raises AssertionError on anything else: a bare
    number, a missing verdict, a certified reading that is not the solver's number, a
    certificate that is not about it, a VERIFIED the gate cannot reproduce, or an abstention
    whose number can be read under a name that does not say ``unverified``.
    """
    records = _records(result)
    assert records, "the path returned no energy at all"
    lam = mh.ground_state_energy() - mh.energy_offset
    raws, verdicts = [], []
    for rec in records:
        assert not isinstance(rec, (float, int, np.floating)), (
            f"bare energy {rec!r}: no certificate and no abstention")
        verdict = getattr(rec, "certificate", None)
        assert isinstance(verdict, Verdict), (
            f"{type(rec).__name__} carries {verdict!r} where the checker's Verdict belongs")
        raw = _raw(rec)
        if verdict.ok:
            assert verdict.line.startswith("VERIFIED"), verdict.line
            assert _gated(rec) == raw, "the certified energy is not the solver's number"
            e = raw - mh.energy_offset
            assert verdict.hi - e <= PAD_MAX, (
                f"certified ceiling {verdict.hi!r} is above the returned energy {e!r}: "
                "the certificate does not prove this number an upper bound")
            assert verdict.lo <= lam <= verdict.hi, (
                f"exact {lam!r} escaped the certified [{verdict.lo!r}, {verdict.hi!r}]")
            _recheck(verdict, mh, directory)
        else:
            assert verdict.line.startswith(("ABSTAIN", "CRASH")), verdict.line
            with pytest.raises(EnergyAbstained):
                _gated(rec)
        raws.append(raw)
        verdicts.append(verdict)
    if not all(v.ok for v in verdicts):
        # Total and electronic frame both: an abstention must not carry its number anywhere
        # but under ``unverified``, including inside the verdict of the refused claim.
        refused = [r for r, v in zip(raws, verdicts) if not v.ok]
        leaks = _leaks(result, refused + [r - mh.energy_offset for r in refused])
        assert not leaks, f"an abstained energy is readable as {leaks}"
    return raws, verdicts


def _certified(name):
    """(mh, result, directory) for one path with the real checker, computed once."""
    if name not in _CACHE:
        directory = OUT / name
        shutil.rmtree(directory, ignore_errors=True)   # only this gate writes here
        mh, result = PATHS[name][1](directory)
        _CACHE[name] = (mh, result, directory)
    return _CACHE[name]


# -- G1 --------------------------------------------------------------------------------------
def _public_callables():
    """Every public function and method defined in the solver and pipeline modules."""
    found = {}
    for module in RETURN_PATH_MODULES[1:]:
        short = module.__name__.rsplit(".", 1)[-1]
        for name, obj in vars(module).items():
            if name.startswith("_") or getattr(obj, "__module__", None) != module.__name__:
                continue
            if inspect.isfunction(obj):
                found[f"{short}.{name}"] = obj
            elif inspect.isclass(obj) and not dataclasses.is_dataclass(obj):
                for attr, member in vars(obj).items():
                    if not attr.startswith("_") and callable(member):
                        found[f"{name}.{attr}"] = member
    return found


def test_G1_the_public_surface_is_pinned_and_reads_no_environment():
    """Every public callable is classified, every energy path's arguments are the pinned
    ones, and nothing on the return path consults the environment."""
    public = _public_callables()
    assert set(public) == set(PATHS) | set(NOT_ENERGY), (
        f"unclassified public callables: {sorted(set(public) - set(PATHS) - set(NOT_ENERGY))}; "
        f"listed but gone: {sorted((set(PATHS) | set(NOT_ENERGY)) - set(public))}")
    for name, (params, _) in PATHS.items():
        got = tuple(p for p in inspect.signature(public[name]).parameters if p != "self")
        assert got == params, f"{name}: arguments {got}, pinned {params}"
    for cls, params in CONSTRUCTORS.items():
        got = tuple(p for p in inspect.signature(cls.__init__).parameters if p != "self")
        assert got == params, f"{cls.__name__}(): arguments {got}, pinned {params}"
    for module in RETURN_PATH_MODULES:
        source = inspect.getsource(module)
        assert "environ" not in source and "getenv" not in source, (
            f"{module.__name__} reads the environment")


# -- G2 --------------------------------------------------------------------------------------
@pytest.mark.parametrize("name", list(PATHS))
def test_G2_every_path_returns_a_certified_energy_the_gate_can_reproduce(name):
    """DEFINITION OF DONE. Default arguments, real checker: every energy comes back VERIFIED,
    is the solver's own number, and its certificate verifies again from the files."""
    mh, result, directory = _certified(name)
    raws, verdicts = _assert_certificate_or_abstention(result, mh, directory)
    assert all(v.ok for v in verdicts), f"{name}: {[v.line for v in verdicts if not v.ok]}"
    # One proof can cover a whole sweep (spec section 3c), but the lowest energy has to be
    # the one it was issued for: its ceiling is that energy plus float padding, not below it.
    low = min(range(len(raws)), key=raws.__getitem__)
    assert 0.0 <= verdicts[low].hi - (raws[low] - mh.energy_offset) <= PAD_MAX, (
        f"{name}: the lowest energy {raws[low]!r} is not the one its certificate is about")


def test_G2_the_certified_paths_return_the_numbers_the_uncertified_door_returns():
    """Certification does not move the energy: QuantumKrylovSolver's certified paths agree
    with the raw Ritz values to 1e-12 Ha (two instances already differ by ~1e-15 through
    ARPACK's random start vector in the default dt)."""
    mh = _h2()
    door, _ = QuantumKrylovSolver(mh).uncertified_eigenstates(M, 2)
    for name in ("QuantumKrylovSolver.solve", "QuantumKrylovSolver.solve_excited",
                 "QuantumKrylovSolver.eigenstates", "QuantumKrylovSolver.convergence",
                 "pipeline.run_geometry"):
        final = _records(_certified(name)[1])[-1]
        assert abs(_gated(final) - door.unverified_energies[0]) <= 1e-12, name
    excited = _certified("QuantumKrylovSolver.solve_excited")[1]
    assert np.allclose(excited.energies, door.unverified_energies, rtol=0, atol=1e-12)


# -- G3 --------------------------------------------------------------------------------------
@pytest.mark.parametrize("name", list(PATHS))
def test_G3_when_the_checker_refuses_no_path_returns_the_number_as_an_energy(name, monkeypatch):
    """Every path, with the checker refusing every certificate: the result is an explicit
    abstention, its number is the one the certified run returned, and it is readable only as
    ``unverified``."""
    monkeypatch.setattr(certification, "check_certificate",
                        lambda *a, **k: (False, "ABSTAIN   refused by the gate"))
    mh, result = PATHS[name][1](None)
    raws, verdicts = _assert_certificate_or_abstention(result, mh)
    assert verdicts and not any(v.ok for v in verdicts)
    assert all(v.line.startswith("ABSTAIN") for v in verdicts)
    _, certified_result, directory = _certified(name)
    certified_raws, _ = _assert_certificate_or_abstention(certified_result, mh, directory)
    assert np.allclose(raws, certified_raws, rtol=0, atol=1e-12), (
        f"{name}: the energies depend on what the checker said")


def test_G3_the_uncertified_door_returns_abstentions_and_starts_no_checker(monkeypatch):
    """The escape hatch is visible in what it returns: abstentions that say they were not
    submitted, the number only under ``unverified``, and no checker process at all."""
    def no_checker(*a, **k):
        raise AssertionError("the uncertified door started a process")
    monkeypatch.setattr(subprocess, "run", no_checker)

    mi, mh_d = _dimer()
    doors = {
        "uncertified_eigenstates": (_h2(), QuantumKrylovSolver(_h2()).uncertified_eigenstates(M, 2)),
        "run_geometry": (_h2(), run_geometry(H2, krylov_dim=M, reference="none",
                                             uncertified_estimates_only=True)),
        "run_geometry(track_convergence=False)": (_h2(), run_geometry(
            H2, krylov_dim=M, reference="none", track_convergence=False,
            uncertified_estimates_only=True)),
        "run_from_integrals": (mh_d, run_from_integrals(
            mi.h1, mi.eri, mi.nelec, mi.e_core, krylov_dim=M, uncertified_estimates_only=True)),
    }
    for name, (mh, result) in doors.items():
        raws, verdicts = _assert_certificate_or_abstention(result, mh)
        assert not any(v.ok for v in verdicts), name
        assert all("not submitted" in v.line for v in verdicts), [v.line for v in verdicts]
    states = doors["uncertified_eigenstates"][1][1]
    assert states.shape == (2, 1 << _h2().num_qubits)        # the states are still returned


def test_G3_a_checker_that_times_out_is_an_abstention_not_a_hang(monkeypatch):
    """A hung checker is killed at CHECK_TIMEOUT_S and reported as a CRASH abstention, and
    reading the energy then raises: certificate-or-abstention holds when the checker never
    answers."""
    def hung(cmd, **kw):
        assert kw.get("timeout") == certification.CHECK_TIMEOUT_S, "the checker runs unbounded"
        raise subprocess.TimeoutExpired(cmd, kw["timeout"])
    monkeypatch.setattr(subprocess, "run", hung)
    step = QuantumKrylovSolver(_h2()).solve(M)
    assert not step.certificate.ok
    assert step.certificate.line.startswith("CRASH  timeout"), step.certificate.line
    with pytest.raises(EnergyAbstained):
        step.energy


def test_G3_a_result_cannot_be_built_without_a_verdict():
    """No energy without certificate-or-abstention, enforced where results are made."""
    from hybrid_quantum_solver.quantum_krylov_solver import ExcitedKrylovStep, KrylovStep
    from hybrid_quantum_solver.skqd import SKQDStep
    for build in (lambda c: KrylovStep(dim=1, rank=1, unverified_energy=-1.0, certificate=c),
                  lambda c: ExcitedKrylovStep(dim=1, rank=1, unverified_energies=[-1.0],
                                              certificate=c),
                  lambda c: SKQDStep(depth=1, n_dets=1, unverified_energy=-1.0, certificate=c)):
        for not_a_verdict in (None, True, "VERIFIED", ("x", "r", True, "VERIFIED", -2.0, -1.0)):
            with pytest.raises(TypeError):
                build(not_a_verdict)


# -- G4 --------------------------------------------------------------------------------------
def test_G4_the_real_checker_refuses_a_sub_ground_state_energy_and_the_read_raises():
    """Shot noise breaks the variational floor. An energy below lambda_min is refused by the
    checker itself -- its verdict, not a crash -- and reading it as an energy raises.
    Non-vacuous by construction: at least one seed has to dip below."""
    mh = _h2()
    lam = mh.ground_state_energy() - mh.energy_offset
    below = 0
    for seed in range(12):
        step = QuantumKrylovSolver(mh, noise_sigma=3e-2, seed=seed).solve(6)
        if not step.certificate.ok:
            _assert_certificate_or_abstention(step, mh)
        if step.unverified_energy - mh.energy_offset < lam:
            below += 1
            assert not step.certificate.ok and step.certificate.line.startswith("ABSTAIN"), (
                f"seed {seed}: below lambda_min yet {step.certificate.line}")
            with pytest.raises(EnergyAbstained, match="ABSTAIN"):
                step.energy
    assert below, "no seed broke the variational floor: G4 tested nothing"


def test_G4_in_a_noisy_sweep_no_energy_borrows_a_certificate_from_above():
    """One proof covers the energies AT OR ABOVE the one it was issued for, never one below.
    Noisy sweeps mix both: every step below lambda_min must come back abstained even though
    a higher step of the same call verified, and every verified step's ceiling must sit at
    or below its own energy (the detector checks that, and re-runs the checker)."""
    mh = _h2()
    lam = mh.ground_state_energy() - mh.energy_offset
    mixed = 0
    for seed in range(4):
        directory = OUT / "noisy_sweep" / f"seed{seed}"
        shutil.rmtree(directory, ignore_errors=True)
        steps = QuantumKrylovSolver(mh, noise_sigma=3e-2, seed=seed).convergence(6, directory)
        _assert_certificate_or_abstention(steps, mh, directory)
        below = [s for s in steps if s.unverified_energy - mh.energy_offset < lam]
        assert not any(s.certificate.ok for s in below), (
            f"seed {seed}: an energy below lambda_min was certified")
        mixed += bool(below) and any(s.certificate.ok for s in steps)
    assert mixed, "no sweep mixed verified and sub-ground-state steps: tested nothing"


# -- G5 --------------------------------------------------------------------------------------
def test_G5_the_detector_rejects_a_bare_energy_and_three_other_mutants(tmp_path, monkeypatch):
    """Non-vacuity, demonstrated: each way of getting an uncertified number past the caller
    turns the detector red."""
    mh = _h2()
    step = _certified("QuantumKrylovSolver.solve")[1]
    e = step.unverified_energy

    # 1. The pre-0pc result: a float called ``energy`` and no certificate.
    with pytest.raises(AssertionError, match="Verdict belongs"):
        _assert_certificate_or_abstention(
            types.SimpleNamespace(dim=M, rank=M, energy=e, certificate=None), mh)

    # 2. An abstention that also exposes the number under an innocent name.
    class Leaky:
        certificate = certification.not_submitted("mutant")
        unverified_energy = e
        estimate = e

        @property
        def energy(self):
            raise EnergyAbstained("mutant")

    with pytest.raises(AssertionError, match="readable as"):
        _assert_certificate_or_abstention(Leaky(), mh)

    # 3. A VERIFIED verdict with nothing behind it: right numbers, no certificate file.
    forged = dataclasses.replace(step, certificate=step.certificate._replace(name="forged"))
    with pytest.raises(AssertionError, match="no certificate file"):
        _assert_certificate_or_abstention(forged, mh, tmp_path)

    # 4. The mutation the spec names: a path that goes back to returning the bare number.
    monkeypatch.setattr(QuantumKrylovSolver, "solve",
                        lambda self, krylov_dim, certificate_dir=None:
                        float(self._ritz(krylov_dim)[0][0]) + self.offset)
    with pytest.raises(AssertionError, match="bare energy"):
        _assert_certificate_or_abstention(PATHS["QuantumKrylovSolver.solve"][1](tmp_path)[1], mh)
    # The pipeline does not wait for the detector: it refuses to wrap the bare number.
    with pytest.raises(TypeError, match="certificate-or-abstention"):
        run_geometry(H2, krylov_dim=M, reference="none", track_convergence=False)


def test_G3_a_sweep_with_an_abstained_step_still_returns_every_step():
    # portfolio-dry7: one abstained step must not hide the rest of the sweep. Each step is
    # gated on its own verdict; only reading the abstained step's energy raises.
    from hybrid_quantum_solver.certification import not_submitted
    from hybrid_quantum_solver.pipeline import PipelineResult
    from hybrid_quantum_solver.quantum_krylov_solver import KrylovStep

    ok = Verdict("v", "temple_inertia", True, "VERIFIED  stub", -2.0, -1.1)
    steps = [KrylovStep(1, 1, -1.10, not_submitted("stub abstention")),
             KrylovStep(2, 2, -1.13, ok)]
    result = PipelineResult(steps=steps, hf_energy=-1.1, reference_energy=None,
                            n_qubits=4, dt=0.1, shots=None)
    conv = result.convergence
    assert [(s.dim, s.rank, s.certificate.ok) for s in conv] == [(1, 1, False), (2, 2, True)]
    assert conv[1].energy == -1.13
    with pytest.raises(EnergyAbstained):
        conv[0].energy
