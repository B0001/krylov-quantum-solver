"""
Acceptance gates G1-G4 for specs/SPEC_certkit_return_paths.md (portfolio-cir).

Claim: every ground-state energy the Krylov solver returns can carry a certkit v1
lambda_min_enclosure certificate whose certified CEILING IS THE RETURNED NUMBER, checked by
the pinned checker out of process -- and the checker refuses one for an energy that is not a
true upper bound on the ground state. Every energy-returning method in hybrid_quantum_solver
takes ``certify=<dir>``; this gate drives each one.

Electronic frame throughout, as in test_certkit_regression_gate_spec.py: the certificate is
about the qubit matrix, so lambda_min is its GLOBAL lowest eigenvalue over every particle
sector, and mh.energy_offset is a trusted addition outside it.

Needs the pinned checker (`uv sync --inexact --extra certkit`). PySCF/qiskit, no block2.
"""
import importlib.util
import os
from pathlib import Path

import pytest

if importlib.util.find_spec("certkit") is None:
    # Same guard as the regression gate: optional locally, and never a silent pass in CI.
    if os.environ.get("CI"):
        raise RuntimeError(
            "certkit is not installed: the return-path certificate gate cannot run. "
            'Install the pinned checker with `uv sync --inexact --extra certkit`.'
        )
    pytest.skip('certkit extra not installed (uv sync --inexact --extra certkit)',
                allow_module_level=True)

import certkit_bridge  # noqa: E402
from hybrid_quantum_solver import (  # noqa: E402
    HardwareKrylovSolver,
    QuantumKrylovSolver,
    TrotterKrylovSolver,
    build_hamiltonian_from_integrals,
    build_molecular_hamiltonian,
    hubbard_chain_integrals,
    run_from_integrals,
    run_geometry,
)
from hybrid_quantum_solver.skqd import SampleKrylovSolver  # noqa: E402

H2 = "H 0 0 0; H 0 0 0.74"
M = 3
# The certified ceiling is the returned energy plus pad_claim's float padding, which is
# 1e-9 relative here; anything past 1e-7 Ha means the claim is no longer about that number.
PAD_MAX = 1e-7
# Written where the sibling gate writes, so a CI failure can upload the certificates.
OUT = Path(__file__).resolve().parent.parent / "certkit_out" / "return_paths"

_CACHE = {}


def _h2():
    if "h2" not in _CACHE:
        _CACHE["h2"] = build_molecular_hamiltonian(atom=H2)
    return _CACHE["h2"]


def _dimer():
    """Hubbard dimer (U=4), for run_from_integrals. e_core=0 and no chemical potential, so
    the global lambda_min lies in another particle sector, 0.17 Ha below the (1,1) energy
    the solver returns: a true but loose upper bound, and the Temple route does not apply."""
    if "dimer" not in _CACHE:
        mi = hubbard_chain_integrals(2, 4.0)
        mh = build_hamiltonian_from_integrals(mi.h1, mi.eri, num_particles=mi.nelec,
                                              energy_offset=mi.e_core)
        _CACHE["dimer"] = (mi, mh)
    return _CACHE["dimer"]


def _paths(certify):
    """path -> (mh, [(returned total energy, certificate)]) for every energy-returning path.

    certify=None runs the same calls uncertified, so G1 can show the flag changes nothing.
    """
    def sub(name):
        return None if certify is None else certify / name

    mh = _h2()
    p = {}
    s = QuantumKrylovSolver(mh).solve(M, sub("qks_solve"))
    p["QuantumKrylovSolver.solve"] = [(s.energy, s.certificate)]
    p["QuantumKrylovSolver.convergence"] = [
        (s.energy, s.certificate)
        for s in QuantumKrylovSolver(mh).convergence(M, certify=sub("qks_convergence"))]
    s = QuantumKrylovSolver(mh).solve_excited(M, n_states=2, certify=sub("qks_excited"))
    p["QuantumKrylovSolver.solve_excited"] = [(s.energies[0], s.certificate)]
    r = QuantumKrylovSolver(mh).eigenstates(M, n_states=2, certify=sub("qks_eigenstates"))
    p["QuantumKrylovSolver.eigenstates"] = [(r[0][0], r[2] if certify else None)]
    p["TrotterKrylovSolver.convergence"] = [
        (s.energy, s.certificate) for s in TrotterKrylovSolver(mh, trotter_reps=2).convergence(
            M, certify=sub("trotter"))]
    s = HardwareKrylovSolver(mh).solve(M, sub("hardware_solve"))
    p["HardwareKrylovSolver.solve"] = [(s.energy, s.certificate)]
    p["HardwareKrylovSolver.convergence"] = [
        (s.energy, s.certificate)
        for s in HardwareKrylovSolver(mh).convergence(M, certify=sub("hardware_convergence"))]
    skqd = SampleKrylovSolver(mh, depth=4, n_shots=50, seed=0)
    s = skqd.solve(sub("skqd_solve"))
    p["SampleKrylovSolver.solve"] = [(s.energy, s.certificate)]
    p["SampleKrylovSolver.convergence"] = [
        (s.energy, s.certificate) for s in SampleKrylovSolver(mh, depth=4, seed=0).convergence(
            (5, 50), certify=sub("skqd_convergence"))]
    r = run_geometry(H2, krylov_dim=M, certify=sub("run_geometry"))
    p["run_geometry"] = [(e, c) for (_, _, e), c in
                         zip(r.convergence, r.certificates or [None] * len(r.convergence))]
    r = run_geometry(H2, krylov_dim=M, track_convergence=False, certify=sub("run_geometry_final"))
    p["run_geometry(track_convergence=False)"] = [
        (r.computed_energy, r.certificates[-1] if certify else None)]
    out = {name: (mh, pairs) for name, pairs in p.items()}

    mi, mh_d = _dimer()
    r = run_from_integrals(mi.h1, mi.eri, mi.nelec, mi.e_core, krylov_dim=M,
                           certify=sub("run_from_integrals"))
    out["run_from_integrals"] = (mh_d, [(e, c) for (_, _, e), c in
                                        zip(r.convergence,
                                            r.certificates or [None] * len(r.convergence))])
    return out


def _certified():
    if "certified" not in _CACHE:
        _CACHE["certified"] = _paths(OUT)
    return _CACHE["certified"]


def _lam(mh):
    return mh.ground_state_energy() - mh.energy_offset


def test_G1_every_path_returns_a_verified_certificate_and_the_same_energy():
    """Every energy every path returns carries a certificate the checker VERIFIED, by the
    route pinned for it, and asking for one does not change the energy returned."""
    plain = _paths(None)
    for name, (mh, pairs) in _certified().items():
        # Not bitwise: the default dt comes from ARPACK with a random start vector, so two
        # solver instances already differ by ~1e-15 Ha before certify= enters at all.
        drift = max(abs(a - b) for (a, _), (b, _) in zip(pairs, plain[name][1]))
        assert len(pairs) == len(plain[name][1]) and drift <= 1e-12, (
            f"{name}: certify= changed the returned energies (by {drift:.1e} Ha)")
        # H2 has 16 states, inside certkit's DENSE_LIMIT; the dimer's sector energy sits too
        # high for a finite Temple bound on the global spectrum (see _dimer).
        route = "gershgorin_rayleigh" if name == "run_from_integrals" else "temple_inertia"
        for e, v in pairs:
            assert v is not None, f"{name}: energy {e!r} came back without a certificate"
            # The checker's own word must agree with its exit code (a crash also exits 1).
            assert v.ok and v.line.startswith("VERIFIED"), f"{name}: {v.line}"
            assert v.rule == route, f"{name}: certified via {v.rule}, pinned {route}"


def test_G2_the_certificate_is_about_the_returned_energy_and_contains_the_ground_state():
    """DEFINITION OF DONE. lambda_min lies in [lo, hi], and hi is the returned energy: the
    checker has proven that number an upper bound on the ground state."""
    for name, (mh, pairs) in _certified().items():
        lam = _lam(mh)
        for e_total, v in pairs:
            e = e_total - mh.energy_offset
            assert v.lo <= lam <= v.hi, (
                f"{name}: exact {lam!r} escaped the certified [{v.lo!r}, {v.hi!r}]")
            assert 0.0 <= v.hi - e <= PAD_MAX, (
                f"{name}: certified ceiling {v.hi!r} is not the returned energy {e!r}")


def test_G3_an_energy_below_the_ground_state_is_refused():
    """The gate can fail. Shot noise on H and S breaks the variational floor; every noisy
    energy that lands below lambda_min must ABSTAIN -- the checker's verdict, not a crash --
    and the ones it does verify must still contain lambda_min. Non-vacuous by construction:
    at least one of the seeds has to dip below."""
    mh = _h2()
    lam = _lam(mh)
    below = 0
    for seed in range(12):
        s = QuantumKrylovSolver(mh, noise_sigma=3e-2, seed=seed).solve(
            6, certify=OUT / "noise" / f"seed{seed}")
        e, v = s.energy - mh.energy_offset, s.certificate
        if e < lam:
            below += 1
            assert not v.ok and v.line.startswith("ABSTAIN"), (
                f"seed {seed}: energy {e!r} is below lambda_min {lam!r} yet got {v.line}")
        elif v.ok:
            assert v.lo <= lam <= v.hi, f"seed {seed}: {lam!r} outside [{v.lo!r}, {v.hi!r}]"
    assert below, "no seed broke the variational floor: G3 tested nothing"


def test_G4_beyond_the_dense_limit_the_gershgorin_route_certifies():
    """N2 CAS(6,6) is 4096 states, past certkit's DENSE_LIMIT, so there is no inertia count
    and the only route is gershgorin_rayleigh: VERIFIED, with a floor 2.7 Ha down."""
    mh = build_molecular_hamiltonian(**certkit_bridge.CASES["N2"])
    s = QuantumKrylovSolver(mh).solve(8, certify=OUT / "N2")
    e, v = s.energy - mh.energy_offset, s.certificate
    # Pinned, and re-derived, in test_certkit_regression_gate_spec.py (G3); not recomputed
    # here because the dense 4096 x 4096 diagonalisation is most of that gate's runtime.
    lam = -11.392972038179622
    assert v.ok and v.line.startswith("VERIFIED"), v.line
    assert v.rule == "gershgorin_rayleigh"
    assert v.lo <= lam <= v.hi
    assert 0.0 <= v.hi - e <= PAD_MAX
