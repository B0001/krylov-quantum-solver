"""
Acceptance gates G1-G4 for specs/SPEC_certkit_return_paths.md (portfolio-cir).

Claim: every ground-state energy the Krylov solver returns carries a certkit v1
lambda_min_enclosure certificate that proves THE RETURNED NUMBER an upper bound on the ground
state, checked by the pinned checker out of process -- and the checker refuses one for an
energy that is not. Every energy-returning method in hybrid_quantum_solver certifies on every
call (portfolio-0pc; that nothing can skip it is tests/test_certificate_or_abstention_spec.py).
This gate drives each path and checks what its certificates say.

Electronic frame throughout, as in test_certkit_regression_gate_spec.py: the certificate is
about the qubit matrix, so lambda_min is its GLOBAL lowest eigenvalue over every particle
sector, and mh.energy_offset is a trusted addition outside it.

The pinned checker is a core dependency, so there is no skip guard here: a missing certkit is
an import error. PySCF/qiskit, no block2.
"""
from pathlib import Path

import certkit_bridge
from hybrid_quantum_solver import (
    HardwareKrylovSolver,
    QuantumKrylovSolver,
    TrotterKrylovSolver,
    build_hamiltonian_from_integrals,
    build_molecular_hamiltonian,
    hubbard_chain_integrals,
    run_from_integrals,
    run_geometry,
)
from hybrid_quantum_solver.skqd import SampleKrylovSolver

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


def _paths(out):
    """path -> (mh, [(returned total energy, certificate)]) for every energy-returning path.

    ``.energy`` is the certified reading: it raises unless the checker VERIFIED, so a path
    that abstained fails here before any assertion below.
    """
    mh = _h2()
    p = {}
    s = QuantumKrylovSolver(mh).solve(M, out / "qks_solve")
    p["QuantumKrylovSolver.solve"] = [(s.energy, s.certificate)]
    p["QuantumKrylovSolver.convergence"] = [
        (s.energy, s.certificate)
        for s in QuantumKrylovSolver(mh).convergence(M, out / "qks_convergence")]
    s = QuantumKrylovSolver(mh).solve_excited(M, n_states=2, certificate_dir=out / "qks_excited")
    p["QuantumKrylovSolver.solve_excited"] = [(s.energies[0], s.certificate)]
    s, _ = QuantumKrylovSolver(mh).eigenstates(M, n_states=2,
                                               certificate_dir=out / "qks_eigenstates")
    p["QuantumKrylovSolver.eigenstates"] = [(s.energies[0], s.certificate)]
    p["TrotterKrylovSolver.convergence"] = [
        (s.energy, s.certificate) for s in TrotterKrylovSolver(mh, trotter_reps=2).convergence(
            M, out / "trotter")]
    s = HardwareKrylovSolver(mh).solve(M, out / "hardware_solve")
    p["HardwareKrylovSolver.solve"] = [(s.energy, s.certificate)]
    p["HardwareKrylovSolver.convergence"] = [
        (s.energy, s.certificate)
        for s in HardwareKrylovSolver(mh).convergence(M, out / "hardware_convergence")]
    skqd = SampleKrylovSolver(mh, depth=4, n_shots=50, seed=0)
    s = skqd.solve(out / "skqd_solve")
    p["SampleKrylovSolver.solve"] = [(s.energy, s.certificate)]
    p["SampleKrylovSolver.convergence"] = [
        (s.energy, s.certificate) for s in SampleKrylovSolver(mh, depth=4, seed=0).convergence(
            (5, 50), out / "skqd_convergence")]
    r = run_geometry(H2, krylov_dim=M, certificate_dir=out / "run_geometry")
    p["run_geometry"] = [(s.energy, s.certificate) for s in r.steps]
    r = run_geometry(H2, krylov_dim=M, track_convergence=False,
                     certificate_dir=out / "run_geometry_final")
    p["run_geometry(track_convergence=False)"] = [(r.computed_energy, r.certificate)]
    paths = {name: (mh, pairs) for name, pairs in p.items()}

    mi, mh_d = _dimer()
    r = run_from_integrals(mi.h1, mi.eri, mi.nelec, mi.e_core, krylov_dim=M,
                           certificate_dir=out / "run_from_integrals")
    paths["run_from_integrals"] = (mh_d, [(s.energy, s.certificate) for s in r.steps])
    return paths


def _certified():
    if "certified" not in _CACHE:
        _CACHE["certified"] = _paths(OUT)
    return _CACHE["certified"]


def _lam(mh):
    return mh.ground_state_energy() - mh.energy_offset


def test_G1_every_path_returns_a_verified_certificate_by_the_pinned_route():
    """Every energy every path returns carries a certificate the checker VERIFIED, by the
    route pinned for it. (That certifying does not move the energy is gated, against the
    uncertified door and against a refusing checker, in test_certificate_or_abstention_spec.)"""
    for name, (mh, pairs) in _certified().items():
        # H2 has 16 states, inside certkit's DENSE_LIMIT; the dimer's sector energy sits too
        # high for a finite Temple bound on the global spectrum (see _dimer).
        route = "gershgorin_rayleigh" if name == "run_from_integrals" else "temple_inertia"
        for e, v in pairs:
            assert v is not None, f"{name}: energy {e!r} came back without a certificate"
            # The checker's own word must agree with its exit code (a crash also exits 1).
            assert v.ok and v.line.startswith("VERIFIED"), f"{name}: {v.line}"
            assert v.rule == route, f"{name}: certified via {v.rule}, pinned {route}"


def test_G2_the_certificate_is_about_the_returned_energy_and_contains_the_ground_state():
    """DEFINITION OF DONE. lambda_min lies in [lo, hi], and hi is not above the returned
    energy: the checker has proven that number an upper bound on the ground state. A sweep
    shares one proof (portfolio-0pc), issued for its lowest energy, so there the ceiling IS
    the returned number; for the steps above it the ceiling is lower, which proves more."""
    for name, (mh, pairs) in _certified().items():
        lam = _lam(mh)
        for e_total, v in pairs:
            e = e_total - mh.energy_offset
            assert v.lo <= lam <= v.hi, (
                f"{name}: exact {lam!r} escaped the certified [{v.lo!r}, {v.hi!r}]")
            assert v.hi - e <= PAD_MAX, (
                f"{name}: certified ceiling {v.hi!r} is above the returned energy {e!r}")
        e_low, v_low = min(pairs, key=lambda pair: pair[0])
        assert 0.0 <= v_low.hi - (e_low - mh.energy_offset) <= PAD_MAX, (
            f"{name}: the lowest energy {e_low!r} is not the one its certificate is about")


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
            6, OUT / "noise" / f"seed{seed}")
        # The raw number: reading ``.energy`` raises on exactly the seeds this gate is about.
        e, v = s.unverified_energy - mh.energy_offset, s.certificate
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
    s = QuantumKrylovSolver(mh).solve(8, OUT / "N2")
    e, v = s.energy - mh.energy_offset, s.certificate
    # Pinned, and re-derived, in test_certkit_regression_gate_spec.py (G3); not recomputed
    # here because the dense 4096 x 4096 diagonalisation is most of that gate's runtime.
    lam = -11.392972038179622
    assert v.ok and v.line.startswith("VERIFIED"), v.line
    assert v.rule == "gershgorin_rayleigh"
    assert v.lo <= lam <= v.hi
    assert 0.0 <= v.hi - e <= PAD_MAX
