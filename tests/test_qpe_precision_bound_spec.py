"""
Acceptance gates G1-G4 for specs/SPEC_qpe_precision_bound.md (chem-0ke, revised chem-bu4): QPE's
precision constant -- measured as ~2.175 on H2 CAS(2,2) by SPEC_qpe_readout_laws and explicitly
left without a derivation there -- is tested here against a specific closed-form prediction,
pi*sin(theta_0) with theta_0 = arccos(E_0/lambda), PRE-REGISTERED in
specs/SPEC_qpe_precision_bound.md sec.3 (and reproducible via
scripts/spec_qpe_precision_bound.py) from the exact spectrum alone, before any QPE sweep runs.

G4 (chem-bu4) gates the EXACT nearest-dyadic-bin ceiling -- a strict per-t bound, unlike
pi*sin(theta_0) which is only its large-t limit -- on both the ordinary (recenter=False) and
recentered (recenter=True, SPEC_ft_identity_shift) walk operator, across all of t=4..14. This is
the gate that would have caught chem-177's landing finding (recentered H2 at t=4 exceeds
pi*sin(theta_0') by 2.9%) as a non-issue instead of a footnote.

Deliberately no new library code: reuses `qpe_walk_readout.run_qpe` and
`qubitization_blueprint.{build_qubit_hamiltonian,pauli_decompose,split_identity}` unmodified.
"""
import numpy as np
import pytest
from pyscf import ao2mo, gto, mcscf, scf

from qpe_walk_readout import run_qpe
from qubitization_blueprint import build_qubit_hamiltonian, pauli_decompose, split_identity

T_RANGE = list(range(4, 15))  # t = 4..14 inclusive, per the bead's pre-registered sweep

# --- Pre-registered prediction (specs/SPEC_qpe_precision_bound.md sec.3), committed BEFORE
# --- this file's sweep gates were run. Recomputed from the exact spectrum in G3 below and
# --- checked for drift against these committed numbers -- not fitted from the sweep itself.
PREREGISTERED = {
    "H2 CAS(2,2)":  dict(theta0=2.327122, bound=2.285077),
    "LiH CAS(2,2)": dict(theta0=2.295789, bound=2.351497),
    "N2 CAS(3,4)":  dict(theta0=2.273942, bound=2.396444),
}

SYSTEMS = {
    "H2 CAS(2,2)":  dict(atom="H 0 0 0; H 0 0 0.74", symmetry=None, norb=2, ne=2),
    "LiH CAS(2,2)": dict(atom="Li 0 0 0; H 0 0 1.6", symmetry=None, norb=2, ne=2),
    "N2 CAS(3,4)":  dict(atom="N 0 0 0; N 0 0 1.10", symmetry="D2h", norb=3, ne=4),
}


def _build(spec):
    kwargs = dict(atom=spec["atom"], basis="sto-3g")
    if spec["symmetry"]:
        kwargs["symmetry"] = spec["symmetry"]
    mol = gto.M(**kwargs)
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    cas = mcscf.CASCI(mf, spec["norb"], spec["ne"])
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), spec["norb"])
    return h1, eri, e_core, spec["norb"], float(cas.e_tot)


@pytest.fixture(scope="module", params=list(SYSTEMS.keys()))
def system(request):
    name = request.param
    h1, eri, e_core, norb, casci = _build(SYSTEMS[name])
    H, n = build_qubit_hamiltonian(h1, eri, norb)
    Ek, Vk = np.linalg.eigh(H)
    ground = Vk[:, 0]
    lam = sum(abs(c) for _, c in pauli_decompose(H, n))
    E0 = float(Ek[0])
    theta0 = float(np.arccos(np.clip(E0 / lam, -1, 1)))
    bound = float(np.pi * np.sin(theta0))
    return dict(name=name, h1=h1, eri=eri, e_core=e_core, norb=norb, casci=casci,
                ground=ground, lam=lam, E0=E0, theta0=theta0, bound=bound)


def _measured_ratios(system):
    ratios = []
    for t in T_RANGE:
        E_est, lam, _ps, _olap = run_qpe(
            system["h1"], system["eri"], system["norb"], system["e_core"], system["ground"], t
        )
        err = abs(E_est - system["casci"])
        ratios.append(err / (lam / 2 ** t))
    return ratios


def test_G3_preregistered_prediction_matches_exact_spectrum(system):
    """Bookkeeping gate: the committed table (specs/SPEC_qpe_precision_bound.md sec.3),
    recomputed here from the exact spectrum, must match to 1e-6 -- guards against the
    committed prediction drifting out of sync with the code that generates it."""
    expected = PREREGISTERED[system["name"]]
    assert system["theta0"] == pytest.approx(expected["theta0"], abs=1e-6)
    assert system["bound"] == pytest.approx(expected["bound"], abs=1e-6)


def test_G1_measured_ratio_never_exceeds_predicted_bound(system):
    """The kill condition: err(t)/(lambda/2^t) must never exceed pi*sin(theta_0), computed
    from the exact spectrum BEFORE this sweep ran, at any t=4..14 on any system."""
    ratios = _measured_ratios(system)
    bound = system["bound"]
    assert max(ratios) <= bound, (system["name"], bound, ratios)


def _exact_nearest_bin_ceiling(theta0, t):
    """Exact per-t ceiling on err(t)/(lambda/2^t) -- chem-bu4's G4. Unlike pi*sin(theta0) (the
    t -> infinity envelope this spec otherwise tests), this is a strict bound at EVERY t: the
    argmax decode rounds the true phase theta0 to the nearest dyadic grid point, at most pi/2^t
    away, and cos is 1-Lipschitz (mean value theorem), so the realized error can never exceed the
    exact worst case at the two neighboring grid points."""
    return 2 ** t * max(abs(np.cos(theta0 + s * np.pi / 2 ** t) - np.cos(theta0)) for s in (1, -1))


def _walk_theta0(h1, eri, norb, recenter):
    """theta0 for the walk variant run_qpe uses internally, matching its own recentering exactly
    (qubitization_blueprint.split_identity) -- so the ceiling is computed on the same spectrum
    run_qpe decodes against, not a separately-reconstructed one."""
    H, n = build_qubit_hamiltonian(h1, eri, norb)
    Ek, _Vk = np.linalg.eigh(H)
    terms = pauli_decompose(H, n)
    if recenter:
        non_identity, c_identity = split_identity(terms)
        lam = sum(abs(c) for _, c in non_identity)
    else:
        c_identity, lam = 0.0, sum(abs(c) for _, c in terms)
    E0 = float(Ek[0]) - c_identity
    return float(np.arccos(np.clip(E0 / lam, -1, 1)))


@pytest.mark.parametrize("recenter", [False, True])
def test_G4_exact_nearest_bin_bound_never_exceeded(system, recenter):
    """chem-bu4's real per-t falsifier: unlike G1's pi*sin(theta0) large-t envelope, this exact
    half-bin ceiling must never be exceeded at ANY t, on BOTH the ordinary (recenter=False) walk
    G1-G3 already cover and the recentered (recenter=True) walk SPEC_ft_identity_shift adds -- the
    landing finding that killed pi*sin(theta0) as a strict bound (recentered H2 at t=4: ratio
    2.6966 > bound 2.6212) does NOT exceed this ceiling (2.7739)."""
    theta0 = _walk_theta0(system["h1"], system["eri"], system["norb"], recenter)
    for t in T_RANGE:
        E_est, lam_r, _ps, _olap = run_qpe(
            system["h1"], system["eri"], system["norb"], system["e_core"], system["ground"], t,
            recenter=recenter,
        )
        err = abs(E_est - system["casci"])
        ratio = err / (lam_r / 2 ** t)
        ceiling = _exact_nearest_bin_ceiling(theta0, t)
        assert ratio <= ceiling + 1e-9, (system["name"], recenter, t, ratio, ceiling)


def test_G2_bound_is_tight_on_at_least_one_system():
    """The second kill condition: if pi*sin(theta_0) is a valid bound but loose by more than
    20% on EVERY system, it is a real bound but does not explain the actual mechanism. This
    check is cross-system (not parametrized) because "at least one" is a claim about the set."""
    tightness = {}
    for name, spec in SYSTEMS.items():
        h1, eri, e_core, norb, casci = _build(spec)
        H, n = build_qubit_hamiltonian(h1, eri, norb)
        Ek, Vk = np.linalg.eigh(H)
        ground = Vk[:, 0]
        lam = sum(abs(c) for _, c in pauli_decompose(H, n))
        E0 = float(Ek[0])
        theta0 = float(np.arccos(np.clip(E0 / lam, -1, 1)))
        bound = float(np.pi * np.sin(theta0))
        sysdict = dict(h1=h1, eri=eri, e_core=e_core, norb=norb, casci=casci, ground=ground,
                       lam=lam, E0=E0, theta0=theta0, bound=bound)
        max_ratio = max(_measured_ratios(sysdict))
        tightness[name] = max_ratio / bound

    assert any(t >= 0.8 for t in tightness.values()), tightness


if __name__ == "__main__":
    # The measured side of SPEC_qpe_precision_bound.md (sec.10): uv run python tests/test_qpe_precision_bound_spec.py
    for name, spec in SYSTEMS.items():
        h1, eri, e_core, norb, casci = _build(spec)
        _Ek, Vk = np.linalg.eigh(build_qubit_hamiltonian(h1, eri, norb)[0])
        ratios = _measured_ratios(dict(h1=h1, eri=eri, e_core=e_core, norb=norb, casci=casci, ground=Vk[:, 0]))
        bound = PREREGISTERED[name]["bound"]
        print(f"{name}: max_ratio(t=4..14)={max(ratios):.6f} at t={T_RANGE[int(np.argmax(ratios))]}  "
              f"bound={bound:.6f}  max/bound={max(ratios) / bound:.1%}")
