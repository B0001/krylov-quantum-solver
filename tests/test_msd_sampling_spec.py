"""
Acceptance gates G1-G4 for specs/SPEC_msd_sampling.md (mirror subspace diagonalization).

Test-first: ``msd`` does not exist yet, so this file is RED until the spec is implemented. MSD
estimates the projected Hamiltonian matrix from overlap measurements at symmetrically shifted
timesteps (a central finite-difference of S(t) = <phi_0|e^{-iHt}|phi_0>) instead of measuring each
Pauli term, so its sampling variance scales with the stencil 1-norm fd1 = ||w||_1/delta rather than
the Hamiltonian 1-norm lambda. With an energy-level shift + a high-order stencil, fd1 << lambda when
lambda/W is large enough -> a real sampling-cost reduction. We validate the construction and the
advantage against the exact (noiseless) KQD result and FCI.

Errors use the MEDIAN over noise seeds (robust to the occasional GEVP blow-up). PySCF/qiskit, no
block2; `make gates` runs it in its own process. (N2 CAS(6,6) builds are dense-diagonalized, so this
gate is on the slower side.)
"""
import numpy as np

from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from msd import build_msd_problem, sample_ground_energy

_CACHE = {}


def _problem(key, atom, order, delta, **kw):
    if key not in _CACHE:
        mh = build_molecular_hamiltonian(atom=atom, **kw)
        _CACHE[key] = (mh, build_msd_problem(mh, n=8, order=order, delta=delta))
    return _CACHE[key]


def _median_err(prob, shots, method, seeds=120):
    e = [abs(sample_ground_energy(prob, shots, method, s) - prob.ref) for s in range(seeds)]
    return float(np.median(e))


N2 = dict(atom="N 0 0 0; N 0 0 1.1", active_electrons=6, active_orbitals=6)
H2 = dict(atom="H 0 0 0; H 0 0 0.74")


def test_G1_construction_correct():
    """Noiseless MSD (energy shift + order-8 stencil) reproduces KQD, and KQD ~ FCI."""
    mh, prob = _problem("n2_o8", N2["atom"], 8, 0.38,
                        active_electrons=6, active_orbitals=6)
    assert prob.msd_bias < 1.6e-3, prob.msd_bias                       # finite-difference bias
    assert abs(prob.ref + prob.offset - mh.ground_state_energy()) < 2e-3  # KQD Krylov ~ FCI


def test_G2_sampling_advantage():
    """DEFINITION OF DONE: on N2 CAS(6,6) the MSD median error beats KQD at matched shots."""
    _, prob = _problem("n2_o8", N2["atom"], 8, 0.38, active_electrons=6, active_orbitals=6)
    mk = _median_err(prob, 100_000, "kqd")
    mm = _median_err(prob, 100_000, "msd")
    assert mm < mk, (mm, mk)
    assert mk / mm > 1.5, (mk / mm)                                   # measured ~3.2x


def test_G3_boundary_no_advantage_at_small_lambda():
    """THE BOUNDARY: H2 (lambda/W ~ 1) gives fd1 > lambda, so MSD is NOT better -- advantage needs
    a large enough lambda/W that a high-order stencil drives fd1 below lambda."""
    _, prob = _problem("h2_o8", H2["atom"], 8, 0.6)
    assert prob.fd1 > prob.lam, (prob.fd1, prob.lam)                  # no-advantage regime
    mk = _median_err(prob, 100_000, "kqd")
    mm = _median_err(prob, 100_000, "msd")
    assert mm > mk, (mm, mk)                                          # MSD worse here


def test_G4_high_order_stencil_shrinks_noise():
    """MECHANISM: a higher-order stencil holds the bias at a larger delta, shrinking fd1 below
    lambda -- that is what turns on the advantage."""
    _, p4 = _problem("n2_o4", N2["atom"], 4, 0.14, active_electrons=6, active_orbitals=6)
    _, p8 = _problem("n2_o8", N2["atom"], 8, 0.38, active_electrons=6, active_orbitals=6)
    assert p4.msd_bias < 1.6e-3 and p8.msd_bias < 1.6e-3, (p4.msd_bias, p8.msd_bias)
    assert p8.fd1 < p4.fd1, (p8.fd1, p4.fd1)                          # order 8 < order 4
    assert p8.fd1 < p8.lam, (p8.fd1, p8.lam)                          # into the advantage regime
