"""
Acceptance gates G1-G4 for specs/SPEC_odmd_spectral.md (ODMD spectroscopy).

Test-first: ``odmd_spectral`` does not exist yet, so this file is RED until the spec is
implemented. Claim: the survival amplitude of a particle-removed/added reference a_i|ref> is the
Lehmann representation in ODMD-ready form -- DMD poles are the (N-/+1)-sector eigenvalues and the
Vandermonde amplitudes are the spectral weights -- giving the photoemission spectrum A(omega)
from the same 1-D signals the stack already measures. Falsified line by line against exact
Lehmann (dense diagonalization, computed independently here), cross-pinned to the capstone's
842.44 meV Nb3I8 gap, and extended to the damping-immunity finding: uniform damping enters
|lambda| only, so INTENSITIES survive along with energies.

Degenerate lines aggregate (their weights add -- physically correct for A(omega)). PySCF/qiskit,
no block2; `make gates` runs it in its own process.
"""
import numpy as np

from hybrid_quantum_solver.model_hamiltonians import ModelIntegrals
from hybrid_quantum_solver.molecular_hamiltonian import build_molecular_hamiltonian
from nb3x8_device_gap import exact_gap
from nb3x8_gaps import NB3X8_LT_BULK, dimer_cluster_integrals
from odmd_spectral import (
    greens_function_lines,
    ladder_operator,
    lines_from_signal,
    photoemission_gap,
    reference_signal,
)

_CACHE = {}


def _h4():
    if "h4" not in _CACHE:
        mh = build_molecular_hamiltonian(atom="H 0 0 0; H 0 0 0.9; H 0 0 1.8; H 0 0 2.7")
        w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
        psi_hf = np.asarray(mh.hf_state().data, dtype=complex)
        pops = np.abs(V.conj().T @ psi_hf) ** 2
        psi0 = V[:, int(np.flatnonzero(pops > 1e-8)[0])]     # exact reachable ground (N=4)
        _CACHE["h4"] = (mh, w, V, psi_hf, psi0)
    return _CACHE["h4"]


def _aggregate(energies, weights, floor):
    """Sum weights of (numerically) degenerate lines; keep those above ``floor``."""
    out = {}
    for e, wt in zip(energies, weights):
        key = round(float(e), 10)
        out[key] = out.get(key, 0.0) + float(wt)
    return [(e, wt) for e, wt in sorted(out.items()) if wt > floor]


def _exact_lines(w, V, psi_raw, floor=1e-3):
    nrm2 = float(np.real(np.vdot(psi_raw, psi_raw)))
    pops = np.abs(V.conj().T @ (psi_raw / np.sqrt(nrm2))) ** 2 * nrm2
    return _aggregate(w, pops, floor)


def _odmd_lines(mh, psi_raw, n=32, mod_window=0.2, damp=None):
    s, tau, mu, nrm2 = reference_signal(mh, psi_raw, n)
    if damp is not None:
        s = s * damp ** np.arange(n)
    return lines_from_signal(s, tau, mu, nrm2, mod_window=mod_window)


def test_G1_lehmann_line_by_line():
    """DEFINITION OF DONE: every visible exact removal line (aggregated weight > 1e-3) of every
    H4 spin-orbital reference is matched by an ODMD pole < 1e-8 Ha with weight error < 1e-4."""
    mh, w, V, psi_hf, _ = _h4()
    n_so = 2 * mh.num_spatial_orbitals
    checked = 0
    for i in range(n_so):
        raw = ladder_operator("-", i, n_so) @ psi_hf
        if np.real(np.vdot(raw, raw)) < 1e-8:
            continue
        poles, wts = _odmd_lines(mh, raw)
        for e_ex, w_ex in _exact_lines(w, V, raw):
            j = int(np.argmin(np.abs(poles - e_ex)))
            assert abs(poles[j] - e_ex) < 1e-8, (i, e_ex, poles[j] - e_ex)
            assert abs(wts[j] - w_ex) < 1e-4, (i, e_ex, wts[j] - w_ex)
            checked += 1
    assert checked >= 8, checked                          # the gate actually saw lines


def test_G2_reference_honesty():
    """Exact-ground-state reference reproduces TRUE Lehmann weights < 1e-4; the HF reference
    deviates by a real but bounded amount on H4 (max line error in (0.005, 0.05))."""
    mh, w, V, psi_hf, psi0 = _h4()
    n_so = 2 * mh.num_spatial_orbitals
    a1 = ladder_operator("-", 1, n_so)
    poles0, wts0 = _odmd_lines(mh, a1 @ psi0)
    for e_ex, w_ex in _exact_lines(w, V, a1 @ psi0):
        j = int(np.argmin(np.abs(poles0 - e_ex)))
        assert abs(wts0[j] - w_ex) < 1e-4, (e_ex, wts0[j] - w_ex)
    poles_hf, wts_hf = _odmd_lines(mh, a1 @ psi_hf)
    true_lines = dict(_exact_lines(w, V, a1 @ psi0, floor=1e-3))
    devs = [abs(wts_hf[int(np.argmin(np.abs(poles_hf - e)))] - wt)
            for e, wt in true_lines.items()]
    assert 0.005 < max(devs) < 0.05, max(devs)


def test_G3_material_spectrum_cross_pinned():
    """Nb3I8 dimer A(omega) (exact N=2 reference): both Hubbard bands present (>= 2 distinct
    lines per side) and min(omega+) - max(omega-) == the capstone charge gap < 0.01 meV."""
    p = NB3X8_LT_BULK["Nb3I8"]
    base = dimer_cluster_integrals(**p)
    mh = ModelIntegrals(base.h1, base.eri, 0.0, (1, 1), 2).to_hamiltonian()
    w, V = np.linalg.eigh(mh.qubit_hamiltonian.to_matrix())
    psi_hf = np.asarray(mh.hf_state().data, dtype=complex)
    pops = np.abs(V.conj().T @ psi_hf) ** 2
    psi0 = V[:, int(np.flatnonzero(pops > 1e-8)[0])]
    rem = greens_function_lines(mh, "-", reference=psi0)
    add = greens_function_lines(mh, "+", reference=psi0)
    assert len({round(line.omega, 2) for line in rem}) >= 2, rem
    assert len({round(line.omega, 2) for line in add}) >= 2, add
    gap = photoemission_gap(mh, reference=psi0)
    assert abs(gap - exact_gap(**p)) < 0.01, (gap, exact_gap(**p))


def test_G4_intensities_are_damping_immune():
    """s -> 0.7^k s (30% loss per step): every pole moves < 1e-6 Ha and every WEIGHT < 1e-6 --
    uniform damping enters |lambda| only; the device immunity extends to the whole spectrum."""
    mh, w, V, psi_hf, _ = _h4()
    raw = ladder_operator("-", 0, 2 * mh.num_spatial_orbitals) @ psi_hf
    poles, wts = _odmd_lines(mh, raw)
    poles_d, wts_d = _odmd_lines(mh, raw, mod_window=2.0, damp=0.7)
    for p, wt in zip(poles, wts):
        j = int(np.argmin(np.abs(poles_d - p)))
        assert abs(poles_d[j] - p) < 1e-6, (p, poles_d[j] - p)
        assert abs(wts_d[j] - wt) < 1e-6, (p, wts_d[j] - wt)
