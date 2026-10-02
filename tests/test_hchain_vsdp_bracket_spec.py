"""Gates for SPEC_hchain_vsdp_bracket (chem-3z8): certified PQG lower bounds on H_n chains."""

import os
import sys
from itertools import product

import numpy as np
import pytest
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import vsdp_hchain as V  # noqa: E402
from benchmark_hchain_tdl import integrals  # noqa: E402


def _ladder(m):
    low = sp.csr_matrix(np.array([[0.0, 1.0], [0.0, 0.0]]))
    Z, Id = sp.diags([1.0, -1.0]), sp.identity(2)
    ops = []
    for p in range(m):
        op = sp.identity(1)
        for M in [Z] * p + [low] + [Id] * (m - p - 1):
            op = sp.kron(op, M, format="csr")
        ops.append(op)
    return ops


def _word_matrix(word, ops):
    M = sp.identity(ops[0].shape[0], format="csr")
    for d, p in word:
        M = M @ (ops[p].T if d else ops[p])
    return M


def _exact(n, t2=False):
    """Brute-force ground state in the (na, nb) sector and the full vector of block entries."""
    h1, eri, (na, nb), e_nuc, _ = integrals(n, localize=True)
    r = h1.shape[0]
    ops = _ladder(2 * r)
    H = sp.csr_matrix(ops[0].shape)
    for s in (0, 1):
        for i, k in product(range(r), repeat=2):
            H = H + h1[i, k] * _word_matrix(((True, s * r + i), (False, s * r + k)), ops)
    for s, t in product((0, 1), repeat=2):
        for i, j, k, l_ in product(range(r), repeat=4):
            if eri[i, k, j, l_]:
                w = ((True, s * r + i), (True, t * r + j), (False, t * r + l_), (False, s * r + k))
                H = H + 0.5 * eri[i, k, j, l_] * _word_matrix(w, ops)
    dim = 2 ** (2 * r)
    # basis index bit (2r-1-p) <-> mode p (kron order)
    occ = lambda x, p: (x >> (2 * r - 1 - p)) & 1  # noqa: E731
    sector = [x for x in range(dim)
              if sum(occ(x, p) for p in range(r)) == na and sum(occ(x, p) for p in range(r, 2 * r)) == nb]
    Hs = H[sector][:, sector].toarray()
    w_, v_ = np.linalg.eigh(Hs)
    psi = np.zeros(dim)
    psi[sector] = v_[:, 0]
    prob = V.build(h1, eri, na, nb, t2)
    x = np.zeros(prob["tot"])
    for k in prob["names"]:
        R, n_, o = prob["rows"][k], prob["size"][k], prob["off"][k]
        for a in range(n_):
            for b in range(n_):
                M = _word_matrix(V.dag(R[a]) + R[b], ops)
                if k.startswith("T2"):
                    M = M + _word_matrix(R[b] + V.dag(R[a]), ops)
                x[o + a + b * n_] = psi @ (M @ psi)
    return prob, x, w_[0] + e_nuc, e_nuc


@pytest.fixture(scope="module")
def h4():
    return _exact(4)


@pytest.mark.parametrize("t2", [False, True])
def test_g1_linear_maps_match_brute_force(t2):
    """Every P/Q/G(/T2) entry, contraction and trace identity holds at the exact H4 ground state."""
    prob, x, e0, e_nuc = _exact(4, t2)
    assert np.max(np.abs(prob["A"] @ x - prob["b"])) < 1e-10
    assert abs(prob["c"] @ x + e_nuc - e0) < 1e-10
    for k in prob["names"]:
        n, o = prob["size"][k], prob["off"][k]
        assert abs(np.trace(x[o:o + n * n].reshape(n, n)) - prob["traces"][k]) < 1e-10, k


def test_g2_h2_pqg_is_exact():
    """N = 2: the D condition is the exact ensemble N-representability condition."""
    _, _, e0, _ = _exact(2)
    prob, e_nuc = V.hchain_problem(2)
    _, y, _ = V.solve(prob)
    lb, _ = V.certify(prob, y, e_nuc)
    assert lb <= e0 and e0 - lb < 1e-6


@pytest.mark.parametrize("n", [4, 6])
def test_g3_certified_bound_below_fci_and_tight(n):
    """The certified bound is below exact E0, and within 1e-5 Ha of the solver's own value."""
    from pyscf import fci

    h1, eri, ne, e_nuc, _ = integrals(n, localize=True)
    e0 = fci.direct_spin1.kernel(h1, eri, n, ne, conv_tol=1e-12)[0] + e_nuc
    prob, e_nuc = V.hchain_problem(n)
    val, y, _ = V.solve(prob)
    lb, _ = V.certify(prob, y, e_nuc)
    assert lb <= e0
    assert val + e_nuc - lb < 1e-5


def test_g4_rigor_survives_garbage_duals(h4):
    """Adversarial: any y (zeros, noise, a perturbed optimum) still gives a bound <= E0."""
    prob, _, e0, e_nuc = h4
    _, y, _ = V.solve(prob)
    rng = np.random.default_rng(0)
    for cand in (np.zeros_like(y), rng.normal(size=y.size), y + 1e-3 * rng.normal(size=y.size), -y):
        lb, _ = V.rigorous_lower_bound(prob, cand, e_nuc)
        assert lb <= e0


def test_g5_t2_tightens_and_stays_valid():
    """PQG+T2 is a subset of PQG: its certified bound is higher, and still below E0."""
    from pyscf import fci

    h1, eri, ne, e_nuc, _ = integrals(4, localize=True)
    e0 = fci.direct_spin1.kernel(h1, eri, 4, ne, conv_tol=1e-12)[0] + e_nuc
    lbs = []
    for t2 in (False, True):
        prob, _ = V.hchain_problem(4, t2)
        lbs.append(V.certify(prob, V.solve(prob)[1], e_nuc)[0])
    assert lbs[0] < lbs[1] <= e0
