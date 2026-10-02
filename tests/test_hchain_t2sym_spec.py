"""Gates for SPEC_hchain_t2sym (chem-7ko): SU(2) x reflection adapted T2 for the certified H_n bracket."""

import os
import sys
from itertools import product

import numpy as np
import pytest
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import vsdp_hchain as V  # noqa: E402
from test_hchain_vsdp_bracket_spec import _ladder, _word_matrix  # noqa: E402


def _op(row, ops):
    """A block row is a word (tuple) or an integer combination of words ({word: coef})."""
    if isinstance(row, dict):
        return sum(c * _word_matrix(w, ops) for w, c in row.items())
    return _word_matrix(row, ops)


def _dag(M):
    return M.T.conj()


@pytest.fixture(scope="module")
def h4():
    """Brute-force ground state of the reflection-symmetrized H4 Hamiltonian, and every block entry."""
    h1, eri, (na, nb), e_nuc = V.sym_integrals(4)
    r = h1.shape[0]
    ops = _ladder(2 * r)
    H = sp.csr_matrix(ops[0].shape)
    for s in (0, 1):
        for i, k in product(range(r), repeat=2):
            H = H + h1[i, k] * _word_matrix(((True, s * r + i), (False, s * r + k)), ops)
    for s, t in product((0, 1), repeat=2):
        for i, j, k, l_ in product(range(r), repeat=4):
            w = ((True, s * r + i), (True, t * r + j), (False, t * r + l_), (False, s * r + k))
            H = H + 0.5 * eri[i, k, j, l_] * _word_matrix(w, ops)
    occ = lambda x, p: (x >> (2 * r - 1 - p)) & 1  # noqa: E731
    sector = [x for x in range(2 ** (2 * r))
              if sum(occ(x, p) for p in range(r)) == na and sum(occ(x, p) for p in range(r, 2 * r)) == nb]
    w_, v_ = np.linalg.eigh(H[sector][:, sector].toarray())
    psi = np.zeros(2 ** (2 * r))
    psi[sector] = v_[:, 0]
    Sp = sum(_word_matrix(((True, i), (False, r + i)), ops) for i in range(r))
    prob = V.build(h1, eri, na, nb, "sym")
    x = np.zeros(prob["tot"])
    for k in prob["names"]:
        R, n_, o = prob["rows"][k], prob["size"][k], prob["off"][k]
        mats = [_op(row, ops) for row in R]
        for a, b in product(range(n_), repeat=2):
            M = _dag(mats[a]) @ mats[b]
            if k.startswith("T2"):
                M = M + mats[b] @ _dag(mats[a])
            x[o + a + b * n_] = psi @ (M @ psi)
    return dict(prob=prob, x=x, e0=w_[0] + e_nuc, gap=w_[1] - w_[0], e_nuc=e_nuc,
                splus_norm=np.linalg.norm(Sp @ psi))


def test_g0_assumptions_hold_at_h4(h4):
    """The adaptation needs a nondegenerate singlet: S+ psi = 0 and a finite gap."""
    assert h4["splus_norm"] < 1e-10 and h4["gap"] > 1e-3


def test_g1_adapted_maps_match_brute_force(h4):
    """Every adapted T2 entry, singlet/reflection equality, objective and trace bound holds at exact psi."""
    prob, x = h4["prob"], h4["x"]
    assert {"T2h+", "T2h-", "T2q+", "T2q-"} <= set(prob["names"])
    assert np.max(np.abs(prob["A"] @ x - prob["b"])) < 1e-10
    assert abs(prob["c"] @ x + h4["e_nuc"] - h4["e0"]) < 1e-10
    for T, ks in prob["groups"].values():
        tr = sum(np.trace(x[prob["off"][k]:prob["off"][k] + prob["size"][k] ** 2].reshape(prob["size"][k], -1))
                 for k in ks)
        assert tr <= T + 1e-10
        if len(ks) == 1:
            assert abs(tr - T) < 1e-10, ks


@pytest.mark.parametrize("n", [4, 6])
def test_g2_block_sizes_shrink(n):
    """Adapted T2 = HW(S=1/2) and S=3/2 blocks, each split by parity: dims add to M=1/2 minus M=3/2."""
    r = n
    full = V.block_rows(r, True)
    prob, _ = V.hchain_problem(n, "sym")
    s = prob["size"]
    assert s["T2h+"] + s["T2h-"] == len(full["T2_+0"]) - len(full["T2_+1"])
    assert s["T2q+"] + s["T2q-"] == len(full["T2_+1"])
    assert max(s[k] for k in s if k.startswith("T2")) <= len(full["T2_+0"]) / 2.5


def test_g3_same_bound_as_full_t2_and_valid():
    """Group averaging: the adapted SDP has the same optimum as Sz-blocked T2, and certifies <= E0."""
    from pyscf import fci

    h1, eri, ne, e_nuc = V.sym_integrals(4)
    e0 = fci.direct_spin1.kernel(h1, eri, 4, ne, conv_tol=1e-12)[0] + e_nuc
    out = {}
    for t2 in (True, "sym"):
        prob = V.build(h1, eri, *ne, t2)
        val, y, _ = V.solve(prob, "SCS", eps_abs=1e-9, eps_rel=1e-9, max_iters=200000)
        out[t2] = (val + e_nuc, V.certify(prob, y, e_nuc)[0])
    assert out["sym"][1] <= e0 and out[True][1] <= e0
    assert abs(out["sym"][0] - out[True][0]) < 2e-6
    assert e0 - out["sym"][1] < 1e-5
