"""Certified lower bounds on H_n chain ground-state energies: PQG 2-RDM SDP + rigorous dual (chem-3z8).

The SDP solver is an untrusted producer. `rigorous_lower_bound` turns ANY dual vector y into a bound
on the true E0 using only weak duality at the true RDM X*:
    E0 = <C, X*> = b.y + sum_j <Z_j, X*_j> >= b.y + sum_j tr(X*_j) * min(0, lambda_min(Z_j)),
with Z_j = C_j - A_j^T y, tr(X*_j) exact for any fixed-(Na, Nb) state, and lambda_min(Z_j) bounded
from below by a floating-point Cholesky (Higham, Accuracy & Stability 2nd ed., Thm 10.3) plus explicit
entrywise rounding bounds -- no directed rounding needed. The bound is for the Hamiltonian defined by
the float64 integrals as given (as in VSDP, Chaykin et al. JCTC 2020).

Every block is a Gram matrix M[a,b] = <O_a^dag O_b>; its entries are reduced to 1-/2-RDM parameters by
`normal_order`, so no P/Q/G formula is hand-typed. tests/test_hchain_vsdp_bracket_spec.py checks the
whole map against brute-force Fock-space expectation values.
"""

import argparse
import json
import math
import os
from itertools import combinations

import numpy as np
import scipy.sparse as sp

U = 2.0**-53


def gamma(k):
    return k * U / (1 - k * U)


# ---- second quantization. A word is a tuple of (is_creator, spin_orbital); alpha = i, beta = r + i.

def _sort_sign(xs, desc):
    xs, sign = list(xs), 1
    for i in range(len(xs)):
        for j in range(len(xs) - 1 - i):
            if (xs[j] < xs[j + 1]) if desc else (xs[j] > xs[j + 1]):
                xs[j], xs[j + 1] = xs[j + 1], xs[j]
                sign = -sign
    return tuple(xs), sign


def canon(word):
    """Normal-ordered word -> ((creators asc), (annihilators desc)), sign; sign 0 if it vanishes."""
    C = [p for d, p in word if d]
    A = [p for d, p in word if not d]
    if len(set(C)) < len(C) or len(set(A)) < len(A):
        return None, 0
    C, s1 = _sort_sign(C, False)
    A, s2 = _sort_sign(A, True)
    return (C, A), s1 * s2


def normal_order(word):
    """{canonical key: integer coef} with word == sum coef * key, via {a_p, a_q^dag} = delta_pq."""
    out, stack = {}, [(tuple(word), 1)]
    while stack:
        w, c = stack.pop()
        for k in range(len(w) - 1):
            if not w[k][0] and w[k + 1][0]:
                stack.append((w[:k] + (w[k + 1], w[k]) + w[k + 2:], -c))
                if w[k][1] == w[k + 1][1]:
                    stack.append((w[:k] + w[k + 2:], c))
                break
        else:
            key, s = canon(w)
            if s:
                out[key] = out.get(key, 0) + c * s
    return {k: v for k, v in out.items() if v}


def dag(op):
    return tuple((not d, p) for d, p in reversed(op))


def pid_of(key):
    """Real wavefunction: <w> == <w^dag>, so a key and its conjugate are one parameter."""
    return min(key, (tuple(reversed(key[1])), tuple(reversed(key[0]))))


def block_rows(r, t2=False):
    """Row operators O_a of each block, split by Sz-change. M[a,b] = <O_a^dag O_b>, except T2 blocks:
    M[a,b] = <{O_a^dag, O_b}> with O = a+_s a_q a_p (the 3-body parts cancel in the anticommutator)."""
    A, B = range(r), range(r, 2 * r)
    c, a = (lambda p: (True, p)), (lambda p: (False, p))
    rows = {
        "D1a": [(a(p),) for p in A],
        "D1b": [(a(p),) for p in B],
        "D2aa": [(a(q), a(p)) for p, q in combinations(A, 2)],
        "D2bb": [(a(q), a(p)) for p, q in combinations(B, 2)],
        "D2ab": [(a(q), a(p)) for p in A for q in B],
        "Q1a": [(c(p),) for p in A],
        "Q1b": [(c(p),) for p in B],
        "Q2aa": [(c(q), c(p)) for p, q in combinations(A, 2)],
        "Q2bb": [(c(q), c(p)) for p, q in combinations(B, 2)],
        "Q2ab": [(c(q), c(p)) for p in A for q in B],
        "G0": [(c(q), a(p)) for p in A for q in A] + [(c(q), a(p)) for p in B for q in B],
        "Gab": [(c(q), a(p)) for p in A for q in B],
        "Gba": [(c(q), a(p)) for p in B for q in A],
    }
    if t2:
        for p, q in combinations(range(2 * r), 2):
            for s_ in range(2 * r):
                dsz = (s_ < r) - (p < r) - (q < r)
                rows.setdefault(f"T2_{dsz:+d}", []).append((c(s_), a(q), a(p)))
    return rows


def block_entry(name, Ra, Rb):
    nf = normal_order(dag(Ra) + Rb)
    if name.startswith("T2"):
        for k, v in normal_order(Rb + dag(Ra)).items():
            nf[k] = nf.get(k, 0) + v
        nf = {k: v for k, v in nf.items() if v}
        assert all(len(k[0]) <= 2 for k in nf), "3-body terms must cancel in T2"
    return nf


D_BLOCKS = ("D1a", "D1b", "D2aa", "D2bb", "D2ab")


def _spin_ok(key, r):
    C, A = key
    return len(C) == len(A) and sum(p < r for p in C) == sum(p < r for p in A)


# ---- symmetry-adapted T2 (chem-7ko): SU(2) highest-weight vectors x chain-reflection parity.
# Valid for a nondegenerate singlet ground state of a reflection-symmetric chain: then S+|psi> = 0 and
# P|psi> = +-|psi>, so <[S+, w]> = 0 and <w> = <P w P^dag> for every word w (imposed on the RDM), and the
# T2 matrix is SU(2)-invariant, hence block-diagonal on (spin S, highest weight M = S) x parity.

def _word(key):
    return tuple((True, p) for p in key[0]) + tuple((False, p) for p in key[1])


def ad_splus(word, r):
    """[S+, w] as {canonical key: coef}, S+ = sum_i a+_{i alpha} a_{i beta}."""
    out = {}
    for i in {p % r for _, p in word}:
        s = ((True, i), (False, r + i))
        for sgn, w in ((1, s + tuple(word)), (-1, tuple(word) + s)):
            for k, v in normal_order(w).items():
                out[k] = out.get(k, 0) + sgn * v
    return {k: v for k, v in out.items() if v}


def reflect(key, r):
    """Site i -> r-1-i on both spins: (canonical key, sign)."""
    f = lambda p: (p // r) * r + r - 1 - p % r  # noqa: E731
    return canon(tuple((True, f(p)) for p in key[0]) + tuple((False, f(p)) for p in key[1]))


def sym_basis(words, r):
    """Orthogonal integer basis of ker(ad S+) within one T2 Sz-sector, split by reflection parity.

    Exact (Fraction) Gram-Schmidt per reflection orbit of spatial labels, so every column is a primitive
    integer vector and the adapted block Q^T M Q has exact integer constraint coefficients.
    Returns {+1: [(orbit, {word: int})], -1: [...]}."""
    from fractions import Fraction
    from math import gcd, lcm

    def spatial(w):
        return w[0][1] % r, tuple(sorted(p % r for _, p in w[1:]))

    groups = {}
    for w in words:
        groups.setdefault(spatial(w), []).append(w)

    def gs(u, basis):
        for b in basis:
            c = sum(x * y for x, y in zip(u, b)) / sum(y * y for y in b)
            u = [x - c * y for x, y in zip(u, b)]
        return u

    out, seen = {1: [], -1: []}, set()
    for g in groups:
        if g in seen:
            continue
        gp = (r - 1 - g[0], tuple(sorted(r - 1 - x for x in g[1])))
        seen |= {g, gp}
        ws = groups[g] + (groups[gp] if gp != g else [])
        idx = {w: i for i, w in enumerate(ws)}
        m = len(ws)
        brows = {}
        for i, w in enumerate(ws):
            for t, v in ad_splus(w, r).items():
                brows.setdefault(t, [Fraction(0)] * m)[i] = Fraction(v)
        rowspace = []
        for b in brows.values():
            b = gs(b, rowspace)
            if any(b):
                rowspace.append(b)
        perm = []
        for w in ws:
            key, s = reflect(canon(w)[0], r)
            perm.append((idx[_word(key)], s))
        for sigma in (1, -1):
            basis = []
            for i in range(m):
                v = gs([Fraction(int(j == i)) for j in range(m)], rowspace)
                u = list(v)
                for j, (pj, s) in enumerate(perm):
                    u[pj] += sigma * s * v[j]
                u = gs(u, basis)
                if any(u):
                    basis.append(u)
            for u in basis:
                L = lcm(*(x.denominator for x in u))
                ints = [int(x * L) for x in u]
                G = gcd(*ints)
                out[sigma].append((len(seen), {ws[j]: x // G for j, x in enumerate(ints) if x}))
    return out


def build(h1, eri, na, nb, t2=False):
    """Standard-form PQG SDP: min c.v  s.t.  A v = b,  v = [vec_F(X_j)], X_j PSD.

    t2=True adds T2 split by Sz; t2="sym" adds the SU(2) x reflection adapted T2 blocks plus the
    singlet/reflection equalities on the RDM (needs exactly reflection-symmetric integrals).
    A and c touch only upper-triangle positions (a <= b) of each block."""
    r = h1.shape[0]
    sym = t2 == "sym"
    rows = block_rows(r, bool(t2) and not sym)
    adapted = {}  # block name -> (sector, [(orbit, combo)])
    if sym:
        assert np.array_equal(h1, h1[::-1, ::-1]) and np.array_equal(eri, eri[::-1, ::-1, ::-1, ::-1])
        full = block_rows(r, True)
        for sec, w in (("h", full["T2_+0"]), ("q", full["T2_+1"])):  # M = 1/2 and M = 3/2 sectors
            for sigma, cols in sym_basis(w, r).items():
                name = f"T2{sec}{'+' if sigma > 0 else '-'}"
                adapted[name] = (sec, cols)
                rows[name] = [c for _, c in cols]
    names = list(rows)
    size = {k: len(rows[k]) for k in names}
    off, tot = {}, 0
    for k in names:
        off[k], tot = tot, tot + size[k] ** 2

    def pos(k, a, b):
        a, b = min(a, b), max(a, b)
        return off[k] + a + b * size[k]

    loc = {}
    for k in D_BLOCKS:
        R = rows[k]
        for a_ in range(len(R)):
            for b_ in range(a_, len(R)):
                ((key, coef),) = normal_order(dag(R[a_]) + R[b_]).items()
                assert coef == 1 and pid_of(key) not in loc
                loc[pid_of(key)] = pos(k, a_, b_)

    def linform(nf):
        terms, const = {}, 0
        for key, c in nf.items():
            if key == ((), ()):
                const += c
            elif _spin_ok(key, r):
                p = loc[pid_of(key)]
                terms[p] = terms.get(p, 0) + c
        return terms, const

    from array import array

    ri, ci, vals, bvec = array("q"), array("q"), array("d"), []  # n=12 sym: 5.15e6 entries (= nnz A)

    def add_row(terms, rhs):
        i = len(bvec)
        for p, c in terms.items():
            if c:
                ri.append(i), ci.append(p), vals.append(c)
        bvec.append(rhs)

    for sec in ("h", "q"):  # adapted entries X[a,b] = sum_{w,w'} Q[w,a] Q[w',b] <{w^dag, w'}>
        blocks = [k for k in adapted if adapted[k][0] == sec]
        if not blocks:
            continue
        norb = 1 + max(o for k in blocks for o, _ in adapted[k][1])
        by_orb = {k: [[] for _ in range(norb)] for k in blocks}
        for k in blocks:
            for a_, (o, combo) in enumerate(adapted[k][1]):
                by_orb[k][o].append((a_, combo))
        for o1 in range(norb):
            for o2 in range(o1, norb):
                F = {}
                for k in blocks:
                    for a_, ca in by_orb[k][o1]:
                        for b_, cb in by_orb[k][o2]:
                            if b_ < a_:
                                continue
                            row, const = {pos(k, a_, b_): 1}, 0
                            for w1, x1 in ca.items():
                                for w2, x2 in cb.items():
                                    if (w1, w2) not in F:
                                        F[w1, w2] = linform(block_entry("T2", w1, w2))
                                    terms, c0 = F[w1, w2]
                                    const += x1 * x2 * c0
                                    for p, c in terms.items():
                                        row[p] = row.get(p, 0) - x1 * x2 * c
                            add_row(row, const)

    for k in names:
        if k in D_BLOCKS or k in adapted:
            continue
        R = rows[k]
        for a_ in range(len(R)):
            for b_ in range(a_, len(R)):
                terms, const = linform(block_entry(k, R[a_], R[b_]))
                row = {p: -c for p, c in terms.items()}
                row[pos(k, a_, b_)] = row.get(pos(k, a_, b_), 0) + 1
                add_row(row, const)

    N = {0: na, 1: nb}
    for s in (0, 1):
        orb = range(s * r, (s + 1) * r)
        add_row({pos(("D1a", "D1b")[s], i, i): 1 for i in range(r)}, N[s])
        for t in (0, 1):
            for i in orb:
                for k in orb:
                    if k < i:
                        continue
                    row = {}
                    for j in range(t * r, (t + 1) * r):
                        terms, const = linform(normal_order(((True, i), (True, j), (False, j), (False, k))))
                        assert const == 0
                        for p, c in terms.items():
                            row[p] = row.get(p, 0) + c
                    (p1, c1), = linform(normal_order(((True, i), (False, k))))[0].items()
                    row[p1] = row.get(p1, 0) - (N[t] - (s == t)) * c1
                    add_row(row, 0)

    if sym:  # singlet: <[S+, w]> = 0 for every lowering w of rank <= 2; reflection: <w> = <P w P^dag>
        seen = set()
        for pid, p in loc.items():
            k2, s_ = reflect(pid, r)
            q = loc[pid_of(k2)]
            if q != p or s_ < 0:
                key = (min(p, q), max(p, q), s_)
                if key not in seen:
                    seen.add(key)
                    add_row({p: 1, q: -s_} if q != p else {p: 1}, 0)
        modes = range(2 * r)
        for nb_ in (1, 2):
            for C in combinations(modes, nb_):
                for A_ in combinations(modes, nb_):
                    if sum(x < r for x in C) - sum(x < r for x in A_) != -1:  # Sz(w) = -1
                        continue
                    terms, const = linform(ad_splus(_word((C, A_[::-1])), r))
                    assert const == 0
                    key = frozenset(terms.items())
                    if terms and key not in seen:
                        seen.add(key)
                        add_row(terms, 0)

    A = sp.csr_matrix((np.frombuffer(vals), (np.frombuffer(ri, np.int64), np.frombuffer(ci, np.int64))),
                      shape=(len(bvec), tot))
    A.sum_duplicates()

    c, cabs, ccnt = np.zeros(tot), np.zeros(tot), np.zeros(tot)

    def addc(word, val):
        key, s = canon(word)
        if s:
            p = loc[pid_of(key)]
            c[p] += s * val
            cabs[p] += abs(val)
            ccnt[p] += 1

    for s in (0, 1):
        for i in range(r):
            for k in range(r):
                addc(((True, s * r + i), (False, s * r + k)), h1[i, k])
    for s in (0, 1):
        for t in (0, 1):
            for i in range(r):
                for j in range(r):
                    for k in range(r):
                        for l_ in range(r):
                            v = eri[i, k, j, l_]
                            if v:
                                addc(((True, s * r + i), (True, t * r + j), (False, t * r + l_), (False, s * r + k)),
                                     0.5 * v)
    cerr = gamma(ccnt + 2) * cabs

    occ = set(range(na)) | set(range(r, r + nb))

    def det_value(nf):  # canonical keys: <a+C a_A> on a determinant is 1 iff A == reversed(C), all occupied
        return sum(v for (C, Ann), v in nf.items() if set(C) == set(Ann) and set(C) <= occ)

    traces = {k: sum(det_value(block_entry(k, R, R)) for R in rows[k]) for k in names if k not in adapted}
    # groups: (upper bound on sum_j tr X*_j, blocks j). An adapted block is Q^T M Q with orthogonal integer
    # columns of one Sz-sector M, so its blocks' traces sum to <= max_a |q_a|^2 * tr M (exact N-only trace).
    groups = {k: (t, [k]) for k, t in traces.items()}
    for sec, k0 in (("h", "T2_+0"), ("q", "T2_+1")):
        ks = [k for k in adapted if adapted[k][0] == sec]
        if ks:
            T = sum(det_value(block_entry("T2", w, w)) for w in block_rows(r, True)[k0])
            maxd = max(sum(x * x for x in cmb.values()) for k in ks for _, cmb in adapted[k][1])
            groups[sec] = (maxd * T, ks)
    return dict(A=A, b=np.array(bvec, float), c=c, cerr=cerr, names=names, size=size, off=off,
                rows=rows, loc=loc, traces=traces, groups=groups, tot=tot, r=r, na=na, nb=nb)


def solve(prob, solver="CLARABEL", **kw):
    import cvxpy as cp

    X = {k: cp.Variable((prob["size"][k],) * 2, PSD=True) for k in prob["names"]}
    v = cp.hstack([cp.vec(X[k], order="F") for k in prob["names"]])
    con = prob["A"] @ v == prob["b"]
    P = cp.Problem(cp.Minimize(prob["c"] @ v), [con])
    P.solve(solver=solver, **kw)
    return P.value, np.asarray(con.dual_value, float), P.status


def solve_scs(prob, **settings):
    """SCS called directly on svec variables (cvxpy's canonicalization does not fit n = 12 in 16 GB).

    Column scaling svec <-> vec does not change the row duals, so y feeds `certify` unchanged."""
    import scs

    cols, scale, sizes = [], [], []
    for k in prob["names"]:
        n, o = prob["size"][k], prob["off"][k]
        a, b = np.triu_indices(n)  # row-major upper == SCS's column-major lower triangle
        cols.append(o + a + b * n)
        scale.append(np.where(a == b, 1.0, 1 / math.sqrt(2)))
        sizes.append(n)
    cols, scale = np.concatenate(cols), np.concatenate(scale)
    m, nv = prob["A"].shape[0], len(cols)
    A = prob["A"].tocsc()[:, cols] @ sp.diags(scale)
    data = dict(A=sp.vstack([A, -sp.identity(nv)]).tocsc(), b=np.concatenate([prob["b"], np.zeros(nv)]),
                c=prob["c"][cols] * scale)
    del A
    sol = scs.SCS(data, dict(z=m, s=sizes), **settings).solve()
    return sol["info"]["pobj"], np.asarray(sol["y"][:m]), sol["info"]["status"]


def _lam_lower(Z):
    """Rigorous lower bound on lambda_min of the float matrix Z (Higham Thm 10.3)."""
    n = Z.shape[0]
    scale = max(1.0, float(np.abs(Z).max()))
    lam = float(np.linalg.eigvalsh(Z)[0])
    margin = 1e-13 * n * scale
    for _ in range(40):
        shift = lam - margin
        M = Z - shift * np.eye(n)
        try:
            R = np.linalg.cholesky(M)
        except np.linalg.LinAlgError:
            margin *= 4
            continue
        s = float(np.sum(R * R)) * (1 + gamma(n * n + 1))
        return shift - gamma(n + 1) * s - U * float(np.abs(np.diag(M)).max()) * 1.01
    raise RuntimeError("Cholesky never succeeded")


def rigorous_lower_bound(prob, y, e_const=0.0):
    """Lower bound on E0 valid for any y; returns (bound, per-block lambda bounds)."""
    A, b = prob["A"], prob["b"]
    w = A.T @ y
    wabs = abs(A).T @ np.abs(y)
    K = int(np.diff(A.tocsc().indptr).max()) + 2
    bound = float(b @ y)
    bound -= gamma(len(b) + 1) * float(np.abs(b) @ np.abs(y))
    lams = {}
    for k in prob["names"]:
        n, o = prob["size"][k], prob["off"][k]
        sl = slice(o, o + n * n)
        W = w[sl].reshape((n, n), order="F")
        Wa = wabs[sl].reshape((n, n), order="F")
        Cm = prob["c"][sl].reshape((n, n), order="F")
        Ce = prob["cerr"][sl].reshape((n, n), order="F")
        Z = (Cm + Cm.T) / 2 - (W + W.T) / 2
        E = gamma(K) * (Wa + Wa.T) / 2 + (Ce + Ce.T) / 2 + U * np.abs(Z)
        lam = _lam_lower(Z) - math.sqrt(float(np.sum(E * E))) * (1 + 1e-6)
        lams[k] = lam
    for T, ks in prob["groups"].values():
        bound += T * min(0.0, *(lams[k] for k in ks)) * (1 + 1e-12)
    return bound + e_const - 1e-14 * (abs(bound) + abs(e_const)), lams


def certify(prob, y, e_const):
    """Best rigorous bound over the two dual-sign conventions (each is valid on its own)."""
    return max((rigorous_lower_bound(prob, s * y, e_const) for s in (1.0, -1.0)), key=lambda t: t[0])


def sym_integrals(n):
    """Löwdin site integrals averaged over the chain reflection, so P H P^dag == H exactly in float."""
    from benchmark_hchain_tdl import integrals

    h1, eri, ne, e_nuc, _ = integrals(n, localize=True)
    return (h1 + h1[::-1, ::-1]) / 2, (eri + eri[::-1, ::-1, ::-1, ::-1]) / 2, ne, e_nuc


def hchain_problem(n, t2=False):
    from benchmark_hchain_tdl import integrals

    if t2 == "sym":
        h1, eri, (na, nb), e_nuc = sym_integrals(n)
    else:
        h1, eri, (na, nb), e_nuc, _ = integrals(n, localize=True)
    return build(h1, eri, na, nb, t2), e_nuc


FCI_MAX_N = 12


def singlet_premise(n):
    """The T2sym premise by converged FCI on the symmetrized integrals: (E0, gap, <S^2> of root 0) over
    the lowest two Sz = 0 roots. Raises unless Davidson converged and root 0 is a nondegenerate singlet.
    pyscf's direct_spin1.kernel defaults (50 iterations, unchecked) reported <S^2> = 0.398 at n = 12 for
    a vector 27 mHa above E0 (chem-7ko, SPEC_hchain_t2sym §6)."""
    from pyscf import fci

    h1, eri, ne, e_nuc = sym_integrals(n)
    s = fci.direct_spin1.FCI()
    s.max_cycle, s.conv_tol, s.nroots = 1000, 1e-10, 2
    e, ci = s.kernel(h1, eri, len(h1), ne)
    gap, s2 = e[1] - e[0], fci.spin_op.spin_square0(ci[0], len(h1), ne)[0]
    if not (np.all(s.converged) and gap > 1e-3 and s2 < 1e-6):
        raise RuntimeError(f"n={n}: singlet premise not verified (converged {s.converged}, gap {gap}, <S^2> {s2})")
    return e[0] + e_nuc, gap, s2


def upper_bound(n, sym=False):
    """Variational upper bound: exact FCI up to n=12, else the DMRG energy from --dmrg-upper."""
    if n <= FCI_MAX_N:
        from benchmark_hchain_tdl import integrals
        from hybrid_quantum_solver.dmrg_reference import fci_energy

        h1, eri, ne, e_nuc = sym_integrals(n) if sym else integrals(n, localize=True)[:4]
        return fci_energy(h1, eri, ne, e_nuc), "FCI"
    path = "data/hchain_vsdp_upper.json"
    ub = json.load(open(path)) if os.path.exists(path) else {}
    return (ub[str(n)]["energy"], ub[str(n)]["method"]) if str(n) in ub else (None, None)


def dmrg_upper(n, bond_dims, threads):
    """Variational DMRG upper bound (final two-site sweep energy, NOT the extrapolation). block2 must
    not share a process with pyscf, so the integrals come from a child process."""
    import subprocess
    import sys

    path = f".dmrg_tmp/hchain_{n}_integrals.npz"
    os.makedirs(".dmrg_tmp", exist_ok=True)
    subprocess.run([sys.executable, __file__, "--dump", str(n), path], check=True)
    d = np.load(path)
    from hybrid_quantum_solver.dmrg_reference import dmrg_energy

    na = int(d["na"])
    e = dmrg_energy(d["h1"], d["eri"], (na, na), float(d["e_nuc"]), bond_dims=bond_dims,
                    noises=(1e-4,) * (len(bond_dims) - 2) + (1e-5, 0.0), n_sweeps=4 * len(bond_dims),
                    n_threads=threads)
    os.makedirs("data", exist_ok=True)
    out = "data/hchain_vsdp_upper.json"
    ub = json.load(open(out)) if os.path.exists(out) else {}
    ub[str(n)] = dict(energy=e, method=f"DMRG SU2 variational, D={'/'.join(map(str, bond_dims))}")
    with open(out, "w") as f:
        json.dump(ub, f, indent=2)
    print(n, ub[str(n)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", default="2,4,6,8,10,12")
    ap.add_argument("--solver", default="SCS", help="SCS; CLARABEL is exact-er but ~(block size)^6")
    ap.add_argument("--eps", type=float, default=1e-8)
    ap.add_argument("--max-iters", type=int, default=200000)
    ap.add_argument("--output", help="default data/hchain_vsdp[_t2|_t2sym].json (rows are keyed by n only)")
    ap.add_argument("--dump", nargs=2, metavar=("N", "PATH"), help=argparse.SUPPRESS)
    ap.add_argument("--dmrg-upper", type=int, metavar="N")
    ap.add_argument("--bond-dims", default="250,500,1000,1500")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--t2", action="store_true", help="add the T2 condition (block ~1.5 r^3: n <= 8 here)")
    ap.add_argument("--sym", action="store_true", help="SU(2) x reflection adapted T2, direct SCS (chem-7ko)")
    ap.add_argument("--table", action="store_true", help="re-pair saved bounds with current upper bounds")
    args = ap.parse_args()
    if args.sym and (args.solver != "SCS" or max(map(int, args.ns.split(","))) > FCI_MAX_N):
        ap.error(f"--sym always runs SCS and needs the FCI-checked singlet premise, so n <= {FCI_MAX_N}")

    if args.dump:
        from benchmark_hchain_tdl import integrals

        h1, eri, (na, _), e_nuc, _ = integrals(int(args.dump[0]), localize=True)
        np.savez(args.dump[1], h1=h1, eri=eri, na=na, e_nuc=e_nuc)
        return
    if args.dmrg_upper:
        dmrg_upper(args.dmrg_upper, tuple(map(int, args.bond_dims.split(","))), args.threads)
        return

    args.output = args.output or f"data/hchain_vsdp{'_t2sym' if args.sym else '_t2' if args.t2 else ''}.json"
    results = {}
    if os.path.exists(args.output):
        with open(args.output) as f:
            results = json.load(f)
    os.makedirs("data", exist_ok=True)
    for n in map(int, args.ns.split(",")):
        premise = singlet_premise(n) if args.sym else None  # seconds, before hours of SDP
        if args.table:
            row = results[str(n)]
        else:
            prob, e_nuc = hchain_problem(n, "sym" if args.sym else args.t2)
            kw = dict(eps_abs=args.eps, eps_rel=args.eps, max_iters=args.max_iters) if args.solver == "SCS" else {}
            if args.sym:
                val, y, status = solve_scs(prob, **kw)
            else:
                val, y, status = solve(prob, args.solver, **kw)
            lb, lams = certify(prob, y, e_nuc)
            cond = "PQG+T2sym" if args.sym else "PQG+T2" if args.t2 else "PQG"
            row = dict(n=n, conditions=cond, solver=args.solver, status=status, sdp_value=val + e_nuc, certified_lb=lb,
                       worst_block_lambda=min(lams.values()),
                       t2_blocks={k: v for k, v in prob["size"].items() if k.startswith("T2")})
        up, how = upper_bound(n, args.sym)
        if premise:
            assert abs(premise[0] - up) < 1e-8, (premise[0], up)  # the premise state is the bounding E0
            row.update(fci_gap=premise[1], fci_s2=premise[2])
        row.update(upper=up, upper_method=how,
                   width_mHa_per_atom=None if up is None else (up - row["certified_lb"]) / n * 1e3)
        results[str(n)] = row
        print(json.dumps(row), flush=True)
        with open(args.output, "w") as f:
            json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
