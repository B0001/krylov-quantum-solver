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


def build(h1, eri, na, nb, t2=False):
    """Standard-form PQG SDP: min c.v  s.t.  A v = b,  v = [vec_F(X_j)], X_j PSD.

    A and c touch only upper-triangle positions (a <= b) of each block."""
    r = h1.shape[0]
    rows = block_rows(r, t2)
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

    ri, ci, vals, bvec = [], [], [], []

    def add_row(terms, rhs):
        i = len(bvec)
        for p, c in terms.items():
            if c:
                ri.append(i), ci.append(p), vals.append(c)
        bvec.append(rhs)

    for k in names:
        if k in D_BLOCKS:
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

    A = sp.csr_matrix((np.array(vals, float), (ri, ci)), shape=(len(bvec), tot))
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

    traces = {k: sum(det_value(block_entry(k, R, R)) for R in rows[k]) for k in names}
    return dict(A=A, b=np.array(bvec, float), c=c, cerr=cerr, names=names, size=size, off=off,
                rows=rows, loc=loc, traces=traces, tot=tot, r=r, na=na, nb=nb)


def solve(prob, solver="CLARABEL", **kw):
    import cvxpy as cp

    X = {k: cp.Variable((prob["size"][k],) * 2, PSD=True) for k in prob["names"]}
    v = cp.hstack([cp.vec(X[k], order="F") for k in prob["names"]])
    con = prob["A"] @ v == prob["b"]
    P = cp.Problem(cp.Minimize(prob["c"] @ v), [con])
    P.solve(solver=solver, **kw)
    return P.value, np.asarray(con.dual_value, float), P.status


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
        bound += prob["traces"][k] * min(0.0, lam) * (1 + 1e-12)
    return bound + e_const - 1e-14 * (abs(bound) + abs(e_const)), lams


def certify(prob, y, e_const):
    """Best rigorous bound over the two dual-sign conventions (each is valid on its own)."""
    return max((rigorous_lower_bound(prob, s * y, e_const) for s in (1.0, -1.0)), key=lambda t: t[0])


def hchain_problem(n, t2=False):
    from benchmark_hchain_tdl import integrals

    h1, eri, (na, nb), e_nuc, _ = integrals(n, localize=True)
    return build(h1, eri, na, nb, t2), e_nuc


FCI_MAX_N = 12


def upper_bound(n):
    """Variational upper bound: exact FCI up to n=12, else the DMRG energy from --dmrg-upper."""
    if n <= FCI_MAX_N:
        from benchmark_hchain_tdl import integrals
        from hybrid_quantum_solver.dmrg_reference import fci_energy

        h1, eri, ne, e_nuc, _ = integrals(n, localize=True)
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
    ap.add_argument("--output", default="data/hchain_vsdp.json")
    ap.add_argument("--dump", nargs=2, metavar=("N", "PATH"), help=argparse.SUPPRESS)
    ap.add_argument("--dmrg-upper", type=int, metavar="N")
    ap.add_argument("--bond-dims", default="250,500,1000,1500")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--t2", action="store_true", help="add the T2 condition (block ~1.5 r^3: n <= 8 here)")
    ap.add_argument("--table", action="store_true", help="re-pair saved bounds with current upper bounds")
    args = ap.parse_args()

    if args.dump:
        from benchmark_hchain_tdl import integrals

        h1, eri, (na, _), e_nuc, _ = integrals(int(args.dump[0]), localize=True)
        np.savez(args.dump[1], h1=h1, eri=eri, na=na, e_nuc=e_nuc)
        return
    if args.dmrg_upper:
        dmrg_upper(args.dmrg_upper, tuple(map(int, args.bond_dims.split(","))), args.threads)
        return

    results = {}
    if os.path.exists(args.output):
        with open(args.output) as f:
            results = json.load(f)
    os.makedirs("data", exist_ok=True)
    for n in map(int, args.ns.split(",")):
        if args.table:
            row = results[str(n)]
        else:
            prob, e_nuc = hchain_problem(n, args.t2)
            kw = dict(eps_abs=args.eps, eps_rel=args.eps, max_iters=args.max_iters) if args.solver == "SCS" else {}
            val, y, status = solve(prob, args.solver, **kw)
            lb, lams = certify(prob, y, e_nuc)
            row = dict(n=n, conditions="PQG+T2" if args.t2 else "PQG", solver=args.solver, status=status, sdp_value=val + e_nuc, certified_lb=lb,
                       worst_block_lambda=min(lams.values()))
        up, how = upper_bound(n)
        row.update(upper=up, upper_method=how,
                   width_mHa_per_atom=None if up is None else (up - row["certified_lb"]) / n * 1e3)
        results[str(n)] = row
        print(json.dumps(row), flush=True)
        with open(args.output, "w") as f:
            json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
