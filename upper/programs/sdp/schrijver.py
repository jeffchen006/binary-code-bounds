"""Schrijver (2005) 3-point SDP bound for A(n,d), block-diagonalized Terwilliger algebra.
Variables x_{i,j}^t, identified by sorted triple of pairwise distances (a<=b<=c).
"""
from fractions import Fraction
from functools import lru_cache
import math, sys, mpmath
from sdpcore import Problem, binom


def allowed(a, d, even):
    if a == 0:
        return True
    if a < d:
        return False
    if even and a % 2:
        return False
    return True


def tri(i, j, t):
    return tuple(sorted((i, j, i + j - 2 * t)))


@lru_cache(None)
def beta(n, i, j, k, t):
    s = 0
    for u in range(0, n + 1):
        b = binom(u, t) * binom(n - 2 * k, u - k) * binom(n - k - u, i - u) * binom(n - k - u, j - u)
        if b:
            s += (1 if (u - t) % 2 == 0 else -1) * b
    return s


def build(n, d, even=None, extra=None, name=None):
    if even is None:
        even = (d % 2 == 0)
    P = Problem(name or 'schr_%d_%d' % (n, d))
    def X(i, j, t):
        """returns None (=0), 'one' (=1), or var key"""
        a, b, c = tri(i, j, t)
        if not (allowed(a, d, even) and allowed(b, d, even) and allowed(c, d, even)):
            return None
        if (a, b, c) == (0, 0, 0):
            return 'one'
        return P.var((a, b, c))
    # objective: sum_i C(n,i) x_{i,0}^0
    P.c0 = Fraction(1)
    for i in range(d, n + 1):
        v = X(i, 0, 0)
        if v is not None:
            P.c[v] = P.c.get(v, 0) + binom(n, i)
    # linear constraints
    for i in range(n + 1):
        for j in range(i, n + 1):
            for t in range(max(0, i + j - n), min(i, j) + 1):
                v = X(i, j, t)
                if v == 'one':
                    continue
                xi, xj = X(i, 0, 0), X(j, 0, 0)
                # x_ij^t <= x_i0^0 ; x_ij^t <= x_j0^0
                if v is not None:
                    for w in (xi, xj):
                        if w is not None and w != 'one' and w != v:
                            P.add_lin(0, {w: 1, v: -1})
                # x_i + x_j <= 1 + x_ij^t
                if i > 0 and j > 0 and xi is not None and xj is not None and xi != 'one':
                    co = {}
                    co[xi] = co.get(xi, 0) - 1
                    co[xj] = co.get(xj, 0) - 1
                    if v is not None:
                        co[v] = co.get(v, 0) + 1
                    co = {a: b for a, b in co.items() if b}
                    P.add_lin(1, co)
    # dedupe linear constraints
    seen = set(); lin = []
    for g0, co in P.lin:
        key = (g0, tuple(sorted(co.items())))
        if key not in seen:
            seen.add(key); lin.append((g0, co))
    P.lin = lin
    # PSD blocks
    for k in range(0, n // 2 + 1):
        idx = list(range(k, n - k + 1))
        for typ in (0, 1):
            rows = []
            ent = {}
            for i in idx:
                for j in idx:
                    if j < i:
                        continue
                    for t in range(max(0, i + j - n), min(i, j) + 1):
                        be = beta(n, i, j, k, t)
                        if be == 0:
                            continue
                        v = X(i, j, t)
                        if typ == 0:
                            terms = [(v, 1)]
                        else:
                            terms = [(X(i + j - 2 * t, 0, 0), 1), (v, -1)]
                        for w, s in terms:
                            if w is None:
                                continue
                            key = None if w == 'one' else w
                            ent.setdefault(key, {})
                            ent[key][(i, j)] = ent[key].get((i, j), 0) + s * be
            # remove zero rows/cols
            used = set()
            for key, mat in ent.items():
                for (i, j), val in mat.items():
                    if val != 0:
                        used.add(i); used.add(j)
            keep = [i for i in idx if i in used]
            if not keep:
                continue
            pos = {i: p for p, i in enumerate(keep)}
            ent2 = {}
            for key, mat in ent.items():
                m2 = {(pos[i], pos[j]): Fraction(val) for (i, j), val in mat.items() if val != 0}
                if m2:
                    ent2[key] = m2
            scale = [1 / mpmath.sqrt(binom(n - 2 * k, i - k)) for i in keep]
            P.add_block(len(keep), ent2, scale)
    if extra:
        extra(P, X, n, d)
    return P


if __name__ == '__main__':
    n, d = int(sys.argv[1]), int(sys.argv[2])
    P = build(n, d)
    print(len(P.vars), 'vars', len(P.blocks), 'blocks', len(P.lin), 'lin')
    info = P.solve('/tmp/work')
    print(info['phase.value'], info['objValPrimal'], info['objValDual'])
