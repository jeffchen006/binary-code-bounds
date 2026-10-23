"""Schrijver SDP + 2-split Terwilliger SDP (Tseng-Lai-Yu) + valid linear constraints for A(n,d).

Variables: x[key], key = sorted triple of pairs ((a1,a2),(b1,b2),(c1,c2)), the split pairwise
distances of a codeword triple (X,Y,Z) on parts T1 (size n1) and T2 (size n2).
Value = average over X in C of fraction of (Y,Z) in C^2 among all pairs with that configuration.
With n1 = 0 this is exactly Schrijver's SDP.
"""
from fractions import Fraction
from functools import lru_cache
import math, sys, itertools
import mpmath
from sdpcore import Problem, binom
import bounds


@lru_cache(None)
def beta(n, i, j, k, t):
    s = 0
    for u in range(max(t, k), n + 1):
        b = binom(u, t) * binom(n - 2 * k, u - k) * binom(n - k - u, i - u) * binom(n - k - u, j - u)
        if b:
            s += (1 if (u - t) % 2 == 0 else -1) * b
    return s


def multi(n, p, q, r):
    if min(p, q, r) < 0 or p + q + r > n:
        return 0
    return math.factorial(n) // (math.factorial(p) * math.factorial(q) * math.factorial(r) * math.factorial(n - p - q - r))


def build(n, d, n1=0, even=None, opts=None, name=None):
    opts = dict(schrijver=True, split=True, lin_basic=True, best18=True, cw22=True,
                mel19=True, mel21=True, best20=True, kt23=True) | (opts or {})
    if even is None:
        even = (d % 2 == 0)
    n2 = n - n1
    parts = (n1, n2) if n1 > 0 else (n,)
    P = Problem(name or 'split_%d_%d_%d' % (n, d, n1))

    def okdist(a):
        if a == 0:
            return True
        if a < d:
            return False
        if even and a % 2:
            return False
        return True

    def key(i, j, t, i2, j2, t2):
        """split triple -> None (zero) / 'one' / var key"""
        if n1 > 0:
            prs = ((i, i2), (j, j2), (i + j - 2 * t, i2 + j2 - 2 * t2))
        else:
            prs = ((i,), (j,), (i + j - 2 * t,))
        for p in prs:
            if not okdist(sum(p)):
                return None
        k = tuple(sorted(prs))
        if all(sum(p) == 0 for p in k):
            return 'one'
        return P.var(k)

    # all configurations
    def confs1(m):
        for i in range(m + 1):
            for j in range(m + 1):
                for t in range(max(0, i + j - m), min(i, j) + 1):
                    yield i, j, t
    if n1 > 0:
        C1 = list(confs1(n1)); C2 = list(confs1(n2))
    else:
        C1 = list(confs1(n)); C2 = [(0, 0, 0)]

    def split_key(c1, c2):
        return key(c1[0], c1[1], c1[2], c2[0], c2[1], c2[2])

    X0 = {}  # (i,i2) -> key of x_{(i,i2),0}
    for i in range(parts[0] + 1):
        for i2 in range((n2 if n1 > 0 else 0) + 1):
            X0[(i, i2)] = key(i, 0, 0, i2, 0, 0)

    # objective
    P.c0 = Fraction(1)
    for (i, i2), v in X0.items():
        if v is not None and v != 'one':
            co = binom(parts[0], i) * (binom(n2, i2) if n1 > 0 else 1)
            P.c[v] = P.c.get(v, 0) + co

    # distance distribution A_a as linear form
    def Aform(a):
        f = {}
        for (i, i2), v in X0.items():
            if i + i2 == a and v is not None:
                co = binom(parts[0], i) * (binom(n2, i2) if n1 > 0 else 1)
                f[v] = f.get(v, 0) + co
        return f

    lin = {}
    def addlin(g0, co):
        co = {a: Fraction(b) for a, b in co.items() if b != 0}
        if 'one' in co:
            g0 = g0 + co.pop('one')
        if not co:
            assert g0 >= 0, 'infeasible constant constraint'
            return
        k = (Fraction(g0), tuple(sorted(co.items())))
        lin[k] = 1

    # ---- basic linear constraints on split variables ----
    if opts['lin_basic']:
        for c1 in C1:
            for c2 in C2:
                v = split_key(c1, c2)
                if v == 'one':
                    continue
                pa = (c1[0], c2[0]) if n1 > 0 else (c1[0], 0)
                pb = (c1[1], c2[1]) if n1 > 0 else (c1[1], 0)
                xa, xb = X0.get(pa), X0.get(pb)
                if v is not None:
                    for w in (xa, xb):
                        if w is not None and w != v:
                            addlin(0, {w: 1, v: -1})
                if sum(pa) > 0 and sum(pb) > 0 and xa is not None and xb is not None:
                    co = {}
                    for w, s in ((xa, -1), (xb, -1), (v, 1)):
                        if w is not None:
                            co[w] = co.get(w, 0) + s
                    addlin(1, co)

    # ---- Schrijver variables as linear forms in split variables ----
    def schr_form(a, b, c):
        """x_{a,b}^c (full length n) as dict"""
        if n1 == 0:
            v = key(a, b, c, 0, 0, 0)
            return {} if v is None else {v: 1}
        f = {}
        den = multi(n, a - c, b - c, c)
        for i in range(max(0, a - n2), min(a, n1) + 1):
            i2 = a - i
            for j in range(max(0, b - n2), min(b, n1) + 1):
                j2 = b - j
                for t in range(max(0, i + j - n1), min(i, j) + 1):
                    t2 = c - t
                    if t2 < max(0, i2 + j2 - n2) or t2 > min(i2, j2):
                        continue
                    v = key(i, j, t, i2, j2, t2)
                    if v is None:
                        continue
                    num = multi(n1, i - t, j - t, t) * multi(n2, i2 - t2, j2 - t2, t2)
                    f[v] = f.get(v, 0) + Fraction(num, den)
        return f

    def schr0(a):
        return schr_form(a, 0, 0)

    if n1 > 0 and opts['lin_basic']:
        for a, b, c in confs1(n):
            fv = schr_form(a, b, c)
            fa, fb = schr0(a), schr0(b)
            if fv:
                for fw in (fa, fb):
                    co = dict(fw)
                    for v, s in fv.items():
                        co[v] = co.get(v, 0) - s
                    addlin(0, co)
            if a > 0 and b > 0:
                co = {}
                for f, s in ((fa, -1), (fb, -1), (fv, 1)):
                    for v, w in f.items():
                        co[v] = co.get(v, 0) + s * w
                addlin(1, co)

    # ---- PSD blocks ----
    def add_psd(rowsets, entry, scale):
        """rowsets: list of row labels; entry(r,s,typ) -> dict var->coef (None key = const)"""
        pass

    def finish_block(idx, ent, scl, meta=None):
        used = set()
        for kk, mat in ent.items():
            for (i, j), val in mat.items():
                if val != 0:
                    used.add(i); used.add(j)
        keep = [r for r in idx if r in used]
        if not keep:
            return
        pos = {r: p for p, r in enumerate(keep)}
        ent2 = {}
        for kk, mat in ent.items():
            m2 = {(pos[i], pos[j]) if pos[i] <= pos[j] else (pos[j], pos[i]): Fraction(val)
                  for (i, j), val in mat.items() if val != 0}
            if m2:
                ent2[kk] = m2
        P.add_block(len(keep), ent2, [scl[r] for r in keep])
        P.blocks[-1]['meta'] = dict(meta or {}, rows=keep)

    def accum(ent, r, s, form, mult):
        for v, cf in form.items():
            kk = None if v == 'one' else v
            ent.setdefault(kk, {})
            ent[kk][(r, s)] = ent[kk].get((r, s), 0) + mult * cf

    if opts['schrijver']:
        for k in range(0, n // 2 + 1):
            idx = list(range(k, n - k + 1))
            for typ in (0, 1):
                ent = {}
                for ii, i in enumerate(idx):
                    for j in idx[ii:]:
                        for t in range(max(0, i + j - n), min(i, j) + 1):
                            be = beta(n, i, j, k, t)
                            if be == 0:
                                continue
                            if n1 == 0:
                                v = key(i, j, t, 0, 0, 0)
                                fv = {} if v is None else {v: 1}
                            else:
                                fv = schr_form(i, j, t)
                            if typ == 0:
                                accum(ent, i, j, fv, be)
                            else:
                                f0 = schr0(i + j - 2 * t) if n1 > 0 else ({} if key(i + j - 2 * t, 0, 0, 0, 0, 0) is None else {key(i + j - 2 * t, 0, 0, 0, 0, 0): 1})
                                accum(ent, i, j, f0, be)
                                accum(ent, i, j, fv, -be)
                scl = {i: 1 / mpmath.sqrt(binom(n - 2 * k, i - k)) for i in idx}
                finish_block(idx, ent, scl, dict(kind='schr', k=k, typ=typ))

    if n1 > 0 and opts['split']:
        for k in range(0, n1 // 2 + 1):
            for k2 in range(0, n2 // 2 + 1):
                idx = [(i, i2) for i in range(k, n1 - k + 1) for i2 in range(k2, n2 - k2 + 1)]
                for typ in (0, 1):
                    ent = {}
                    for ii, r in enumerate(idx):
                        for s in idx[ii:]:
                            i, i2 = r; j, j2 = s
                            for t in range(max(0, i + j - n1), min(i, j) + 1):
                                b1 = beta(n1, i, j, k, t)
                                if b1 == 0:
                                    continue
                                for t2 in range(max(0, i2 + j2 - n2), min(i2, j2) + 1):
                                    b2 = beta(n2, i2, j2, k2, t2)
                                    if b2 == 0:
                                        continue
                                    v = key(i, j, t, i2, j2, t2)
                                    be = b1 * b2
                                    if typ == 0:
                                        if v is not None:
                                            accum(ent, r, s, {v: 1}, be)
                                    else:
                                        w = X0.get((i + j - 2 * t, i2 + j2 - 2 * t2))
                                        if w is not None:
                                            accum(ent, r, s, {w: 1}, be)
                                        if v is not None:
                                            accum(ent, r, s, {v: 1}, -be)
                    scl = {(i, i2): 1 / mpmath.sqrt(binom(n1 - 2 * k, i - k) * binom(n2 - 2 * k2, i2 - k2)) for (i, i2) in idx}
                    finish_block(idx, ent, scl, dict(kind='split', k=k, k2=k2, typ=typ))

    # ---- linear constraints on distance distribution ----
    A = {a: Aform(a) for a in range(n + 1)}
    A[0] = {'one': 1}
    def comb_forms(terms):
        f = {}
        for coef, a in terms:
            if a < 0 or a > n:
                continue
            for v, c in A[a].items():
                f[v] = f.get(v, 0) + coef * c
        return f
    h = d // 2 if d % 2 == 0 else None
    if opts['best18']:
        for kk in range(1, n):
            try:
                ub = bounds.A_upper(n - kk, d)
            except Exception:
                continue
            f = comb_forms([(-binom(n - a, kk), a) for a in range(n + 1)])
            addlin(binom(n, kk) * ub, f)
    if opts['cw22']:
        for a in range(d, n + 1):
            f = comb_forms([(-1, a)])
            addlin(bounds.cw_upper(n, d, a), f)
    if h is not None and opts['mel19']:
        q = (2 * n) // d
        f = comb_forms([(-1, n - h)] + [(-q, n - i) for i in range(h)])
        addlin(q, f)
    if h is not None and opts['mel21']:
        for i in range(1, h):
            S = bounds.cw_upper(n, d, h + i)
            Tt = bounds.cw_upper(n - h + i, d, h + i)
            f = comb_forms([(-1, n - h - i), (-(S - Tt), n - h + i)] + [(-S, n - h + j) for j in range(i + 1, h + 1)])
            addlin(S, f)
    if h is not None and opts['best20']:
        Pp = bounds.cw_upper(n - 1, d, h + 1)
        Q = bounds.cw_upper(n - h, d, h + 1)
        R = bounds.cw_upper(n - h + 2, d, h + 2)
        terms = [(-(h + 2), n - h - 2), (-h * (Pp - Q), n - h), (-(n * Pp - (h + 2) * R), n - h + 2)]
        terms += [(-n * Pp, a) for a in range(n - h + 3, n + 1)]
        addlin(n * Pp, comb_forms(terms))
    if opts['kt23']:
        # x_{i,j}^t <= T(t,i,j-t,n-i,d)/(C(i,t)C(n-i,j-t)) x_{i,0}^0  (full-length Schrijver variables)
        for i in range(d, n + 1):
            for j in range(d, n + 1):
                for t in range(max(0, i + j - n), min(i, j) + 1):
                    c = i + j - 2 * t
                    if not (okdist(i) and okdist(j) and okdist(c)) or c == 0:
                        continue
                    Tv = bounds.T_upper(t, i, j - t, n - i, d)
                    tot = binom(i, t) * binom(n - i, j - t)
                    if Tv >= tot:
                        continue
                    fv = schr_form(i, j, t)
                    f0 = schr0(i)
                    co = {v: Fraction(Tv, tot) * w for v, w in f0.items()}
                    for v, w in fv.items():
                        co[v] = co.get(v, 0) - w
                    addlin(0, co)
    P.lin = []
    for (g0, co) in lin.keys():
        mx = max([abs(c) for _, c in co] + [abs(g0)])
        P.lin.append((g0 / mx, {v: c / mx for v, c in co}))
    return P


if __name__ == '__main__':
    import time
    n, d, n1 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    S = '/tmp/work'
    t0 = time.time()
    P = build(n, d, n1)
    print(len(P.vars), 'vars', len(P.blocks), 'blocks', max(b['size'] for b in P.blocks), 'maxblk', len(P.lin), 'lin', '%.1fs' % (time.time() - t0), flush=True)
    info = P.solve(S)
    print(n, d, n1, info['phase.value'], info['bound'], '%.1fs' % (time.time() - t0))
