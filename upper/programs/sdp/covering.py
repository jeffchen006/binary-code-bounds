"""Schrijver SDP for A(N, 2e+1) (odd distance, no parity restriction) plus
 (a) linear constraints of the extended even code (length N+1, distance 2e+2) mapped to the odd code,
 (b) Johnson-type covering constraints with parameter T (valid for every code with M >= T):
     V(N,e) + |U|/M <= 2^N/T, where U = words at distance exactly e+1 from the code, and
     |U| >= sum_v [alpha k_v + beta k_v(k_v-1)] for (alpha,beta) = (2/(j+1), -1/(j(j+1))), j>=1, or (1/K, 0),
     sum_{v in U} k_v      = M (C(N,e+1) - C(2e+1,e) A_{2e+1})
     sum_{v in U} k_v(k_v-1) = M (C(2e+2,e+1) A_{2e+2} - T3)
     T3 = per-codeword # ordered pairs (c',c'') with d(c,c')=2e+2, d(c,c'')=2e+1, d(c',c'')=2e+1
 If the SDP maximum is < T then A(N,2e+1) = A(N+1,2e+2) <= T-1.
"""
from fractions import Fraction
import math, sys, time
import split, bounds
from sdpcore import binom


def xkey(i, j, t):
    return tuple(sorted(((i,), (j,), (i + j - 2 * t,))))


def build(N, dodd, T=None, opts=None, ext=True, cover=True, js=None, name=None):
    assert dodd % 2 == 1
    e = (dodd - 1) // 2
    n = N + 1
    deven = dodd + 1
    P = split.build(N, dodd, 0, even=False, opts=dict(opts or {}), name=name or 'cov_%d_%d_%s' % (N, dodd, T))
    newlin = []

    def X(i, j, t):
        k = xkey(i, j, t)
        if all(a[0] == 0 for a in k):
            return 'one'
        return k if k in P.vidx else None

    def A(i):  # odd-code distance distribution form
        if i == 0:
            return {'one': Fraction(1)}
        v = X(i, 0, 0)
        return {} if v is None else {v: Fraction(binom(N, i))}

    def Aext(a):  # extended code (length n, even) distribution: A'_a = A_{a-1} + A_a (a even), 0 for a odd
        f = {}
        if a % 2:
            return f
        for i in (a - 1, a):
            if 0 <= i <= N:
                for v, c in A(i).items():
                    f[v] = f.get(v, 0) + c
        return f

    def add(g0, co):
        co = {v: Fraction(c) for v, c in co.items() if c != 0}
        g0 = Fraction(g0)
        if 'one' in co:
            g0 += co.pop('one')
        if not co:
            assert g0 >= 0, 'constant constraint violated'
            return
        mx = max([abs(c) for c in co.values()] + [abs(g0)])
        newlin.append((g0 / mx, {v: c / mx for v, c in co.items()}))

    def comb(terms):
        f = {}
        for coef, form in terms:
            for v, c in form.items():
                f[v] = f.get(v, 0) + coef * c
        return f

    if ext:
        h = deven // 2
        # Best Lemma 7 for the extended code
        for k in range(1, n):
            try:
                ub = bounds.A_upper(n - k, deven)
            except Exception:
                continue
            add(binom(n, k) * ub, comb([(-binom(n - a, k), Aext(a)) for a in range(0, n + 1, 2)]))
        for a in range(deven, n + 1, 2):
            add(bounds.cw_upper(n, deven, a), comb([(-1, Aext(a))]))
        q = (2 * n) // deven
        add(q, comb([(-1, Aext(n - h))] + [(-q, Aext(n - i)) for i in range(h)]))
        for i in range(1, h):
            S = bounds.cw_upper(n, deven, h + i)
            Tt = bounds.cw_upper(n - h + i, deven, h + i)
            add(S, comb([(-1, Aext(n - h - i)), (-(S - Tt), Aext(n - h + i))] + [(-S, Aext(n - h + j)) for j in range(i + 1, h + 1)]))
        Pp = bounds.cw_upper(n - 1, deven, h + 1)
        Q = bounds.cw_upper(n - h, deven, h + 1)
        R = bounds.cw_upper(n - h + 2, deven, h + 2)
        terms = [(-(h + 2), Aext(n - h - 2)), (-h * (Pp - Q), Aext(n - h)), (-(n * Pp - (h + 2) * R), Aext(n - h + 2))]
        terms += [(-n * Pp, Aext(a)) for a in range(n - h + 3, n + 1)]
        add(n * Pp, comb(terms))

    if cover and T is not None:
        V = sum(binom(N, i) for i in range(e + 1))
        K = N // (e + 1)
        S1 = comb([(Fraction(binom(N, e + 1)), {'one': 1}), (-binom(2 * e + 1, e), A(2 * e + 1))])
        i, j, t = 2 * e + 2, 2 * e + 1, e + 1
        gam = split.multi(N, i - t, j - t, t)
        v3 = X(i, j, t)
        T3 = {} if v3 is None else {v3: Fraction(gam)}
        S2 = comb([(binom(2 * e + 2, e + 1), A(2 * e + 2)), (-1, T3)])
        rhs = Fraction(2 ** N, T) - V
        # 1 >= k/K
        add(rhs, comb([(-Fraction(1, K), S1)]))
        for jj in (js or range(1, K)):
            al = Fraction(2, jj + 1); be = Fraction(-1, jj * (jj + 1))
            add(rhs, comb([(-al, S1), (-be, S2)]))
    P.lin = P.lin + newlin
    return P


if __name__ == '__main__':
    N, dodd = int(sys.argv[1]), int(sys.argv[2])
    T = int(sys.argv[3]) if len(sys.argv) > 3 else None
    t0 = time.time()
    P = build(N, dodd, T)
    info = P.solve('agents/sdp/runs/cov', threads=4)
    print(N, dodd, T, info['phase.value'], info['bound'], '%.0fs' % (time.time() - t0))
