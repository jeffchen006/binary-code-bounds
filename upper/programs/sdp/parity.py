"""Exact-size exclusion tests (d even, even-weight code of length n, size exactly M).
S_u = sum_c (-1)^{u.c};  sum_{|u|=k} S_u^2 = M^2 B_k = M * sum_i A_i K_k(i).
 * M odd:        S_u odd                     => sum_i A_i K_k(i) >= C(n,k)/M
 * M = 2 mod 4:  s = sum_c c (mod 2), |s| = w even (weights even); S_u = 2 mod 4 when <u,s>=0
                 => sum_i A_i K_k(i) >= 2 (C(n,k) + K_k(w)) / M      (case split over even w)
 * M = 0 mod 4:  S_u = 2 mod 4 when <u,s>=1  => sum_i A_i K_k(i) >= 2 (C(n,k) - K_k(w)) / M
                 (w = 0 gives nothing)
If max sum A_i (SDP with these constraints) < M for every allowed case, no code of size exactly M exists.
"""
from fractions import Fraction
import sys, time
import split
from sdpcore import binom


def kraw(n, k, x):
    return sum((-1) ** j * binom(x, j) * binom(n - x, k - j) for j in range(k + 1))


def cases(M, n):
    if M % 2 == 1:
        return [('odd', None)]
    if M % 4 == 2:
        return [('2mod4', w) for w in range(0, n + 1, 2)]
    return [('0mod4', w) for w in range(2, n + 1, 2)] + [('0mod4', 0)]


def build(n, d, n1, M, case, w, opts=None, name=None):
    P = split.build(n, d, n1, opts=opts, name=name or 'par_%d_%d_%d_%d_%s_%s' % (n, d, n1, M, case, w))
    # distance distribution forms
    X0 = {}
    for v in P.vars:
        # x_{(a),0}: key ((0..),(a..),(a..)) with one pair all-zero
        if all(s == 0 for s in v[0]) and v[1] == v[2]:
            a = sum(v[1])
            co = 1
            if n1 > 0:
                co = binom(n1, v[1][0]) * binom(n - n1, v[1][1])
            else:
                co = binom(n, a)
            X0.setdefault(a, {})[v] = Fraction(co)
    new = []
    for k in range(0, n + 1):
        if case == 'odd':
            rhs = Fraction(binom(n, k), M)
        elif case == '2mod4':
            rhs = Fraction(2 * (binom(n, k) + kraw(n, k, w)), M)
        else:
            rhs = Fraction(2 * (binom(n, k) - kraw(n, k, w)), M)
        # sum_i A_i K_k(i) - rhs >= 0 ; A_0 = 1
        g0 = Fraction(kraw(n, k, 0)) - rhs
        co = {}
        for a, f in X0.items():
            K = kraw(n, k, a)
            for v, c in f.items():
                co[v] = co.get(v, 0) + c * K
        co = {v: c for v, c in co.items() if c != 0}
        if not co:
            continue
        mx = max([abs(c) for c in co.values()] + [abs(g0)])
        new.append((g0 / mx, {v: c / mx for v, c in co.items()}))
    P.lin = P.lin + new
    return P


if __name__ == '__main__':
    n, d, n1, M = map(int, sys.argv[1:5])
    for case, w in cases(M, n):
        t = time.time()
        P = build(n, d, n1, M, case, w)
        info = P.solve('agents/sdp/runs/par', threads=4)
        print(n, d, n1, M, case, w, info['phase.value'], info['bound'], 'EXCLUDED' if (info['bound'] or 1e18) < M else 'feasible', '%.0fs' % (time.time() - t), flush=True)


def build_split2mod4(n, d, M, w, opts=None, name=None):
    """M = 2 mod 4, |s| = w, split T1 = supp(s) (n1 = w), T2 = rest.
    For u with |u cap T1| = k1 even: S_u = 2 mod 4 => S_u^2 >= 4. Per class (k1,k2):
      sum_{a1,a2} A_{a1,a2} K^{w}_{k1}(a1) K^{n-w}_{k2}(a2) >= 4 C(w,k1) C(n-w,k2) / M   (k1 even)."""
    assert M % 4 == 2 and w % 2 == 0 and 0 < w < n
    n1, n2 = w, n - w
    P = split.build(n, d, n1, opts=opts, name=name or 'psplit_%d_%d_%d_w%d' % (n, d, M, w))
    X0 = {}
    for v in P.vars:
        if v[0] == (0, 0) and v[1] == v[2]:
            a1, a2 = v[1]
            X0[(a1, a2)] = (v, Fraction(binom(n1, a1) * binom(n2, a2)))
    new = []
    for k1 in range(0, n1 + 1, 2):
        for k2 in range(0, n2 + 1):
            rhs = Fraction(4 * binom(n1, k1) * binom(n2, k2), M)
            g0 = Fraction(kraw(n1, k1, 0) * kraw(n2, k2, 0)) - rhs
            co = {}
            for (a1, a2), (v, c) in X0.items():
                K = kraw(n1, k1, a1) * kraw(n2, k2, a2)
                if K:
                    co[v] = co.get(v, 0) + c * K
            if not co:
                continue
            mx = max([abs(c) for c in co.values()] + [abs(g0)])
            new.append((g0 / mx, {v: c / mx for v, c in co.items()}))
    P.lin = P.lin + new
    return P
