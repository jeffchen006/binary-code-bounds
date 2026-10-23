#!/usr/bin/env python3
"""[novaB B1 extension of agents/sdp/check_sdp_cert.py: adds OA rows (case 'oa'), see METHOD_B1.md]
Independent exact checker for Schrijver-SDP dual certificates (binary codes, even d, even-weight codes).

Usage: check_sdp_cert.py cert.json [cert2.json ...]

Each certificate proves: every even-weight binary code C of length n, min distance >= d
(and, if 'M' is given, |C| = M with the stated parity case) satisfies |C| <= bound < target.

Relaxation (Schrijver 2005). For a code C, x(i,j,t) = (1/|C|) * #{(X,Y,Z) in C^3 : d(X,Y)=i, d(X,Z)=j,
|(X+Y) & (X+Z)| = t} / (n!/(t!(i-t)!(j-t)!(n-i-j+t)!)).  Facts used (all standard, Schrijver 2005):
  x depends only on the sorted triple of distances (i, j, i+j-2t);  x(0,0,0) = 1;
  x = 0 if one of the three distances is in {1..d-1} or odd (even-weight code);  0 <= x <= 1;
  x(i,j,t) <= x(i,0,0);  x(i,0,0) + x(j,0,0) <= 1 + x(i,j,t);
  for k = 0..n/2 the matrices  [sum_t beta(n,i,j,k,t) x(i,j,t)]_{i,j=k..n-k}  and
  [sum_t beta(n,i,j,k,t) (x(i+j-2t,0,0) - x(i,j,t))]_{i,j=k..n-k}  are PSD;
  |C| = sum_i C(n,i) x(i,0,0).
Parity rows (A_i = C(n,i) x(i,0,0), A_0 = 1, K_k = Krawtchouk; S_u = sum_c (-1)^{u.c}):
  sum_{|u|=k} S_u^2 = |C| sum_i A_i K_k(i).
  M odd: S_u odd => sum_i A_i K_k(i) >= C(n,k)/M.
  M = 2 mod 4, w = |s|, s = sum of codewords mod 2 (w even since all weights are even):
     S_u = M - 2#{c: u.c = 1} = 2 mod 4 whenever u.s = 0, and #{|u|=k, u.s=0} = (C(n,k)+K_k(w))/2
     => sum_i A_i K_k(i) >= 2 (C(n,k) + K_k(w)) / M.
Certificate: PSD matrices Y = V V^T / 4^K (exact), multipliers z / 2^K >= 0 on valid linear rows.
Weak duality: |C| <= 1 + sum_v c_v x_v + sum <Y, A(x)> + sum z g(x) = K0 + sum_v r_v x_v <= K0 + sum_v max(r_v, 0).
"""
import json, sys, math
from fractions import Fraction


def binom(a, b):
    return math.comb(a, b) if 0 <= b <= a else 0


def beta(n, i, j, k, t):
    s = 0
    for u in range(t, n + 1):           # binom(u,t) = 0 for u < t
        s += (1 if (u - t) % 2 == 0 else -1) * binom(u, t) * binom(n - 2 * k, u - k) * binom(n - k - u, i - u) * binom(n - k - u, j - u)
    assert isinstance(s, int)
    return s


def kraw(n, k, x):
    s = sum((1 if j % 2 == 0 else -1) * binom(x, j) * binom(n - x, k - j) for j in range(k + 1))
    assert isinstance(s, int)
    return s


def check(path):
    C = json.load(open(path))
    pr = C['problem']
    n, d = pr['n'], pr['d']
    M, case, w = pr.get('M'), pr.get('case'), pr.get('w')
    K = C['K']
    assert d % 2 == 0 and d >= 2

    def allowed(a):
        return a == 0 or (a >= d and a % 2 == 0)

    def x(i, j, t):
        """linear form (dict key->coef, key None = constant) of x(i,j,t)"""
        assert 0 <= t <= min(i, j) and i + j - t <= n
        tri = tuple(sorted((i, j, i + j - 2 * t)))
        if not all(allowed(a) for a in tri):
            return {}
        if tri == (0, 0, 0):
            return {None: Fraction(1)}
        return {tri: Fraction(1)}

    def add(acc, form, c):
        for k_, v in form.items():
            acc[k_] = acc.get(k_, 0) + c * v

    L = {}  # Lagrangian form
    # objective
    L[None] = Fraction(1)
    for i in range(1, n + 1):
        add(L, x(i, 0, 0), binom(n, i))
    # PSD blocks
    D2 = Fraction(1, 4 ** K)
    for blk in C['blocks']:
        k, typ, rows, V = blk['k'], blk['typ'], blk['rows'], blk['V']
        assert 0 <= k <= n // 2 and typ in (0, 1)
        assert len(set(rows)) == len(rows) and all(k <= r <= n - k for r in rows)
        assert len(V) == len(rows)
        m = len(rows)
        Y = [[sum(a * b for a, b in zip(V[p], V[q])) for q in range(m)] for p in range(m)]
        for p in range(m):
            for q in range(m):
                if Y[p][q] == 0:
                    continue
                i, j = rows[p], rows[q]
                ent = {}
                for t in range(max(0, i + j - n), min(i, j) + 1):
                    be = beta(n, i, j, k, t)
                    if be == 0:
                        continue
                    if typ == 0:
                        add(ent, x(i, j, t), be)
                    else:
                        add(ent, x(i + j - 2 * t, 0, 0), be)
                        add(ent, x(i, j, t), -be)
                add(L, ent, Y[p][q] * D2)
    # valid linear rows
    def norm(g0, co):
        co = {k_: Fraction(v) for k_, v in co.items() if v != 0}
        g0 = Fraction(g0) + co.pop(None, 0)
        if not co:
            return None
        mx = max([abs(v) for v in co.values()] + [abs(g0)])
        return (g0 / mx, frozenset((k_, v / mx) for k_, v in co.items()))
    valid = set()
    for i in range(n + 1):
        for j in range(n + 1):
            for t in range(max(0, i + j - n), min(i, j) + 1):
                v = x(i, j, t)
                a, b = x(i, 0, 0), x(j, 0, 0)
                if v and None not in v:
                    for f in (a, b):
                        if f and f != v:
                            r = {}
                            add(r, f, 1); add(r, v, -1)
                            valid.add(norm(0, r))
                if i > 0 and j > 0 and a and b:
                    r = {}
                    add(r, a, -1); add(r, b, -1); add(r, v, 1)
                    valid.add(norm(1, r))
    if M is not None and case != 'oa':
        for kk in range(n + 1):
            if case == 'odd':
                assert M % 2 == 1
                rhs = Fraction(binom(n, kk), M)
            elif case == '2mod4':
                assert M % 4 == 2 and w % 2 == 0 and 0 <= w <= n
                rhs = Fraction(2 * (binom(n, kk) + kraw(n, kk, w)), M)
            else:
                raise ValueError(case)
            r = {}
            for i in range(1, n + 1):
                add(r, x(i, 0, 0), binom(n, i) * kraw(n, kk, i))
            valid.add(norm(kraw(n, kk, 0) - rhs, r))
    valid.discard(None)

    def oa_row(blocks):
        # V = span of r vectors with disjoint supports of sizes blocks (dim r); m_k = #nonzero v in V of weight k.
        # coset counts of V^perp are r-bit integers summing to M => sum_{v in V, v != 0} S_v^2 >= rho (2^r - rho),
        # rho = M mod 2^r. Averaging over S_n: sum_k m_k / C(n,k) * sum_{|u|=k} S_u^2 / 1 ... i.e.
        # sum_k m_k / C(n,k) * sum_i A_i K_k(i) >= rho (2^r - rho) / M   (A_0 = 1).
        assert all(isinstance(a, int) and a >= 1 for a in blocks) and sum(blocks) <= n
        r = len(blocks)
        m = {}
        for mask in range(1, 2 ** r):
            wt = sum(blocks[b] for b in range(r) if mask >> b & 1)
            m[wt] = m.get(wt, 0) + 1
        rho = M % 2 ** r
        rhs = Fraction(rho * (2 ** r - rho), M)
        g0 = -rhs
        r_ = {}
        for kk, mk in m.items():
            wgt = Fraction(mk, binom(n, kk))
            g0 += wgt * kraw(n, kk, 0)
            for i in range(1, n + 1):
                add(r_, x(i, 0, 0), wgt * binom(n, i) * kraw(n, kk, i))
        return norm(g0, r_)
    def oa_row_cols(r, mult):
        # V = row space of the r x l generator matrix whose columns are the nonzero c in F2^r with multiplicity mult[c-1].
        # weight of message x = sum of mult[c] over c with x.c odd; dim V = r requires every weight > 0.
        assert len(mult) == 2 ** r - 1 and all(isinstance(a, int) and a >= 0 for a in mult) and sum(mult) <= n
        m = {}
        nz = [(c, a) for c, a in enumerate(mult, 1) if a]
        for xm in range(1, 2 ** r):
            wt = sum(a for c, a in nz if bin(xm & c).count('1') % 2 == 1)
            assert wt > 0, 'V not of dimension r'
            m[wt] = m.get(wt, 0) + 1
        rho = M % 2 ** r
        g0 = -Fraction(rho * (2 ** r - rho), M)
        r_ = {}
        for kk, mk in m.items():
            wgt = Fraction(mk, binom(n, kk))
            g0 += wgt * kraw(n, kk, 0)
            for i in range(1, n + 1):
                add(r_, x(i, 0, 0), wgt * binom(n, i) * kraw(n, kk, i))
        return norm(g0, r_)

    def cap_row(Ap, eps=Fraction(0)):
        # PREMISE A(n-1,d) <= Ap. Shortening at coordinate i: #{c: c_i = b} <= Ap for b = 0,1,
        # so |S_{e_i}| = |M - 2 w_i| <= 2Ap - M, hence sum_i S_{e_i}^2 = M sum_j A_j K_1(j) <= n (2Ap - M)^2.
        assert 2 * Ap >= M
        r_ = {}
        for i in range(1, n + 1):
            add(r_, x(i, 0, 0), -binom(n, i) * kraw(n, 1, i))
        return norm(Fraction(n * (2 * Ap - M) ** 2, M) + eps - kraw(n, 1, 0), r_)
    D1 = Fraction(1, 2 ** K)
    for row in C['lin']:
        z = row['z']
        assert z >= 0
        g0 = Fraction(row['g0'])
        co = {tuple(int(s) for s in key.split(',')): Fraction(v) for key, v in row['coef'].items()}
        nr = norm(g0, co)
        if row.get('oa_blocks'):
            assert case == 'oa'
            assert nr == oa_row(row['oa_blocks']), ('OA row mismatch', row['oa_blocks'])
        elif row.get('oa_cols'):
            assert case == 'oa'
            assert nr == oa_row_cols(row['oa_cols']['r'], row['oa_cols']['mult']), ('OA row mismatch', row['oa_cols'])
        elif row.get('cap') is not None:
            assert case == 'oa' and row['cap'] == pr['premise_A_n1_d_le']
            eps = Fraction(row.get('cap_eps', '0'))
            assert eps >= 0       # g <= cap + eps is weaker than g <= cap, so still valid
            assert nr == cap_row(row['cap'], eps), ('cap row mismatch', row['cap'])
        else:
            assert nr in valid, ('invalid linear row', row)
        zz = z * D1
        L[None] = L.get(None, 0) + zz * g0
        add(L, co, zz)
    bound = L.pop(None, Fraction(0))
    for key, r in L.items():
        if r > 0:
            bound += r        # 0 <= x <= 1
    return pr, bound


if __name__ == '__main__':
    ok = True
    for p in sys.argv[1:]:
        pr, b = check(p)
        target = pr.get('M') if pr.get('M') is not None else pr.get('target')
        verdict = (b < target) if target is not None else None
        print('%s: n=%d d=%d M=%s case=%s premise A(n-1,d)<=%s  bound=%.9f  target=%s  %s' % (
            p, pr['n'], pr['d'], pr.get('M'), pr.get('case'), pr.get('premise_A_n1_d_le'), float(b), target,
            {True: 'OK (bound < target)', False: 'FAIL', None: ''}[verdict]))
        ok = ok and verdict is not False
    print('ALL OK' if ok else 'SOME FAILED')
