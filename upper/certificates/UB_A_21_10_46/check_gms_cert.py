#!/usr/bin/env python3
"""Independent exact checker for 4-point (Gijswijt-Mittelmann-Schrijver 2012) SDP dual certificates
for binary codes with even minimum distance d (all distances even), homogenized form.

Usage:
  check_gms_cert.py cert.json [...]     check certificates
  check_gms_cert.py --selftest          compare every block family with full 2^n-dimensional matrices (n <= 7)

Variables (all are G-averaged code moments; for a code C and a set S of words, x(S) = Pr_g[g(S) subset C],
g uniform in the isometry group of {0,1}^n):
  p          = x({v})                         (so |C| = 2^n p)
  T3 (a,b,c) = x({0,u,v}) with sorted pairwise distances (a,b,c), a = 0 allowed (2-sets: (0,m,m))
  Q  (8 counts) = x({w0,w1,w2,w3}) of a 4-set, canonical column-type vector (see canon4)
All x(S) = 0 when two words of S have distance in {1..d-1} or odd distance. 0 <= x(S) <= p <= 1.
PSD constraints (all derived from M_S(x) of GMS (3)/(4) for S = empty and S = {0,u}, and Schrijver's R, R'):
  'schr' blocks  (k, typ): Schrijver 2005, [sum_t beta^t_{ijk} x(i,j,t)] and [sum_t beta (x(i+j-2t,0,0) - x(i,j,t))]
  'gms'  blocks  (m, k, l): M_{0,u}, |u| = m, Terwilliger(m) x Terwilliger(n-m) block diagonalization
  'gms0' blocks  (m, k, l): M_empty via GMS (39)-(41) (Krawtchouk coefficients), unordered pairs (i' even)
  'gms0inv'      : invariant block of M_empty on span{empty, chi_a}, corner x(empty) = 1
Linear rows (each must be in the generated valid family, after the normalization used here):
  x(T) <= x(2-subset);  x(i,0,0)+x(j,0,0) <= p + x(i,j,t);  Q <= its 3-subsets;  x <= p;
  parity rows (M odd / M = 2 mod 4 with |s| = w) homogenized: sum_i A_i K_k(i) >= c_k / M with A_i = C(n,i) x_i / p.
Every certificate is a claim about codes of size EXACTLY Mt (Mt = problem M, or Mt = target for a 'base'
certificate; a code with more than Mt words contains a subcode of exactly Mt words, so 'no code of size Mt' implies
A(n,d) <= Mt - 1). For such a code p = Mt/2^n and 0 <= x(S) <= p for every S, so
Weak duality: 2^n p <= K0 + sum_v max(r_v, 0) * Mt/2^n  must be < Mt for a contradiction.
"""
import json, sys, math, itertools
from fractions import Fraction
from functools import lru_cache


def C(a, b):
    return math.comb(a, b) if 0 <= b <= a else 0


@lru_cache(None)
def beta(n, i, j, k, t):
    s = 0
    for u in range(t, n + 1):
        s += (1 if (u - t) % 2 == 0 else -1) * C(u, t) * C(n - 2 * k, u - k) * C(n - k - u, i - u) * C(n - k - u, j - u)
    return s


@lru_cache(None)
def K(N, r, x):
    return sum((1 if i % 2 == 0 else -1) * C(x, i) * C(N - x, r - i) for i in range(r + 1))


PERMS = list(itertools.permutations(range(4)))


def canon4(cnt):
    """cnt[b1*4+b2*2+b3] = number of coordinates where (w1,w2,w3) = (b1,b2,b3) and w0 = 0 (columns mod complement)."""
    best = None
    for pm in PERMS:
        c = [0] * 8
        for q in range(8):
            bits = (0, q >> 2 & 1, q >> 1 & 1, q & 1)
            nb = [bits[pm[r]] for r in range(4)]
            if nb[0]:
                nb = [1 - b for b in nb]
            c[nb[1] * 4 + nb[2] * 2 + nb[3]] += cnt[q]
        c = tuple(c)
        if best is None or c < best:
            best = c
    return best


class Model:
    def __init__(self, n, d):
        self.n, self.d = n, d

    def ok(self, a):
        return a == 0 or (a >= self.d and a % 2 == 0)

    def setform(self, cnt):
        return dict(self._setform(tuple(cnt)))

    @lru_cache(None)
    def _setform(self, cnt):
        return tuple(self._setform0(list(cnt)).items())

    def _setform0(self, cnt):
        """cnt: 8 column-type counts of 4 words (w0=0,w1,w2,w3). Returns form {key: coef} for x({w0..w3}) (as a set)."""
        def dist(a, b):
            ma, mb = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)][a], [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)][b]
            mk = tuple((x + y) % 2 for x, y in zip(ma, mb))
            return sum(cnt[q] for q in range(8) if ((q >> 2 & 1) * mk[0] + (q >> 1 & 1) * mk[1] + (q & 1) * mk[2]) % 2)
        ds = []
        for a in range(4):
            if all(dist(a, b) != 0 for b in ds):
                ds.append(a)
        for a, b in itertools.combinations(ds, 2):
            if not self.ok(dist(a, b)):
                return {}
        if len(ds) == 1:
            return {'p': 1}
        if len(ds) == 2:
            a = dist(ds[0], ds[1])
            return {('T', 0, a, a): 1}
        if len(ds) == 3:
            a, b, c = ds
            return {('T',) + tuple(sorted((dist(a, b), dist(a, c), dist(b, c)))): 1}
        return {('Q',) + canon4(cnt): 1}

    def x3(self, i, j, t):
        """x({0,y,z}), |y|=i, |z|=j, |y&z|=t."""
        cnt = [0] * 8
        n = self.n
        cnt[0] = n - i - j + t
        cnt[4] = i - t   # w1 only  -> (1,0,0) with w2 = w1? use words (0, y, z, z): pattern (y,z,z)
        cnt = [0] * 8
        # words w0 = 0, w1 = y, w2 = z, w3 = z  (repeat allowed)
        cnt[0b000] = n - i - j + t
        cnt[0b100] = i - t
        cnt[0b011] = j - t
        cnt[0b111] = t
        return self.setform(cnt)

    def x4(self, m, i, j, t, i2, j2, t2):
        """x({0,u,v,w}), |u| = m, v: i ones in supp(u), i2 outside; w: j, j2; overlaps t, t2."""
        n = self.n
        cnt = [0] * 8
        cnt[0b100] = m - i - j + t
        cnt[0b110] = i - t
        cnt[0b101] = j - t
        cnt[0b111] = t
        mm = n - m
        cnt[0b000] = mm - i2 - j2 + t2
        cnt[0b010] = i2 - t2
        cnt[0b001] = j2 - t2
        cnt[0b011] = t2
        return self.setform(cnt)

    @lru_cache(None)
    def xlam(self, lam):
        """ordered pairs ((v,w),(v',w')); coordinate class t=(t1,t2,t3) = (v+w, v'+w', w+w'); lam[t1*4+t2*2+t3].
        words: w = 0, v = t1, w' = t3, v' = t2+t3."""
        cnt = [0] * 8
        for t in range(8):
            t1, t2, t3 = t >> 2 & 1, t >> 1 & 1, t & 1
            cnt[t1 * 4 + t3 * 2 + (t2 ^ t3)] += lam[t]
        return self.setform(cnt)

    # ---- block entry forms ----
    def schr_entry(self, k, typ, i, j):
        n = self.n
        f = {}
        for t in range(max(0, i + j - n), min(i, j) + 1):
            b = beta(n, i, j, k, t)
            if b == 0:
                continue
            if typ == 0:
                add(f, self.x3(i, j, t), b)
            else:
                add(f, self.x3(i + j - 2 * t, 0, 0), b)
                add(f, self.x3(i, j, t), -b)
        return f

    def gms_entry(self, m, k, l, r, s):
        n = self.n; mm = n - m
        (i, i2), (j, j2) = r, s
        f = {}
        for t in range(max(0, i + j - m), min(i, j) + 1):
            b1 = beta(m, i, j, k, t)
            if b1 == 0:
                continue
            for t2 in range(max(0, i2 + j2 - mm), min(i2, j2) + 1):
                b2 = beta(mm, i2, j2, l, t2)
                if b2:
                    add(f, self.x4(m, i, j, t, i2, j2, t2), b1 * b2)
        return f

    def theta(self, al, be):
        key = (al, be)
        if not hasattr(self, '_th'):
            self._th = {}
        if key in self._th:
            return self._th[key]
        lp = tuple(al[q] + be[q] for q in range(4))
        f = {}
        for rr in itertools.product(*[range(lp[q] + 1) for q in range(4)]):
            co = 1
            for q in range(4):
                co *= K(lp[q], rr[q], be[q])
                if co == 0:
                    break
            if co == 0:
                continue
            lam = [0] * 8
            for q, (t1, t2) in enumerate([(0, 0), (1, 0), (0, 1), (1, 1)]):
                lam[t1 * 4 + t2 * 2 + 1] = rr[q]
                lam[t1 * 4 + t2 * 2] = lp[q] - rr[q]
            add(f, self.xlam(tuple(lam)), co)
        self._th[key] = f
        return f

    def gms0_entry(self, m, k, l, r, s):
        n = self.n; mm = n - m
        (ia, ib), (ja, jb) = r, s
        f = {}
        for ta in range(max(0, ia + ja - m), min(ia, ja) + 1):
            b1 = beta(m, ia, ja, k, ta)
            if b1 == 0:
                continue
            al = (m - ia - ja + ta, ia - ta, ja - ta, ta)        # P order 00,10,01,11 ; c = (t1,t2)
            for tb in range(max(0, ib + jb - mm), min(ib, jb) + 1):
                b2 = beta(mm, ib, jb, l, tb)
                if b2 == 0:
                    continue
                be = (mm - ib - jb + tb, ib - tb, jb - tb, tb)
                add(f, self.theta(al, be), b1 * b2)
        return f

    def inv_entry(self, r, s):
        """invariant block rows 'e' (empty) and pair distances a; unscaled: (e,e) = 1, (e,a) = 2^n C(n,a) x2(a),
        (a,b) = 2^n sum_{lambda: d(v,w)=a, d(v',w')=b} multinomial(lambda) x(lambda)."""
        n = self.n
        f = {}
        if r == 'e' and s == 'e':
            return {None: 1}
        if r == 'e' or s == 'e':
            a = s if r == 'e' else r
            add(f, self.x3(a, 0, 0) if a > 0 else {'p': 1}, 2 ** n * C(n, a))
            return f
        a, b = r, s
        # lambda over 8 classes with d(v,w) = #t1 = a, d(v',w') = #t2 = b
        for lp in itertools.product(range(n + 1), repeat=3):
            if sum(lp) > n:
                continue
            lpf = (n - sum(lp),) + lp
            if lpf[1] + lpf[3] != a or lpf[2] + lpf[3] != b:
                continue
            for rr in itertools.product(*[range(lpf[q] + 1) for q in range(4)]):
                lam = [0] * 8
                for q, (t1, t2) in enumerate([(0, 0), (1, 0), (0, 1), (1, 1)]):
                    lam[t1 * 4 + t2 * 2 + 1] = rr[q]
                    lam[t1 * 4 + t2 * 2] = lpf[q] - rr[q]
                mult = math.factorial(n)
                for q in lam:
                    mult //= math.factorial(q)
                add(f, self.xlam(tuple(lam)), 2 ** n * mult)
        return f

    def block_entry(self, meta, r, s):
        kd = meta['kind']
        if kd == 'schr':
            return self.schr_entry(meta['k'], meta['typ'], r, s)
        if kd == 'gms':
            return self.gms_entry(meta['m'], meta['k'], meta['l'], tuple(r), tuple(s))
        if kd == 'gms0':
            assert meta['m'] % 1 == 0 and not (meta['m'] == self.n and meta['k'] == 0 and meta['l'] == 0)
            for q in (r, s):
                assert q[1] % 2 == 0, 'gms0 rows must have i2 even (tau-symmetric part)'
            return self.gms0_entry(meta['m'], meta['k'], meta['l'], tuple(r), tuple(s))
        if kd == 'gms0inv':
            return self.inv_entry(r, s)
        raise ValueError(kd)

    def valid_rows(self, M=None, case=None, w=None):
        """normalized-then-homogenized valid linear rows (frozensets of (key, coef))."""
        n = self.n
        rows = set()
        def norm_h(g0, co):
            co = {k: Fraction(v) for k, v in co.items() if v != 0}
            g0 = Fraction(g0) + co.pop('p', 0)      # in normalized coordinates p <-> constant 1
            if not co:
                return None
            mx = max([abs(v) for v in co.values()] + [abs(g0)])
            h = {k: v / mx for k, v in co.items()}
            if g0 != 0:
                h['p'] = h.get('p', 0) + g0 / mx
            return frozenset(h.items())
        for i in range(n + 1):
            for j in range(n + 1):
                for t in range(max(0, i + j - n), min(i, j) + 1):
                    v = self.x3(i, j, t); a = self.x3(i, 0, 0); b = self.x3(j, 0, 0)
                    if v and 'p' not in v:
                        for f in (a, b):
                            if f and f != v:
                                r = {}; add(r, f, 1); add(r, v, -1); rows.add(norm_h(0, r))
                    if i > 0 and j > 0 and a and b:
                        r = {}; add(r, a, -1); add(r, b, -1); add(r, v, 1); rows.add(norm_h(1, r))
        if M is not None:
            for kk in range(n + 1):
                if case == 'odd':
                    assert M % 2 == 1
                    rhs = Fraction(C(n, kk), M)
                elif case == '2mod4':
                    assert M % 4 == 2 and w % 2 == 0 and 0 <= w <= n
                    rhs = Fraction(2 * (C(n, kk) + K(n, kk, w)), M)
                else:
                    raise ValueError(case)
                r = {}
                for i in range(1, n + 1):
                    add(r, self.x3(i, 0, 0), C(n, i) * K(n, kk, i))
                rows.add(norm_h(K(n, kk, 0) - rhs, r))
        rows.discard(None)
        return rows


def add(f, g, c):
    for k, v in g.items():
        f[k] = f.get(k, 0) + c * v


def parse_key(s):
    if s == 'p':
        return 'p'
    p = s.split(',')
    if p[0] == 'T':
        return ('T',) + tuple(int(x) for x in p[1:])
    if p[0] == 'Q':
        return ('Q',) + tuple(int(x) for x in p[1:])
    raise ValueError(s)


def check(path):
    J = json.load(open(path))
    pr = J['problem']; n, d = pr['n'], pr['d']; M, case, w = pr.get('M'), pr.get('case'), pr.get('w')
    Kd = J['K']
    md = Model(n, d)
    L = {'p': Fraction(2 ** n)}         # objective 2^n p
    D2 = Fraction(1, 4 ** Kd)
    for blk in J['blocks']:
        meta = blk['meta']; rows = [tuple(r) if isinstance(r, list) else r for r in meta['rows']]
        V = blk['V']; m = len(rows); assert len(V) == m and len(set(rows)) == m
        for p_ in range(m):
            for q_ in range(m):
                y = sum(a * b for a, b in zip(V[p_], V[q_]))
                if y == 0:
                    continue
                add(L, md.block_entry(meta, rows[p_], rows[q_]), y * D2)
    valid = md.valid_rows(M, case, w)
    # Q <= its 3-subsets, x <= p : generated on demand
    D1 = Fraction(1, 2 ** Kd)
    for row in J['lin']:
        z = row['z']; assert z >= 0
        co = {parse_key(k): Fraction(v) for k, v in row['coef'].items()}
        assert Fraction(row['g0']) == 0
        fs = frozenset(co.items())
        ok = fs in valid
        if not ok:
            # x(Q) <= x(3-subset)  or  x <= p
            ks = list(co.items())
            if len(ks) == 2:
                (k1, c1), (k2, c2) = ks
                if c1 == -c2:
                    big, small = (k1, k2) if c1 > 0 else (k2, k1)
                    if big == 'p' and small != 'p':
                        ok = True
                    elif big[0] == 'T' and small[0] == 'Q':
                        cnt = small[1:]
                        subs = set()
                        for drop in range(4):
                            keep = [a for a in range(4) if a != drop]
                            # distances among kept words
                            def dist(a, b):
                                ma = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
                                mk = tuple((x + y) % 2 for x, y in zip(ma[a], ma[b]))
                                return sum(cnt[q] for q in range(8) if ((q >> 2 & 1) * mk[0] + (q >> 1 & 1) * mk[1] + (q & 1) * mk[2]) % 2)
                            subs.add(('T',) + tuple(sorted(dist(a, b) for a, b in itertools.combinations(keep, 2))))
                        ok = big in subs
        assert ok, ('invalid linear row', row)
        add(L, co, z * D1)
    Mt = M if M is not None else pr['target']
    ub = Fraction(Mt, 2 ** n)
    bound = L.pop(None, Fraction(0))
    for k, r in L.items():
        if r > 0:
            bound += r * ub
    return pr, bound


def selftest():
    import numpy as np
    rng = np.random.default_rng(5)
    for (n, d) in [(4, 2), (5, 2), (6, 2), (6, 4), (7, 4)]:
        md = Model(n, d)
        vals = {}
        def val(f):
            s = 0.0
            for k, c in f.items():
                if k is None:
                    s += c
                else:
                    if k not in vals:
                        vals[k] = rng.normal()
                    s += c * vals[k]
            return s
        pc = lambda z: bin(z).count('1')
        def xset(ws):
            ws = list(dict.fromkeys(ws))
            for a, b in itertools.combinations(ws, 2):
                if not md.ok(pc(a ^ b)):
                    return 0.0
            ws = ws + [ws[-1]] * (4 - len(ws))
            w0 = ws[0]
            cnt = [0] * 8
            for bit in range(n):
                b = [(x ^ w0) >> bit & 1 for x in ws]
                cnt[b[1] * 4 + b[2] * 2 + b[3]] += 1
            return val(md.setform(cnt))
        N = 2 ** n
        # full M_empty (ordered pairs + empty row) vs gms0 + inv blocks
        rows = [(v, w) for v in range(N) for w in range(N) if md.ok(pc(v ^ w))]
        F = np.zeros((len(rows) + 1, len(rows) + 1))
        F[0, 0] = 1.0
        for a, (v, w) in enumerate(rows):
            F[0, a + 1] = F[a + 1, 0] = xset([v, w]) if v != w else val({'p': 1})
            for b, (v2, w2) in enumerate(rows):
                F[a + 1, b + 1] = xset([v, w, v2, w2])
        full = np.linalg.eigvalsh(F)
        ev = []
        Dset = [a for a in range(n + 1) if md.ok(a)]
        for m in range(n + 1):
            for k in range(m // 2 + 1):
                for l in range((n - m) // 2 + 1):
                    if m == n and k == 0 and l == 0:
                        continue
                    rr = [(i, i2) for i in range(k, m - k + 1) for i2 in range(l, n - m - l + 1) if i2 % 2 == 0 and md.ok(i + i2)]
                    if not rr:
                        continue
                    A = np.array([[val(md.gms0_entry(m, k, l, r, s)) / math.sqrt(C(m - 2 * k, r[0] - k) * C(n - m - 2 * l, r[1] - l) * C(m - 2 * k, s[0] - k) * C(n - m - 2 * l, s[1] - l)) for s in rr] for r in rr])
                    ev.append(np.linalg.eigvalsh(A))
        rr = ['e'] + Dset
        sc = [1.0] + [1 / math.sqrt(2 ** n * C(n, a)) for a in Dset]
        A = np.array([[val(md.inv_entry(r, s)) * sc[p_] * sc[q_] for q_, s in enumerate(rr)] for p_, r in enumerate(rr)])
        ev.append(np.linalg.eigvalsh(A))
        be = np.concatenate(ev)
        okb = all(np.min(np.abs(full - x)) < 1e-7 * max(1, abs(x)) for x in be)
        miss = [x for x in full if np.min(np.abs(be - x)) > 1e-7 * max(1, abs(x))]
        print('M_empty n=%d d=%d: block eigs in full: %s ; full eigs missing from blocks: %d (max |.| %.1e, must be ~0)' % (n, d, okb, len(miss), max([abs(x) for x in miss] + [0])))
        # M_{0,u} for each even m >= d
        for m in range(d, n + 1, 2):
            u = (1 << m) - 1
            Np = [v for v in range(N) if md.ok(pc(v)) and pc(v) > 0 and md.ok(pc(u ^ v))]
            Fm = np.array([[xset([0, u, v, w]) for w in Np] for v in Np])
            fe = np.linalg.eigvalsh(Fm)
            ev = []
            for k in range(m // 2 + 1):
                for l in range((n - m) // 2 + 1):
                    rr = [(i, i2) for i in range(k, m - k + 1) for i2 in range(l, n - m - l + 1) if md.ok(i + i2) and i + i2 > 0 and md.ok(m - i + i2)]
                    if not rr:
                        continue
                    A = np.array([[val(md.gms_entry(m, k, l, r, s)) / math.sqrt(C(m - 2 * k, r[0] - k) * C(n - m - 2 * l, r[1] - l) * C(m - 2 * k, s[0] - k) * C(n - m - 2 * l, s[1] - l)) for s in rr] for r in rr])
                    ev.append(np.linalg.eigvalsh(A))
            be = np.concatenate(ev)
            okb = all(np.min(np.abs(fe - x)) < 1e-7 * max(1, abs(x)) for x in be)
            okf = all(np.min(np.abs(be - x)) < 1e-7 * max(1, abs(x)) for x in fe)
            print('  M_{0,u} n=%d d=%d m=%d: blocks<->full %s %s' % (n, d, m, okb, okf))


if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        selftest(); sys.exit(0)
    allok = True
    for p in sys.argv[1:]:
        pr, b = check(p)
        T = pr.get('target') if pr.get('M') is None else pr['M']
        ok = b < T
        allok &= ok
        print('%s: n=%d d=%d M=%s case=%s w=%s  bound(2^n p)=%.9f  target=%s  %s' % (p, pr['n'], pr['d'], pr.get('M'), pr.get('case'), pr.get('w'), float(b), T, 'OK' if ok else 'FAIL'))
    print('ALL OK' if allok else 'SOME FAILED')
