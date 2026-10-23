"""Light GMS (Gijswijt-Mittelmann-Schrijver 2012) 4-point SDP, in Schrijver-normalized coordinates.

Even-weight codes, d even. For a code C and a set T of words, X(T) := 2^{n-1} * #{T' in orbit(T): T' subset C} / (|C| |orbit(T)|)
(orbit under the even-weight-preserving isometry group). X(single) = 1, X({0,u}) = Schrijver x(m,0,0), X(3-set) = Schrijver
x(i,j,t), X(4-set) = new variable.

Constraints: Schrijver's (split.build, n1=0) + for every even m in [d,n], S={0,u}, |u|=m:
   M_S = [X(S u {v,w})]_{v,w in N'}  PSD, N' = {v: d(0,v) in [d,n], d(u,v) in {0}u[d,n]},
 block-diagonalized via Terwilliger(m) x Terwilliger(n-m) (same beta formulas as the split SDP).
 + 0 <= X4 <= X(3-subsets).
"""
import itertools, math, sys, time
from fractions import Fraction
from functools import lru_cache
import mpmath
import split
from sdpcore import binom

PERMS4 = list(itertools.permutations(range(4)))


def _pat_index(bits):
    # bits: tuple of 4 bits (w0,w1,w2,w3) for one coordinate; normalize w0 = 0 by complement
    if bits[0]:
        bits = tuple(1 - b for b in bits)
    return bits[1] * 4 + bits[2] * 2 + bits[3]  # 0..7


# precompute pattern permutation tables: for perm p, pattern index q -> new index
_PT = []
for p in PERMS4:
    tab = []
    for q in range(8):
        bits = (0, (q >> 2) & 1, (q >> 1) & 1, q & 1)
        nb = tuple(bits[p[r]] for r in range(4))
        tab.append(_pat_index(nb))
    _PT.append(tab)


def canon4(counts):
    """counts: tuple of 8 counts indexed by pattern (u,v,w) bits with word0 = 0. Return canonical tuple."""
    best = None
    for tab in _PT:
        c = [0] * 8
        for q in range(8):
            c[tab[q]] += counts[q]
        c = tuple(c)
        if best is None or c < best:
            best = c
    return best


def dists4(c):
    """pairwise distances among words (0,u,v,w) from pattern counts"""
    d = {}
    def wt(mask):
        return sum(c[q] for q in range(8) if ((q >> 2 & 1) * mask[0] + (q >> 1 & 1) * mask[1] + (q & 1) * mask[2]) % 2)
    d['0u'] = wt((1, 0, 0)); d['0v'] = wt((0, 1, 0)); d['0w'] = wt((0, 0, 1))
    d['uv'] = wt((1, 1, 0)); d['uw'] = wt((1, 0, 1)); d['vw'] = wt((0, 1, 1))
    return d


def build(n, d, opts=None, ms=None, name=None, quad_lin=True):
    assert d % 2 == 0
    opts = dict(opts or {})
    P = split.build(n, d, 0, opts=opts, name=name or 'gms_%d_%d' % (n, d))

    def okd(a):
        return a == 0 or (a >= d and a % 2 == 0)

    def key3(a, b, c):
        k = tuple(sorted((a, b, c)))
        if not all(okd(x) for x in k):
            return None
        if k == (0, 0, 0):
            return 'one'
        kk = tuple((x,) for x in k)
        return kk if kk in P.vidx else P.var(kk)

    q4 = {}

    def setkey(m, i, j, t, i2, j2, t2):
        """X of set {0,u,v,w}: u = 1^m 0^(n-m); v: i ones inside, i2 outside; w: j, j2; overlaps t, t2."""
        c = [0] * 8
        # inside supp(u): pattern (1, v, w)
        c[4 + 0] = m - i - j + t       # (1,0,0)
        c[4 + 2] = i - t               # (1,1,0)
        c[4 + 1] = j - t               # (1,0,1)
        c[4 + 3] = t                   # (1,1,1)
        mm = n - m
        c[0] = mm - i2 - j2 + t2
        c[2] = i2 - t2
        c[1] = j2 - t2
        c[3] = t2
        assert min(c) >= 0
        dd = dists4(c)
        # distinct words: 0,u distinct (m>=d). v==w iff vw==0 ; v==u iff uv==0 ; w==u iff uw==0 ; v==0/w==0 excluded
        for k_ in ('0u', '0v', '0w', 'uv', 'uw', 'vw'):
            if not okd(dd[k_]):
                return None
        words = ['0', 'u', 'v', 'w']
        same = {('v', 'w'): dd['vw'] == 0, ('u', 'v'): dd['uv'] == 0, ('u', 'w'): dd['uw'] == 0,
                ('0', 'v'): dd['0v'] == 0, ('0', 'w'): dd['0w'] == 0}
        distinct = ['0', 'u']
        if not (same[('u', 'v')] or same[('0', 'v')]):
            distinct.append('v')
        wdup = same[('u', 'w')] or same[('0', 'w')] or (same[('v', 'w')])
        if not wdup:
            distinct.append('w')
        if len(distinct) == 2:
            return key3(m, 0, m) if m > 0 else 'one'      # x(m,0,0): sorted (0,m,m)
        if len(distinct) == 3:
            x3 = distinct[2]
            a = dd['0u']; b = dd['0' + x3]; cc = dd['u' + x3]
            return key3(a, b, cc)
        cc = canon4(tuple(c))
        k = ('Q',) + cc
        if k not in P.vidx:
            P.var(k)
            q4[k] = c
        return k

    def okrow(i, i2, m):
        return okd(i + i2) and (i + i2) > 0 and okd(m - i + i2)

    t0 = time.time()
    for m in (ms or range(d, n + 1, 2)):
        mm = n - m
        for k in range(0, m // 2 + 1):
            for l in range(0, mm // 2 + 1):
                rows = [(i, i2) for i in range(k, m - k + 1) for i2 in range(l, mm - l + 1) if okrow(i, i2, m)]
                if not rows:
                    continue
                ent = {}
                for a, r in enumerate(rows):
                    for s in rows[a:]:
                        i, i2 = r; j, j2 = s
                        for t in range(max(0, i + j - m), min(i, j) + 1):
                            b1 = split.beta(m, i, j, k, t)
                            if b1 == 0:
                                continue
                            for t2 in range(max(0, i2 + j2 - mm), min(i2, j2) + 1):
                                b2 = split.beta(mm, i2, j2, l, t2)
                                if b2 == 0:
                                    continue
                                v = setkey(m, i, j, t, i2, j2, t2)
                                if v is None:
                                    continue
                                kk = None if v == 'one' else v
                                ent.setdefault(kk, {})
                                ent[kk][(r, s)] = ent[kk].get((r, s), 0) + b1 * b2
                used = set()
                for kk, mat in ent.items():
                    for (r, s), val in mat.items():
                        if val:
                            used.add(r); used.add(s)
                keep = [r for r in rows if r in used]
                if not keep:
                    continue
                pos = {r: p for p, r in enumerate(keep)}
                ent2 = {}
                for kk, mat in ent.items():
                    m2 = {(pos[r], pos[s]): Fraction(val) for (r, s), val in mat.items() if val}
                    if m2:
                        ent2[kk] = m2
                scl = [1 / mpmath.sqrt(binom(m - 2 * k, i - k) * binom(mm - 2 * l, i2 - l)) for (i, i2) in keep]
                P.add_block(len(keep), ent2, scl)
                P.blocks[-1]['meta'] = dict(kind='gms', m=m, k=k, l=l, rows=keep)
    # linear: 0 <= X4 <= X3 for each 3-subset
    if quad_lin:
        for k4, c in q4.items():
            # the 4 words 0,u,v,w ; 3-subsets: drop one word
            dd = dists4(c)
            trip = [(dd['uv'], dd['uw'], dd['vw']), (dd['0v'], dd['0w'], dd['vw']), (dd['0u'], dd['0w'], dd['uw']), (dd['0u'], dd['0v'], dd['uv'])]
            for (a, b, cc) in trip:
                k3 = key3(a, b, cc)
                if k3 is None or k3 == 'one':
                    continue
                P.lin.append((Fraction(0), {k3: Fraction(1), k4: Fraction(-1)}))
    print('gms build n=%d d=%d: vars=%d (4-pt %d) blocks=%d maxblk=%d lin=%d %.0fs' % (
        n, d, len(P.vars), len(q4), len(P.blocks), max(b['size'] for b in P.blocks), len(P.lin), time.time() - t0), flush=True)
    return P


if __name__ == '__main__':
    n, d = int(sys.argv[1]), int(sys.argv[2])
    P = build(n, d)
    t = time.time()
    info = P.solve('agents/sdp/runs/gms', threads=4)
    print(n, d, info['phase.value'], info['bound'], '%.0fs' % (time.time() - t))
