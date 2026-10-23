"""M_empty block of the GMS 4-point SDP (Gijswijt-Mittelmann-Schrijver 2012, Sec. 5), normalized coordinates.

M_empty: rows/cols = {empty} u {singletons} u {pairs {v,w}, d(v,w) in D}; entry (S',S'') = x(S' u S'') (0 if not a code).
Normalized: X(T) = x(T)/p for T nonempty (p = |C|/2^n), X(empty) = 2^n/M  (here: 2^n/T for a target T; a larger
corner keeps PSD, so corner 2^n/T is valid for every code with |C| >= T).

Ordered-pair matrix X^ord_{(v,w),(v',w')} = X({v,w,v',w'}) = sum_lambda X(lambda) B_lambda, lambda in Lambda^n_T,
class of coordinate = (v+w, v'+w', w+w').  Block diagonalization (GMS (39)-(41)):
 psi_{m,k,l}(B_lambda) = sum_{alpha in Lambda^m_P, beta in Lambda^{n-m}_P, alpha+beta = lambda'}
                          prod_{c in P} K^{lambda'(c)}_{lambda(c1)}(beta(c)) * Gamma_{alpha,k} (x) Gamma_{beta,l}
 with Gamma_{alpha,k} = Schrijver beta^{alpha(11)}_{i_alpha,j_alpha,k}(m) (up to diagonal scaling).
Block (n,0,0) is replaced by the invariant block on span{empty, chi_a (a = pair distance)}.
"""
import itertools, math, sys, time
from fractions import Fraction
from functools import lru_cache
import mpmath
import split, gms
from sdpcore import binom

P4 = [(0, 0), (1, 0), (0, 1), (1, 1)]          # c = (t1, t2)


@lru_cache(None)
def kraw(N, r, x):
    return sum((1 if i % 2 == 0 else -1) * binom(x, i) * binom(N - x, r - i) for i in range(r + 1))


def add_mempty(P, n, d, T, key3, setkey_words, extra_corner=None, debug_ordered=False):
    """Add M_empty blocks to problem P. key3(a,b,c) -> var key for 3-set by distances; setkey_words(list of words as
    bit-count patterns) handled inside. T = target size (corner 2^n/T)."""
    def okd(a):
        return a == 0 or (a >= d and a % 2 == 0)
    Dset = [a for a in range(n + 1) if okd(a)]

    # ---- variable for type lambda: counts over classes t=(t1,t2,t3) -> words w=0, v=t1, w'=t3, v'=t2+t3
    lam_cache = {}

    def lam_var(lam):
        """lam: tuple of 8 counts indexed by t = t1*4 + t2*2 + t3. Returns var key / 'one' / None."""
        if lam in lam_cache:
            return lam_cache[lam]
        # pattern counts for words (w, v, w', v') with w = 0 -> use gms 4-word pattern (u,v,w) bits w.r.t. word0 = w
        # words: W0 = w = 0, W1 = v (bit t1), W2 = w' (bit t3), W3 = v' (bit t2^t3)
        c = [0] * 8
        for t in range(8):
            if lam[t] == 0:
                continue
            t1, t2, t3 = (t >> 2) & 1, (t >> 1) & 1, t & 1
            b1, b2, b3 = t1, t3, t2 ^ t3
            c[b1 * 4 + b2 * 2 + b3] += lam[t]
        res = gms_setvar(c)
        lam_cache[lam] = res
        return res

    def gms_setvar(c):
        # words 0, A, B, Cc with pattern counts c[(a,b,cc)]; compute distinct set and its key
        def wt(mask):
            return sum(c[q] for q in range(8) if (((q >> 2) & 1) * mask[0] + ((q >> 1) & 1) * mask[1] + (q & 1) * mask[2]) % 2)
        names = ['0', 'A', 'B', 'C']
        masks = {'0': (0, 0, 0), 'A': (1, 0, 0), 'B': (0, 1, 0), 'C': (0, 0, 1)}
        def dist(x, y):
            m1, m2 = masks[x], masks[y]
            return wt(tuple((a + b) % 2 for a, b in zip(m1, m2)))
        distinct = []
        for x in names:
            if all(dist(x, y) != 0 for y in distinct):
                distinct.append(x)
        for x, y in itertools.combinations(distinct, 2):
            if not okd(dist(x, y)):
                return None
        if len(distinct) == 1:
            return 'one'
        if len(distinct) == 2:
            a = dist(distinct[0], distinct[1])
            return key3(0, a, a)
        if len(distinct) == 3:
            x, y, z = distinct
            return key3(dist(x, y), dist(x, z), dist(y, z))
        # 4 distinct: canonical 4-set type; re-express with word0 = '0'
        k = ('Q',) + gms.canon4(tuple(c))
        if k not in P.vidx:
            P.var(k)
        return k

    # lists of nonzero lambda for each lambda' (P-counts)
    lp_cache = {}

    def lam_list(lp):
        if lp in lp_cache:
            return lp_cache[lp]
        out = []
        for r in itertools.product(*[range(lp[ci] + 1) for ci in range(4)]):
            lam = [0] * 8
            for ci, (t1, t2) in enumerate(P4):
                lam[t1 * 4 + t2 * 2 + 1] = r[ci]
                lam[t1 * 4 + t2 * 2 + 0] = lp[ci] - r[ci]
            v = lam_var(tuple(lam))
            if v is not None:
                out.append((r, v))
        lp_cache[lp] = out
        return out

    nblk0 = len(P.blocks)
    t0 = time.time()
    for m in range(0, n + 1):
        mm = n - m
        # alpha: (a00,a10,a01,a11) on m coords; i = a10+a11, j = a01+a11, t = a11
        alphas = []
        for i in range(m + 1):
            for j in range(m + 1):
                for t in range(max(0, i + j - m), min(i, j) + 1):
                    alphas.append((m - i - j + t, i - t, j - t, t))
        betas = []
        for i in range(mm + 1):
            if i % 2 and not debug_ordered:
                continue        # tau-symmetric rows: i' even
            for j in range(0, mm + 1, 1 if debug_ordered else 2):
                for t in range(max(0, i + j - mm), min(i, j) + 1):
                    betas.append((mm - i - j + t, i - t, j - t, t))
        # theta(alpha,beta) for allowed rows/cols
        theta = {}
        for al in alphas:
            ia, ja = al[1] + al[3], al[2] + al[3]
            for be in betas:
                ib, jb = be[1] + be[3], be[2] + be[3]
                if not (okd(ia + ib) and okd(ja + jb)):
                    continue
                lp = tuple(al[q] + be[q] for q in range(4))
                f = {}
                for r, v in lam_list(lp):
                    coef = 1
                    for ci in range(4):
                        coef *= kraw(lp[ci], r[ci], be[ci])
                        if coef == 0:
                            break
                    if coef:
                        kk = None if v == 'one' else v
                        f[kk] = f.get(kk, 0) + coef
                f = {a: b for a, b in f.items() if b}
                if f:
                    theta[(al, be)] = f
        # distribute to blocks (k,l)
        for k in range(0, m // 2 + 1):
            for l in range(0, mm // 2 + 1):
                if m == n and k == 0 and l == 0 and not debug_ordered:
                    continue
                ent = {}
                for (al, be), f in theta.items():
                    ia, ja, ta = al[1] + al[3], al[2] + al[3], al[3]
                    ib, jb, tb = be[1] + be[3], be[2] + be[3], be[3]
                    if not (k <= ia <= m - k and k <= ja <= m - k and l <= ib <= mm - l and l <= jb <= mm - l):
                        continue
                    r, s = (ia, ib), (ja, jb)
                    if s < r:
                        continue          # use upper triangle (matrix is symmetric)
                    b1 = split.beta(m, ia, ja, k, ta)
                    if b1 == 0:
                        continue
                    b2 = split.beta(mm, ib, jb, l, tb)
                    if b2 == 0:
                        continue
                    for kk, cf in f.items():
                        ent.setdefault(kk, {})
                        ent[kk][(r, s)] = ent[kk].get((r, s), 0) + b1 * b2 * cf
                used = set()
                for kk, mat in ent.items():
                    for (r, s), val in mat.items():
                        if val:
                            used.add(r); used.add(s)
                keep = sorted(used)
                if not keep:
                    continue
                pos = {r: p for p, r in enumerate(keep)}
                ent2 = {}
                for kk, mat in ent.items():
                    m2 = {}
                    for (r, s), val in mat.items():
                        if val:
                            a, b = pos[r], pos[s]
                            if a > b:
                                a, b = b, a
                            m2[(a, b)] = m2.get((a, b), 0) + Fraction(val)
                    m2 = {q: v for q, v in m2.items() if v}
                    if m2:
                        ent2[kk] = m2
                scl = [1 / mpmath.sqrt(binom(m - 2 * k, i - k) * binom(mm - 2 * l, i2 - l)) for (i, i2) in keep]
                P.add_block(len(keep), ent2, scl)
                P.blocks[-1]['meta'] = dict(kind='gms0', m=m, k=k, l=l, rows=keep)
    if debug_ordered:
        return P
    # invariant block: rows empty + a in Dset (pair distance, a=0 -> singletons)
    # U(e,e) = 2^n/T ; U(e,a) = 2^n C(n,a) X2(a) ; U(a,b) = 2^n sum_{lambda: d(v,w)=a, d(v',w')=b} multinom(lambda) X(lambda)
    ent = {}
    def acc(kk, r, s, val):
        if val == 0:
            return
        kk = None if kk == 'one' else kk
        ent.setdefault(kk, {})
        ent[kk][(r, s)] = ent[kk].get((r, s), 0) + val
    rows = ['e'] + Dset
    acc('one', 0, 0, Fraction(2 ** n, T))
    for p, a in enumerate(Dset):
        v = 'one' if a == 0 else key3(0, a, a)
        if v is not None:
            acc(v, 0, p + 1, 2 ** n * binom(n, a))
    # enumerate lambda' with i_{lambda'} = a, j = b and then lambda
    for lp in itertools.product(range(n + 1), repeat=3):
        if sum(lp) > n:
            continue
        lpf = (n - sum(lp),) + lp          # counts for P4 order (00,10,01,11)
        a = lpf[1] + lpf[3]; b = lpf[2] + lpf[3]
        if a not in Dset or b not in Dset or b < a:
            continue
        pa, pb = Dset.index(a) + 1, Dset.index(b) + 1
        for r, v in lam_list(lpf):
            lam = [0] * 8
            for ci, (t1, t2) in enumerate(P4):
                lam[t1 * 4 + t2 * 2 + 1] = r[ci]
                lam[t1 * 4 + t2 * 2] = lpf[ci] - r[ci]
            mult = math.factorial(n)
            for q in lam:
                mult //= math.factorial(q)
            acc(v, pa, pb, 2 ** n * mult)
    ent2 = {kk: {q: Fraction(v) for q, v in mat.items() if v} for kk, mat in ent.items()}
    ent2 = {kk: m for kk, m in ent2.items() if m}
    scl = [mpmath.mpf(1)] + [1 / mpmath.sqrt(mpmath.mpf(2) ** n * binom(n, a)) for a in Dset]
    P.add_block(len(rows), ent2, scl)
    P.blocks[-1]['meta'] = dict(kind='gms0inv', rows=rows)
    print('M_empty: %d blocks, max size %d, %.0fs' % (len(P.blocks) - nblk0, max(b['size'] for b in P.blocks[nblk0:]), time.time() - t0), flush=True)
    return P


def build_full(n, d, T, opts=None, name=None, light=True, quad_lin=True):
    """Schrijver + (light GMS M_S blocks if light) + M_empty with corner 2^n/T."""
    P = gms.build(n, d, opts=opts, name=name or 'gmsf_%d_%d_%d' % (n, d, T), quad_lin=quad_lin) if light else split.build(n, d, 0, opts=opts, name=name)
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
    add_mempty(P, n, d, T, key3, None)
    return P


if __name__ == '__main__':
    import sdpcore
    n, d, T = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    light = len(sys.argv) <= 4 or sys.argv[4] != 'nolight'
    off = dict(best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)
    t = time.time()
    P = build_full(n, d, T, opts=off, light=light)
    print('vars', len(P.vars), 'blocks', len(P.blocks), 'lin', len(P.lin), '%.0fs' % (time.time() - t), flush=True)
    info = sdpcore.solve_csdp(P, 'agents/sdp/runs/gms0', threads=4, cpus=None)
    print(n, d, T, info['phase.value'], 'bound', info['bound'], 'pbound', info['pbound'], '%.0fs' % (time.time() - t))


def homogenize(P, n, fix_p=None, xle_p=True):
    """GMS unnormalized form: x(S) = p X(S). Constants (X=1) become variable 'p'; M_empty corner becomes 1.
    Objective: max 2^n p. If fix_p is given (Fraction), p is fixed (exact-size tests)."""
    pk = ('p',)
    P.var(pk, 0, 1)
    for b in P.blocks:
        inv = b.get('meta', {}).get('kind') == 'gms0inv'
        const = b['ent'].pop(None, {})
        for (i, j), val in const.items():
            if inv and (i, j) == (0, 0):
                b['ent'].setdefault(None, {})[(0, 0)] = Fraction(1)
            else:
                b['ent'].setdefault(pk, {})
                b['ent'][pk][(i, j)] = b['ent'][pk].get((i, j), 0) + val
    new = []
    for g0, co in P.lin:
        co = dict(co)
        if g0 != 0:
            co[pk] = co.get(pk, 0) + g0
        new.append((Fraction(0), co))
    P.lin = new
    # box: x <= p for all X-variables (X <= 1)
    if xle_p:
        for v in P.vars:
            if v != pk:
                P.lin.append((Fraction(0), {pk: Fraction(1), v: Fraction(-1)}))
    P.c = {pk: Fraction(2 ** n)}
    P.c0 = Fraction(0)
    if fix_p is not None:
        P.lin.append((Fraction(0) + fix_p, {pk: Fraction(-1)}))
        P.lin.append((-fix_p, {pk: Fraction(1)}))
    return P


def build_gms(n, d, opts=None, name=None, light=True):
    P = build_full(n, d, 1, opts=opts, name=name, light=light)
    return homogenize(P, n)
