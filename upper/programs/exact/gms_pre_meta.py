#!/usr/bin/env python3
"""Gijswijt-Mittelmann-Schrijver quadruple SDP bound A4(n,d) (IEEE IT 58 (2012), arXiv:1005.4959),
with optional extra valid constraints.  d even; all distances restricted to D = {0} u [d,n]_even.

Variables: y(S) = 2^n x(S) for each G-orbit of codes S with 1 <= |S| <= 4 (x(empty) = 1).
For a real code C, x(S) = Pr[g(S) subset of C] over g in G, so y({v}) = |C| = objective.

PSD constraints:
 * M_empty(x): block diagonalized as in GMS Sec. 5 (ordered pairs, Sym^n(A), psi_{m,k,l}), with the
   Terwilliger blocks Gamma_{alpha,k} replaced by Schrijver's integer beta matrices (diagonally congruent),
   rows restricted to distance in D and i' even (unordered pairs); block (n,0,0) replaced by the
   G-invariant block including the empty set.
 * M_S(x), S = {0,u}, |u| = m in D: Sym^m (x) Sym^(n-m) blocks (GMS Sec. 4).
 * y >= 0 (box), y <= ub.
"""
import sys, os, itertools, math, time, argparse
from fractions import Fraction
from functools import lru_cache
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'sdp'))
pass
import mpmath
from sdpcore import Problem

comb = math.comb
S2RED = os.environ.get('GMS_S2', '1') == '1'
S2CHECK = os.environ.get('GMS_S2CHECK', '0') == '1'


def C(a, b):
    return comb(a, b) if 0 <= b <= a else 0


@lru_cache(None)
def Kr(N, k, t):
    """Krawtchouk K_k^N(t) = sum_i (-1)^i C(t,i) C(N-t,k-i)."""
    return sum((-1) ** i * C(t, i) * C(N - t, k - i) for i in range(k + 1))


@lru_cache(None)
def beta(m, i, j, k, t):
    """Schrijver 2005 block coefficient (unnormalized)."""
    s = 0
    for u in range(0, m + 1):
        b = C(u, t) * C(m - 2 * k, u - k) * C(m - k - u, i - u) * C(m - k - u, j - u)
        if b:
            s += (-1) ** (u - t) * b
    return s


P = [(0, 0), (0, 1), (1, 0), (1, 1)]
T = [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]


def comps(n, k):
    if k == 1:
        yield (n,)
        return
    for i in range(n + 1):
        for r in comps(n - i, k - 1):
            yield (i,) + r


_canon_cache = {}


def canon(cols):
    """cols: tuple of (column_bits_tuple, count) for r words (with possible duplicates).
    Returns orbit key of the SET of distinct words under coordinate permutations + translations."""
    if cols in _canon_cache:
        return _canon_cache[cols]
    r = len(cols[0][0])
    act = [(c, k) for c, k in cols if k > 0]
    # distinct words
    reps = []
    for a in range(r):
        dup = False
        for b in reps:
            if all(c[a] == c[b] for c, k in act):
                dup = True
                break
        if not dup:
            reps.append(a)
    best = None
    for perm in itertools.permutations(reps):
        dct = {}
        for c, k in act:
            col = tuple(c[p] ^ c[perm[0]] for p in perm)
            dct[col] = dct.get(col, 0) + k
        key = (len(reps), tuple(sorted(dct.items())))
        if best is None or key < best:
            best = key
    _canon_cache[cols] = best
    return best


def dists_ok(cols, D):
    r = len(cols[0][0])
    for a in range(r):
        for b in range(a + 1, r):
            dd = sum(k for c, k in cols if c[a] != c[b])
            if dd not in D:
                return False
    return True


def build(n, d, ub, name=None, verbose=True, skip_MS=False):
    t0 = time.time()
    D = set([0] + [j for j in range(d, n + 1) if j % 2 == 0])
    Dp = sorted(D - {0})
    Pr = Problem(name or 'gms_%d_%d' % (n, d))
    one = canon((((0,), n),))
    def V(key):
        return Pr.var(key, 0, ub)
    V(one)
    Pr.c = {one: Fraction(1)}
    Pr.c0 = Fraction(0)

    # ---------------- M_empty ----------------
    lam_valid = []
    for lam in comps(n, 8):
        cols = tuple(((t1, 0, t2 ^ t3, t3), l) for (t1, t2, t3), l in zip(T, lam))
        if not dists_ok(cols, D):
            continue
        lam_valid.append((lam, canon(cols)))
    if verbose:
        print('valid lambda', len(lam_valid), 'vars so far', len(Pr.vars), '%.1fs' % (time.time() - t0), flush=True)
    blocks = {}   # (m,k,l) -> {var: {(rowkey,colkey): int}}
    Zacc = {}     # (a,b) -> {var: Fraction}
    for lam, key in lam_valid:
        V(key)
        lt = dict(zip(T, lam))
        lp = {c: lt[(c[0], c[1], 0)] + lt[(c[0], c[1], 1)] for c in P}
        l1 = {c: lt[(c[0], c[1], 1)] for c in P}
        a_ = lp[(1, 0)] + lp[(1, 1)]
        b_ = lp[(0, 1)] + lp[(1, 1)]
        mult = math.factorial(n)
        for l in lam:
            mult //= math.factorial(l)
        Zacc.setdefault((a_, b_), {})
        Zacc[(a_, b_)][key] = Zacc[(a_, b_)].get(key, 0) + Fraction(mult, C(n, a_) * C(n, b_))
        ranges = [range(lp[c] + 1) for c in P]
        for al in itertools.product(*ranges):
            alpha = dict(zip(P, al))
            bet = {c: lp[c] - alpha[c] for c in P}
            th = 1
            for c in P:
                th *= Kr(lp[c], l1[c], bet[c])
                if th == 0:
                    break
            if th == 0:
                continue
            m = sum(al)
            i, j, t = alpha[(1, 0)] + alpha[(1, 1)], alpha[(0, 1)] + alpha[(1, 1)], alpha[(1, 1)]
            ip, jp, tp = bet[(1, 0)] + bet[(1, 1)], bet[(0, 1)] + bet[(1, 1)], bet[(1, 1)]
            if ip % 2 or jp % 2:
                continue
            for k in range(0, min(i, j, m - i, m - j) + 1):
                if m == n and k == 0:
                    continue  # block (n,0,0) replaced by Z block
                bk = beta(m, i, j, k, t)
                if bk == 0:
                    continue
                for l in range(0, min(ip, jp, n - m - ip, n - m - jp) + 1):
                    bl = beta(n - m, ip, jp, l, tp)
                    if bl == 0:
                        continue
                    blk = blocks.setdefault((m, k, l), {})
                    ent = blk.setdefault(key, {})
                    rc = ((i, ip), (j, jp))
                    ent[rc] = ent.get(rc, 0) + th * bk * bl
    if verbose:
        print('M_empty accumulated: blocks', len(blocks), '%.1fs' % (time.time() - t0), flush=True)
    nblk = 0
    for (m, k, l), blk in sorted(blocks.items()):
        rows = sorted(set(rc[0] for ent in blk.values() for rc in ent) | set(rc[1] for ent in blk.values() for rc in ent))
        pos = {r: p for p, r in enumerate(rows)}
        ent2 = {}
        for key, ent in blk.items():
            mat = {}
            for (r1, r2), val in ent.items():
                if val == 0:
                    continue
                a, b = pos[r1], pos[r2]
                if a <= b:
                    mat[(a, b)] = mat.get((a, b), 0) + val
                # symmetry check
                if ent.get((r2, r1), 0) != val:
                    raise RuntimeError('asymmetric block %s' % ((m, k, l),))
            mat = {ab: Fraction(v) for ab, v in mat.items() if v != 0}
            if mat:
                ent2[key] = mat
        if not ent2:
            continue
        scale = [1 / mpmath.sqrt(C(m - 2 * k, r[0] - k) * C(n - m - 2 * l, r[1] - l)) for r in rows]
        Pr.add_block(len(rows), ent2, scale)
        nblk += 1
    # Z block: rows empty, O_a (a in D)
    Drows = sorted(D)
    zent = {None: {(0, 0): Fraction(1)}}
    for p, a in enumerate(Drows):
        key = one if a == 0 else canon((((0, 0), n - a), ((0, 1), a)))
        V(key)
        zent.setdefault(key, {})[(0, p + 1)] = Fraction(1)
    for (a, b), dct in Zacc.items():
        pa, pb = Drows.index(a), Drows.index(b)
        if pa > pb:
            continue
        for key, val in dct.items():
            zent.setdefault(key, {})
            zent[key][(pa + 1, pb + 1)] = zent[key].get((pa + 1, pb + 1), 0) + val
    Pr.add_block(len(Drows) + 1, zent, None)
    nblk += 1
    if verbose:
        print('M_empty blocks', nblk, 'vars', len(Pr.vars), '%.1fs' % (time.time() - t0), flush=True)

    # ---------------- M_S, S = {0,u}, |u| = m ----------------
    if not skip_MS:
        for m in Dp:
            n2 = n - m
            blk = {}
            Lm = [a for a in comps(m, 4)]
            Ln = [b for b in comps(n2, 4)]
            # row validity
            def rowok(i, ip):
                return (i + ip) in D and (i + ip) != 0 and (m - i + ip) in D
            for al in Lm:
                alpha = dict(zip(P, al))
                i, j, t = alpha[(1, 0)] + alpha[(1, 1)], alpha[(0, 1)] + alpha[(1, 1)], alpha[(1, 1)]
                for bl_ in Ln:
                    bet = dict(zip(P, bl_))
                    ip, jp, tp = bet[(1, 0)] + bet[(1, 1)], bet[(0, 1)] + bet[(1, 1)], bet[(1, 1)]
                    if not (rowok(i, ip) and rowok(j, jp)):
                        continue
                    dvw = alpha[(1, 0)] + alpha[(0, 1)] + bet[(1, 0)] + bet[(0, 1)]
                    if dvw not in D:
                        continue
                    cols = tuple([((0, 1, c[0], c[1]), alpha[c]) for c in P] + [((0, 0, c[0], c[1]), bet[c]) for c in P])
                    if not dists_ok(cols, D):
                        continue
                    key = canon(cols)
                    V(key)
                    for k in range(0, min(i, j, m - i, m - j) + 1):
                        bk = beta(m, i, j, k, t)
                        if bk == 0:
                            continue
                        for l in range(0, min(ip, jp, n2 - ip, n2 - jp) + 1):
                            bb = beta(n2, ip, jp, l, tp)
                            if bb == 0:
                                continue
                            e = blk.setdefault((k, l), {}).setdefault(key, {})
                            rc = ((i, ip), (j, jp))
                            e[rc] = e.get(rc, 0) + bk * bb
            for (k, l), bd in sorted(blk.items()):
                rows = sorted(set(rc[0] for ent in bd.values() for rc in ent) | set(rc[1] for ent in bd.values() for rc in ent))
                if not S2RED or (k == 0 and l == 0):
                    pos = {r: p for p, r in enumerate(rows)}
                    ent2 = {}
                    for key, ent in bd.items():
                        mat = {}
                        for (r1, r2), val in ent.items():
                            if ent.get((r2, r1), 0) != val:
                                raise RuntimeError('asymmetric M_S block')
                            a, b = pos[r1], pos[r2]
                            if a <= b and val:
                                mat[(a, b)] = Fraction(val)
                        if mat:
                            ent2[key] = mat
                    if ent2:
                        scale = [1 / mpmath.sqrt(C(m - 2 * k, r[0] - k) * C(n2 - 2 * l, r[1] - l)) for r in rows]
                        Pr.add_block(len(rows), ent2, scale)
                        nblk += 1
                    continue
                # S2 reduction: sigma: v -> v+u acts on rows as (i,i') -> (m-i,i') with sign eps(i) = sign beta(m,i,m-i,k,0)
                rows = [(i, ip) for i in range(k, m - k + 1) for ip in range(l, n2 - l + 1) if rowok(i, ip)]
                rset = set(rows)
                cols = {+1: [], -1: []}   # list of [(row, coef), ...]
                done = set()
                for r in rows:
                    if r in done:
                        continue
                    i, ip = r
                    rs = (m - i, ip)
                    if rs not in rset:
                        raise RuntimeError('row set not sigma-invariant')
                    eb = beta(m, i, m - i, k, 0)
                    if eb == 0:
                        raise RuntimeError('zero sigma entry')
                    eps = 1 if eb > 0 else -1
                    if rs == r:
                        cols[eps].append([(r, 1)])
                    else:
                        cols[+1].append([(r, 1), (rs, eps)])
                        cols[-1].append([(r, 1), (rs, -eps)])
                        done.add(rs)
                    done.add(r)
                for sgn in (+1, -1):
                    cl = cols[sgn]
                    if not cl:
                        continue
                    ent2 = {}
                    for key, ent in bd.items():
                        mat = {}
                        for a in range(len(cl)):
                            for b in range(a, len(cl)):
                                v = 0
                                for (ra, ca) in cl[a]:
                                    for (rb, cb) in cl[b]:
                                        v += ca * cb * ent.get((ra, rb), 0)
                                if v:
                                    mat[(a, b)] = Fraction(v)
                        if mat:
                            ent2[key] = mat
                    if ent2:
                        scale = [1 / mpmath.sqrt(len(c) * C(m - 2 * k, c[0][0][0] - k) * C(n2 - 2 * l, c[0][0][1] - l)) for c in cl]
                        Pr.add_block(len(cl), ent2, scale)
                        nblk += 1
                # check the cross part vanishes
                if S2CHECK and cols[+1] and cols[-1]:
                    for key, ent in bd.items():
                        for ca_ in cols[+1]:
                            for cb_ in cols[-1]:
                                v = sum(ca * cb * ent.get((ra, rb), 0) for (ra, ca) in ca_ for (rb, cb) in cb_)
                                if v:
                                    raise RuntimeError('S2 cross term nonzero m=%d k=%d l=%d' % (m, k, l))
            if verbose:
                print('M_S m=%d done, blocks %d vars %d %.1fs' % (m, nblk, len(Pr.vars), time.time() - t0), flush=True)
    return Pr


# ---------------- extra valid linear constraints ----------------
def pair_key(n, a):
    return canon((((0,), n),)) if a == 0 else canon((((0, 0), n - a), ((0, 1), a)))


def triple_key(n, i, j, t):
    """words 0, u (|u|=i), w (|w|=j), |u&w| = t"""
    return canon((((0, 1, 1), t), ((0, 1, 0), i - t), ((0, 0, 1), j - t), ((0, 0, 0), n - i - j + t)))


def add_extra(Pr, n, d, opts, verbose=True):
    import bounds
    D = set([0] + [j for j in range(d, n + 1) if j % 2 == 0])
    Dp = sorted(D - {0})
    one = pair_key(n, 0)
    cnt = 0
    def linA(g0, terms):
        """g0 + sum coef*A_a >= 0, A_a = C(n,a) y2(a)/y1 (A_0 = 1); homogenized by y1 > 0."""
        nonlocal cnt
        co = {one: Fraction(g0)}
        for coef, a in terms:
            if a == 0:
                co[one] = co.get(one, 0) + Fraction(coef)
            elif a in D:
                k = pair_key(n, a)
                if k in Pr.vidx:
                    co[k] = co.get(k, 0) + Fraction(coef) * C(n, a)
        co = {k: v for k, v in co.items() if v != 0}
        if co:
            Pr.add_lin(0, co); cnt += 1
    h = d // 2
    if opts.get('cw'):
        for a in Dp:
            linA(bounds.cw_upper(n, d, a), [(-1, a)])
    if opts.get('best18'):
        for kk in range(1, n):
            try:
                ub = bounds.A_upper(n - kk, d)
            except Exception:
                continue
            linA(C(n, kk) * ub, [(-C(n - a, kk), a) for a in range(n + 1)])
    if opts.get('mel'):
        q = (2 * n) // d
        linA(q, [(-1, n - h)] + [(-q, n - i) for i in range(h)])
        for i in range(1, h):
            S = bounds.cw_upper(n, d, h + i)
            Tt = bounds.cw_upper(n - h + i, d, h + i)
            linA(S, [(-1, n - h - i), (-(S - Tt), n - h + i)] + [(-S, n - h + j) for j in range(i + 1, h + 1)])
        Pp = bounds.cw_upper(n - 1, d, h + 1)
        Q = bounds.cw_upper(n - h, d, h + 1)
        R = bounds.cw_upper(n - h + 2, d, h + 2)
        terms = [(-(h + 2), n - h - 2), (-h * (Pp - Q), n - h), (-(n * Pp - (h + 2) * R), n - h + 2)]
        terms += [(-n * Pp, a) for a in range(n - h + 3, n + 1)]
        linA(n * Pp, terms)
    if opts.get('kt'):
        for i in Dp:
            for j in Dp:
                for t in range(max(0, i + j - n), min(i, j) + 1):
                    c = i + j - 2 * t
                    if c not in D or c == 0:
                        continue
                    Tv = bounds.T_upper(t, i, j - t, n - i, d)
                    tot = C(i, t) * C(n - i, j - t)
                    if Tv >= tot:
                        continue
                    k3 = triple_key(n, i, j, t); k2 = pair_key(n, i)
                    if k3 not in Pr.vidx:
                        continue
                    Pr.add_lin(0, {k2: Fraction(Tv, tot), k3: Fraction(-1)}); cnt += 1
    if verbose:
        print('extra linear constraints:', cnt, opts, flush=True)


def _add_norm(Pr, g0, co):
    mx = max([abs(Fraction(v)) for v in co.values()] + [abs(Fraction(g0))])
    Pr.add_lin(Fraction(g0) / mx, {k: Fraction(v) / mx for k, v in co.items()})


def add_mod4min(Pr, n, d, M):
    """M = 2 mod 4, all cases w at once: row k uses min over even w of 2(C(n,k)+K_k(w)) (valid for every w)."""
    assert M % 4 == 2
    D = [j for j in range(d, n + 1) if j % 2 == 0]
    one = pair_key(n, 0)
    Pr.add_lin(-M, {one: Fraction(1)})
    Pr.add_lin(M, {one: Fraction(-1)})
    s = Pr.var('s_parity', -1000, 1000)
    for k in range(1, n + 1):
        rhs = min(2 * (C(n, k) + Kr(n, k, w)) for w in range(0, n + 1, 2))
        co = {s: Fraction(-rhs)} if rhs else {}
        for i in D:
            kk = pair_key(n, i)
            if kk in Pr.vidx:
                co[kk] = co.get(kk, 0) + Fraction(C(n, i) * Kr(n, k, i))
        _add_norm(Pr, M * Kr(n, k, 0), co)
    Pr.c = {s: Fraction(1)}
    Pr.c0 = Fraction(0)


def add_mod4(Pr, n, d, M, w):
    """|C| = M = 2 mod 4, u = sum of codewords (mod 2) has weight w (invariant under translations/perms).
    S_a = sum_c (-1)^{a.c} = 2 mod 4 when a.u = 0, so sum_{|a|=k} S_a^2 >= 4 #{a: |a|=k, a.u=0} = 2(C(n,k)+K_k(w)).
    Rows: M*K_k(0) + sum_i C(n,i) K_k(i) y2(i) >= s * 2(C(n,k)+K_k(w)); objective max s; s < 1 excludes (M, w)."""
    assert M % 4 == 2
    D = [j for j in range(d, n + 1) if j % 2 == 0]
    one = pair_key(n, 0)
    Pr.add_lin(-M, {one: Fraction(1)})
    Pr.add_lin(M, {one: Fraction(-1)})
    s = Pr.var('s_parity', -1000, 1000)
    for k in range(1, n + 1):
        rhs = 2 * (C(n, k) + Kr(n, k, w))
        co = {s: Fraction(-rhs)} if rhs else {}
        for i in D:
            kk = pair_key(n, i)
            if kk in Pr.vidx:
                co[kk] = co.get(kk, 0) + Fraction(C(n, i) * Kr(n, k, i))
        _add_norm(Pr, M * Kr(n, k, 0), co)
    Pr.c = {s: Fraction(1)}
    Pr.c0 = Fraction(0)


def add_parity(Pr, n, d, M):
    """|C| = M fixed (odd).  For odd M and every a with |a| = k: |sum_c (-1)^{a.c}| >= 1, so
    M * sum_i A_i K_k(i) = sum_{|a|=k} (sum_c (-1)^{a.c})^2 >= C(n,k).
    With y1 = M:  M*K_k(0) + sum_i C(n,i) K_k(i) y2(i) >= s*C(n,k); objective max s."""
    assert M % 2 == 1
    D = [j for j in range(d, n + 1) if j % 2 == 0]
    one = pair_key(n, 0)
    Pr.add_lin(-M, {one: Fraction(1)})
    Pr.add_lin(M, {one: Fraction(-1)})
    s = Pr.var('s_parity', -1000, 1000)
    for k in range(1, n + 1):
        co = {s: Fraction(-C(n, k))}
        for i in D:
            kk = pair_key(n, i)
            if kk in Pr.vidx:
                co[kk] = co.get(kk, 0) + Fraction(C(n, i) * Kr(n, k, i))
        _add_norm(Pr, M * Kr(n, k, 0), co)
    Pr.c = {s: Fraction(1)}
    Pr.c0 = Fraction(0)


def add_b1(Pr, n, d, M, rmax=64, extra_wd=()):
    """|C| = M fixed.  B1 coset-count rows (novaB METHOD_B1): for linear V (dim r, weight counts m_k of nonzero words),
    sum_{v in V, v != 0} S_v^2 >= rho(2^r - rho), rho = M mod 2^r.  S_n-average:
    sum_k m_k/C(n,k) * R_k >= s * rho(2^r - rho),  R_k = M K_k(0) + sum_i C(n,i) K_k(i) y2(i) = sum_{|u|=k} S_u^2.
    V = block subspaces (all partitions, r <= rmax) + extra weight enumerators. r = 1 rows = parity rows. Objective max s."""
    sys.path.insert(0, 'agents/novaB')
    from t3 import block_rows
    D = [j for j in range(d, n + 1) if j % 2 == 0]
    one = pair_key(n, 0)
    Pr.add_lin(-M, {one: Fraction(1)})
    Pr.add_lin(M, {one: Fraction(-1)})
    s = Pr.var('s_parity', -1000, 1000)
    Rk = {}
    for k in range(1, n + 1):
        co = {}
        for i in D:
            kk = pair_key(n, i)
            if kk in Pr.vidx:
                co[kk] = co.get(kk, 0) + Fraction(C(n, i) * Kr(n, k, i))
        Rk[k] = (Fraction(M * Kr(n, k, 0)), co)
    rows = [(r, wd) for r, wd in block_rows(n, M, rmin=1) if r <= rmax] + list(extra_wd)
    seen = set()
    cnt = 0
    for r, wd in rows:
        key = (r, tuple(sorted(wd.items())))
        if key in seen:
            continue
        seen.add(key)
        rho = M % 2 ** r
        rhs = rho * (2 ** r - rho)
        if rhs == 0:
            continue
        g0 = Fraction(0); co = {s: Fraction(-rhs)}
        for k, m in wd.items():
            wgt = Fraction(m, C(n, k))
            g0 += wgt * Rk[k][0]
            for v, c in Rk[k][1].items():
                co[v] = co.get(v, 0) + wgt * c
        co = {v: c for v, c in co.items() if c}
        _add_norm(Pr, g0, co)
        cnt += 1
    Pr.c = {s: Fraction(1)}
    Pr.c0 = Fraction(0)
    print('B1 rows:', cnt, flush=True)
    return cnt


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('n', type=int); ap.add_argument('d', type=int)
    ap.add_argument('--ub', type=float, default=None)
    ap.add_argument('--work', default='/tmp/work')
    ap.add_argument('--skipMS', action='store_true')
    ap.add_argument('--threads', type=int, default=10)
    ap.add_argument('--param', default=None, help='SDPA parameter file (e.g. param_dd_loose.sdpa)')
    ap.add_argument('--noMS', action='store_true', help='alias of --skipMS')
    ap.add_argument('--extra', default='', help='comma list: cw,best18,mel,kt')
    ap.add_argument('--mod4min', type=int, default=None, help='M = 2 mod 4, min over w per row (one run for all w)')
    ap.add_argument('--mod4', type=int, nargs=2, default=None, metavar=('M', 'w'), help='M = 2 mod 4 exclusion for sum-weight w')
    ap.add_argument('--b1', type=int, default=None, help='fix |C| = M, B1 block rows (incl. parity rows), maximize s; s<1 excludes M')
    ap.add_argument('--b1rmax', type=int, default=64)
    ap.add_argument('--parity', type=int, default=None, help='fix |C| = M (odd) and maximize s in sum_i A_i K_k(i) >= s*C(n,k)/M; s<1 excludes M')
    a = ap.parse_args()
    ub = Fraction(a.ub) if a.ub else Fraction(2 ** a.n)
    Pr = build(a.n, a.d, ub, skip_MS=(a.skipMS or a.noMS))
    if a.extra:
        add_extra(Pr, a.n, a.d, {k: True for k in a.extra.split(',')})
    if a.parity:
        add_parity(Pr, a.n, a.d, a.parity)
    if a.b1:
        add_b1(Pr, a.n, a.d, a.b1, a.b1rmax)
    if a.mod4min:
        add_mod4min(Pr, a.n, a.d, a.mod4min)
    if a.mod4:
        add_mod4(Pr, a.n, a.d, a.mod4[0], a.mod4[1])
    print('vars', len(Pr.vars), 'blocks', len(Pr.blocks), 'maxblk', max(b['size'] for b in Pr.blocks), flush=True)
    info = Pr.solve(a.work, param=(os.path.abspath(a.param) if a.param else None), threads=a.threads)
    print(a.n, a.d, info['phase.value'], 'bound', info['bound'], 'pbound', info['pbound'], flush=True)


