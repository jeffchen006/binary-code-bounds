"""Cutting-plane loop for B1 rows with GENERAL linear V.
V is given by r and a column multiset mult[c], c in F2^r \ 0 (index c-1), length l = sum mult <= n.
Generator matrix columns = the c's; message x != 0 has weight w(x) = sum_{c: x.c odd} mult[c]; V has dim r iff w(x) > 0 for all x.
Row: sum_x beta_{w(x)} >= rho(2^r - rho), rho = M mod 2^r  (beta_k = avg_{|u|=k} S_u^2).
usage: cp.py n d M Ap rounds"""
import sys, os, random, json
from fractions import Fraction
from math import comb
from oa import *
from t3 import block_rows
import sdpcore

def wd_of(r, mult):
    wd = {}
    for x in range(1, 2 ** r):
        w = sum(m for c, m in enumerate(mult, 1) if m and bin(x & c).count('1') & 1)
        if w == 0: return None
        wd[w] = wd.get(w, 0) + 1
    return wd

import numpy as np
_PAR = {}
def par(r):
    if r not in _PAR:
        X = np.arange(1, 2 ** r)
        _PAR[r] = np.array([[bin(x & c).count('1') & 1 for c in X] for x in X], dtype=np.int64)
    return _PAR[r]

def viol(beta, n, M, r, mult):
    if sum(mult) > n: return None
    w = par(r) @ np.asarray(mult, dtype=np.int64)
    if (w == 0).any(): return None
    rho = M % 2 ** r
    b = np.asarray([beta[k] for k in range(n + 1)])
    v = b[w].sum() - rho * (2 ** r - rho)
    return v, None

def separate(beta, n, M, rs=None, tries=None, steps=400, rng=random):
    rs = rs or range(2, int(os.environ.get("RMAX", "8")) + 1)
    tries = tries or int(os.environ.get("TRIES", "120"))
    found = {}
    for r in rs:
        rho = M % 2 ** r
        if rho == 0: continue
        T = 2 ** r - 1
        for _ in range(tries):
            # random start: r independent unit columns + random extras
            mult = [0] * T
            for b in range(r): mult[(1 << b) - 1] += 1
            for _ in range(rng.randint(0, n - r)):
                mult[rng.randrange(T)] += 1
            if sum(mult) > n: continue
            cur = viol(beta, n, M, r, mult)
            if cur is None: continue
            for _ in range(steps):
                m2 = list(mult)
                mv = rng.random()
                i = rng.randrange(T)
                if mv < 0.4: m2[i] += 1
                elif mv < 0.8:
                    if m2[i] == 0: continue
                    m2[i] -= 1
                else:
                    j = rng.randrange(T)
                    if m2[i] == 0: continue
                    m2[i] -= 1; m2[j] += 1
                v = viol(beta, n, M, r, m2)
                if v is None: continue
                if v[0] <= cur[0]:
                    mult, cur = m2, v
            if cur[0] < -1e-6:
                cur = (cur[0], wd_of(r, mult))
                key = (r, tuple(sorted(cur[1].items())))
                if key not in found or found[key][0] > cur[0]:
                    found[key] = (cur[0], r, list(mult), cur[1])
    return sorted(found.values())

def build(n, d, M, Ap, extra):
    opts = dict(lin_basic=True, best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)
    P = split.build(n, d, 0, opts=opts, name='cp_%d_%d_%d' % (n, d, M))
    meta = [None] * len(P.lin)
    seen = set()
    for r, wd in block_rows(n, M, rmin=1):
        seen.add((r, tuple(sorted(wd.items()))))
        # block V as column multiset: block b -> column e_b with multiplicity a_b; recover blocks from wd not needed: store blocks
    from t3 import partitions, block_wd
    for tot in range(1, n + 1):
        for part in partitions(tot):
            r = len(part)
            if M % 2 ** r == 0: continue
            wd = block_wd(part); key = (r, tuple(sorted(wd.items())))
            if key in [m[0] for m in meta if m] : pass
            mult = [0] * (2 ** r - 1)
            for b, a in enumerate(part): mult[(1 << b) - 1] = a
            if add_oa_row(P, n, M, wd, r): meta.append(((r, tuple(sorted(wd.items()))), r, mult))
    for (r, mult) in extra:
        wd = wd_of(r, mult)
        if add_oa_row(P, n, M, wd, r): meta.append(((r, tuple(sorted(wd.items()))), r, mult))
    if Ap:
        c0, f = shell_form(P, n, 1)
        eps = Fraction(os.environ.get('CAPEPS', '0'))
        P.lin.append((Fraction(n * (2 * Ap - M) ** 2, M) + eps - c0, {v: -c for v, c in f.items()}))
        meta.append(('cap', Ap, eps))
    return P, meta

if __name__ == '__main__':
    n, d, M, Ap, rounds = map(int, sys.argv[1:6])
    random.seed(12345)
    extra = []; keys = set()
    work = "agents/novaB/runs/cp"
    for it in range(rounds):
        P, meta = build(n, d, M, Ap, extra)
        info = sdpcore.solve_csdp(P, work, threads=1, cpus=os.environ.get("CPU"))
        A, beta = betas(P, n, M, info['x'])
        print('round', it, 'rows', len(P.lin), 'bound', info['bound'], 'EXCL' if info['bound'] < M else 'feas', flush=True)
        if info['bound'] < M: break
        cuts = separate(beta, n, M)
        new = 0
        for v, r, mult, wd in cuts:
            key = (r, tuple(sorted(wd.items())))
            if key in keys: continue
            keys.add(key); extra.append((r, mult)); new += 1
        print('  cuts found', len(cuts), 'new', new, 'best', cuts[0][:3] if cuts else None, flush=True)
        if new == 0: break
    json.dump({'n': n, 'd': d, 'M': M, 'Ap': Ap, 'extra': extra}, open('runs/cp_%d_%d_%d.json' % (n, d, M), 'w'))
