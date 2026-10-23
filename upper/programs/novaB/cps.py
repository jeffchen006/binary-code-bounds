"""B1s: B1 rows on the two shortened subcodes C^(i,b) = {c in C: c_i = b}, averaged over i and over S_{n-1} acting on V.
For V <= F2^{n-1} (dim r, weight counts m_k), premise A(n-1,d) <= Ap (so M-Ap <= |C^(i,b)| <= Ap):
  (1/n) sum_a A_a (n-a) sum_k m_k K^{(n-1)}_k(a) / C(n-1,k)  >=  F/M,     A_0 = 1,
  F = min over M0+M1 = M, M-Ap <= M0 <= Ap of f_r(M0) + f_r(M1),  f_r(m) = rho(2^r - rho), rho = m mod 2^r.
Also the full-code B1 rows (cp.build) and the cap row."""
import sys, os, random, json
from fractions import Fraction
from math import comb
from oa import *
import oa
from cp import build as build_full, separate, wd_of, viol
import sdpcore

def f_r(m, r):
    rho = m % 2 ** r
    return rho * (2 ** r - rho)

def F(M, Ap, r):
    return min(f_r(M0, r) + f_r(M - M0, r) for M0 in range(M - Ap, Ap + 1))

def short_form(P, n, k, X0):
    """const, coeffs of (1/n) sum_a A_a (n-a) K^{(n-1)}_k(a)   (A_0 = 1)"""
    c0 = Fraction(n * kraw(n - 1, k, 0), n)
    co = {}
    for a, f in X0.items():
        K = kraw(n - 1, k, a) * (n - a)
        if K == 0: continue
        for v, c in f.items():
            co[v] = co.get(v, 0) + c * Fraction(K, n)
    return c0, co

def add_short_row(P, n, M, Ap, wd, r, SF):
    rhs = Fraction(F(M, Ap, r), M)
    if rhs == 0: return False
    g0 = -rhs; co = {}
    for k, m in wd.items():
        c0, f = SF[k]
        w = Fraction(m, comb(n - 1, k))
        g0 += w * c0
        for v, c in f.items(): co[v] = co.get(v, 0) + w * c
    co = {v: c for v, c in co.items() if c}
    mx = max([abs(c) for c in co.values()] + [abs(g0)])
    P.lin.append((g0 / mx, {v: c / mx for v, c in co.items()}))
    return True

def build(n, d, M, Ap, pool_full, pool_short):
    P, meta = build_full(n, d, M, Ap, pool_full)
    X0 = X0forms(P, n)
    SF = {k: short_form(P, n, k, X0) for k in range(n)}
    from t3 import partitions, block_wd
    for tot in range(1, n):
        for part in partitions(tot):
            r = len(part)
            if r > 12: continue
            mult = [0] * (2 ** r - 1)
            for b, a in enumerate(part): mult[(1 << b) - 1] = a
            if add_short_row(P, n, M, Ap, block_wd(part), r, SF): meta.append(('short', r, mult))
    for r, mult in pool_short:
        if add_short_row(P, n, M, Ap, wd_of(r, mult), r, SF): meta.append(('short', r, mult))
    return P, meta, SF

def short_betas(P, n, M, x, SF):
    vals = dict(zip(P.vars, x))
    out = {}
    for k in range(n):
        c0, f = SF[k]
        out[k] = M * (float(c0) + sum(float(c) * float(vals[v]) for v, c in f.items())) / comb(n - 1, k)
    return out

def separate_short(bs, n, M, Ap, rs, tries, steps=400, rng=random):
    import numpy as np
    from cp import par
    b = np.asarray([bs[k] for k in range(n)] + [1e18])
    found = {}
    for r in rs:
        Fv = F(M, Ap, r)
        if Fv == 0: continue
        T = 2 ** r - 1
        def val(mult):
            if sum(mult) > n - 1: return None
            w = par(r) @ np.asarray(mult)
            if (w == 0).any(): return None
            return b[w].sum() - Fv
        for _ in range(tries):
            mult = [0] * T
            for q in range(r): mult[(1 << q) - 1] += 1
            for _ in range(rng.randint(0, n - 1 - r)): mult[rng.randrange(T)] += 1
            cur = val(mult)
            if cur is None: continue
            for _ in range(steps):
                m2 = list(mult); i = rng.randrange(T); mv = rng.random()
                if mv < 0.4: m2[i] += 1
                elif m2[i] == 0: continue
                elif mv < 0.8: m2[i] -= 1
                else: m2[i] -= 1; m2[rng.randrange(T)] += 1
                v = val(m2)
                if v is not None and v <= cur: mult, cur = m2, v
            if cur < -1e-6:
                key = (r, tuple(mult))
                found[key] = (cur, r, mult)
    return sorted(found.values())

if __name__ == '__main__':
    n, d, M, Ap, rounds = map(int, sys.argv[1:6])
    random.seed(int(os.environ.get('SEED', '5')))
    pf, ps, keys = [], [], set()
    if len(sys.argv) > 6: pf = [tuple(e) for e in json.load(open(sys.argv[6]))['extra']]
    rmax = int(os.environ.get('RMAX', '8'))
    while True:
        ok = False
        for it in range(rounds):
            P, meta, SF = build(n, d, M, Ap, pf, ps)
            info = sdpcore.solve_csdp(P, "agents/novaB/runs/cps_%d" % n, threads=1, cpus=os.environ.get("CPU"))
            if info['bound'] is not None and info['bound'] < M: ok = True; break
            A, beta = betas(P, n, M, info['x'])
            bs = short_betas(P, n, M, info['x'], SF)
            new = 0
            for v, r, mult, wd in separate(beta, n, M, rs=range(2, rmax + 1), tries=100):
                if ('f', r, tuple(mult)) in keys: continue
                keys.add(('f', r, tuple(mult))); pf.append((r, mult)); new += 1
            for v, r, mult in separate_short(bs, n, M, Ap, range(1, rmax + 1), 100):
                if ('s', r, tuple(mult)) in keys: continue
                keys.add(('s', r, tuple(mult))); ps.append((r, mult)); new += 1
            if new == 0: break
        print(n, d, M, 'Ap', Ap, 'rounds', it + 1, 'pools', len(pf), len(ps), 'bound', info['bound'], 'EXCL' if ok else 'feas', flush=True)
        json.dump({'n': n, 'd': d, 'Ap': Ap, 'extra': pf, 'short': ps}, open('runs/cps_%d_%d_c%s.json' % (n, d, os.environ.get('CPU')), 'w'))
        if not ok: break
        M -= 1
