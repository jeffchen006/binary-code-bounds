#!/usr/bin/env python3
"""B1 separation at a finished gms.py solution (|C| = M fixed, y1 = M).  A_i = C(n,i) y2(i)/M,
beta_k = (M K_k(0) + sum_i C(n,i) K_k(i) y2(i)) / C(n,k) = average of S_a^2 over |a| = k.
B1 row (V, dim r, weight counts m_k): sum_k m_k beta_k >= rho (2^r - rho), rho = M mod 2^r.
usage: sep_gms.py n d M ub sdpa_out (--parity | --b1 RMAX | --mod4 w)"""
import sys, os, json, random, time, argparse
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, 'agents/novaB')
import gms
import cp
ap = argparse.ArgumentParser()
ap.add_argument('n', type=int); ap.add_argument('d', type=int); ap.add_argument('M', type=int); ap.add_argument('ub', type=int)
ap.add_argument('out'); ap.add_argument('--parity', action='store_true'); ap.add_argument('--b1', type=int); ap.add_argument('--cw', action='store_true')
a = ap.parse_args()
n, d, M = a.n, a.d, a.M
t = time.time()
P = gms.build(n, d, Fraction(a.ub), verbose=False)
if a.parity:
    gms.add_parity(P, n, d, M)
else:
    gms.add_b1(P, n, d, M, a.b1)
if a.cw:
    gms.add_cw(P, n, d, M)
print('built %.0fs' % (time.time() - t), flush=True)
info = P.read_result(a.out)
val = {v: float(info['x'][i]) for i, v in enumerate(P.vars)}
print('phase', info['phase.value'], 'objP', info['objValPrimal'], 'objD', info['objValDual'], 'y1 =', val[gms.pair_key(n, 0)])
D = [j for j in range(d, n + 1) if j % 2 == 0]
A = {0: 1.0}
for i in D:
    k = gms.pair_key(n, i)
    A[i] = gms.C(n, i) * val.get(k, 0.0) / M
beta = {k: M * sum(A[i] * gms.Kr(n, k, i) for i in A) / gms.C(n, k) for k in range(n + 1)}
print('A:', {i: round(v, 5) for i, v in A.items()})
print('beta:', {k: round(v, 4) for k, v in beta.items()})
# block rows, exhaustive
def partitions(total, maxpart=None):
    maxpart = total if maxpart is None else maxpart
    if total == 0:
        yield (); return
    for p in range(min(total, maxpart), 0, -1):
        for rest in partitions(total - p, p):
            yield (p,) + rest
res = []
for tot in range(1, n + 1):
    for pt in partitions(tot):
        r = len(pt)
        rho = M % 2 ** r
        if rho == 0:
            continue
        wd = {0: 1}
        for x in pt:
            nw = dict(wd)
            for w_, c in wd.items():
                nw[w_ + x] = nw.get(w_ + x, 0) + c
            wd = nw
        wd.pop(0)
        lhs = sum(m * beta[k] for k, m in wd.items())
        rhs = rho * (2 ** r - rho)
        res.append((lhs / rhs, r, pt))
res.sort()
print('block rows: %d; smallest lhs/rhs (< 1 = violated):' % len(res))
for z in res[:15]:
    print('  %.5f r=%d %s' % z)
random.seed(1)
cuts = cp.separate(beta, n, M, rs=range(2, 11), tries=80, steps=500)
print('general V with violation (lhs - rhs < 0):', len(cuts))
for c in cuts[:10]:
    print('  %.4f r=%d mult=%s' % (c[0], c[1], c[2]))
json.dump({'A': A, 'beta': beta, 'block': [(z[0], z[1], list(z[2])) for z in res[:300]], 'cuts': [(c[0], c[1], c[2]) for c in cuts[:300]]},
          open(a.out + '.sep.json', 'w'))
print('done %.0fs' % (time.time() - t))
