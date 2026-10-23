"""Ablation: Schrijver SDP + OLD rows only. usage: abl.py n d M Ap case w
case: odd (r=1 rows), 2mod4 (all r=2 rows + BBMOS Thm8/Kim-Toan split for w), 0mod4 (split rows for w; w=0: none)."""
import sys, os
from fractions import Fraction
from math import comb
from oa import *
import sdpcore
n, d, M, Ap = map(int, sys.argv[1:5]); case = sys.argv[5]; w = int(sys.argv[6])
opts = dict(lin_basic=True, best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)
P = split.build(n, d, 0, opts=opts, name='abl_%d_%d_%d_%s_%d' % (n, d, M, case, w))
X0 = X0forms(P, n)
def shell_row(k, rhs):
    c0, f = shell_form(P, n, k, X0)
    g0 = c0 - rhs
    if not f: return
    mx = max([abs(c) for c in f.values()] + [abs(g0)])
    P.lin.append((g0 / mx, {v: c / mx for v, c in f.items()}))
if case == 'odd':
    assert M % 2 == 1
    for k in range(1, n + 1): shell_row(k, Fraction(comb(n, k), M))
elif case == '2mod4':
    assert M % 4 == 2
    for k in range(1, n + 1): shell_row(k, Fraction(2 * (comb(n, k) + kraw(n, k, w)), M))
    for p in range(n + 1):
        for q in range(n + 1 - p):
            for s in range(n + 1 - p - q):
                a, b, c = p + s, q + s, p + q
                if min(a, b, c) == 0: continue
                wd = {}
                for x in (a, b, c): wd[x] = wd.get(x, 0) + 1
                add_oa_row(P, n, M, wd, 2)
elif case == '0mod4':
    assert M % 4 == 0
    for k in range(1, n + 1):
        rhs = Fraction(2 * (comb(n, k) - kraw(n, k, w)), M)
        if rhs: shell_row(k, rhs)
if Ap:
    c0, f = shell_form(P, n, 1, X0)
    P.lin.append((Fraction(n * (2 * Ap - M) ** 2, M) - c0, {v: -c for v, c in f.items()}))
info = sdpcore.solve_csdp(P, 'agents/novaB/abl/w%s' % os.environ.get('CPU'), threads=1, cpus=os.environ.get('CPU'))
print(n, d, M, Ap, case, w, info['phase.value'][:24].replace(' ', '_'), info['bound'], flush=True)
