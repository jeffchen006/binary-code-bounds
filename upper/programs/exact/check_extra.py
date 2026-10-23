#!/usr/bin/env python3
"""Evaluate the extra valid rows of gms.add_extra (cw, best18, mel, kt) at a finished gms.py solution.
usage: check_extra.py n d M ub out (--parity | --b1 R)"""
import sys, os, argparse
from fractions import Fraction
sys.path.append('agents/sdp')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gms
ap = argparse.ArgumentParser()
ap.add_argument('n', type=int); ap.add_argument('d', type=int); ap.add_argument('M', type=int); ap.add_argument('ub', type=int)
ap.add_argument('out'); ap.add_argument('--parity', action='store_true'); ap.add_argument('--b1', type=int); ap.add_argument('--cw', action='store_true')
a = ap.parse_args()
P = gms.build(a.n, a.d, Fraction(a.ub), verbose=False)
if a.parity:
    gms.add_parity(P, a.n, a.d, a.M)
else:
    gms.add_b1(P, a.n, a.d, a.M, a.b1)
if a.cw:
    gms.add_cw(P, a.n, a.d, a.M)
info = P.read_result(a.out)
val = {v: float(info['x'][i]) for i, v in enumerate(P.vars)}
n0 = len(P.lin)
for opt in ['cw', 'best18', 'mel', 'kt']:
    k0 = len(P.lin)
    try:
        gms.add_extra(P, a.n, a.d, {opt: True}, verbose=False)
    except Exception as e:
        print(opt, 'failed', repr(e)); continue
    worst = None
    for g0, co in P.lin[k0:]:
        v = float(g0) + sum(float(c) * val.get(x, 0.0) for x, c in co.items())
        sc = max([abs(float(c)) for c in co.values()] + [abs(float(g0)), 1e-30])
        if worst is None or v / sc < worst[0]:
            worst = (v / sc, v)
    print(opt, 'rows', len(P.lin) - k0, 'min normalized slack', worst)
