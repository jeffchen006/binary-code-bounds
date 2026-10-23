#!/usr/bin/env python3
"""Schrijver 3-point SDP (agents/sdp/schrijver.py, read-only import) for an even (n,M,d) code with
M >= Mfix: max / min of S = sum_j j*A_j.  A coordinate with max part m gives sum_i w_i(M-w_i) and, if
every coordinate has max part <= m, then S >= 2*n*m*(M-m)/M.  If maxS < that, the code must have a
coordinate with max part > m."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'sdp'))
pass
from fractions import Fraction
import schrijver
from sdpcore import binom

n, d, Mfix, sense = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
work = sys.argv[5]
only = int(sys.argv[6]) if len(sys.argv) > 6 else None  # objective: A_only instead of sum j A_j
tproj = int(os.environ.get('TPROJ', '0'))  # if >0: objective sum_j A_j C(n-j,t)

def extra(P, X, n, d):
    # M >= Mfix : c0 + sum c_v x_v - Mfix >= 0   (c is the size objective)
    co = dict(P.c)
    P.add_lin(P.c0 - Mfix, co)
    P.add_lin(Mfix - P.c0, {v: -a for v, a in co.items()})
    # new objective: +/- sum_j j*C(n,j) x_{j,0,0}
    newc = {}
    for j in range(d, n + 1):
        v = X(j, 0, 0)
        if v is not None and v != 'one':
            if only is not None and j != only: continue
            wj = 1 if only is not None else j
            if tproj: wj = binom(n - j, tproj)
            newc[v] = newc.get(v, 0) + Fraction(wj * binom(n, j)) * (1 if sense == 'max' else -1)
    P.c = newc
    P.c0 = Fraction(0)

P = schrijver.build(n, d, extra=extra, name='avg_%d_%d_%d_%s_%s_t%d' % (n, d, Mfix, sense, only, tproj))
os.makedirs(work, exist_ok=True)
info = P.solve(work)
print(n, d, Mfix, sense, info.get('phase.value'), info.get('objValPrimal'), info.get('objValDual'), flush=True)
