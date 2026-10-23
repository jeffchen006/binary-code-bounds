"""scan M downward: plain Schrijver + block OA rows + column cap row (needs proven A(n-1,d) <= Ap):
   M - Ap <= w_i <= Ap  =>  |S_{e_i}| <= 2Ap - M  => sum_{|u|=1} S_u^2 <= n (2Ap-M)^2
   i.e.  sum_i A_i K_1(i) <= n (2Ap - M)^2 / M   (A_0 = 1).
usage: scan2.py n d M0 Ap"""
import sys, os
from fractions import Fraction
from oa import *
from t3 import block_rows
import sdpcore
n, d, M, Ap = map(int, sys.argv[1:5])
while True:
    opts = dict(lin_basic=True, best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)
    P = split.build(n, d, 0, opts=opts, name='scan2_%d_%d_%d' % (n, d, M))
    nr = sum(add_oa_row(P, n, M, wd, r) for r, wd in block_rows(n, M, rmin=1))
    c0, f = shell_form(P, n, 1)
    cap = Fraction(n * (2 * Ap - M) ** 2, M)
    # cap - c0 - f(x) >= 0
    P.lin.append((cap - c0, {v: -c for v, c in f.items()}))
    info = sdpcore.solve_csdp(P, "agents/novaB/runs/scan2", threads=1, cpus=os.environ.get("CPU"))
    ex = info['bound'] is not None and info['bound'] < M
    print(n, d, M, 'Ap', Ap, nr, info['phase.value'][:10], info['bound'], 'EXCL' if ex else 'feas', flush=True)
    if not ex: break
    M -= 1
