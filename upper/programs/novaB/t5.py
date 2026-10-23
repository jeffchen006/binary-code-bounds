import sys, os
from fractions import Fraction
from oa import *
from t3 import block_rows
import sdpcore
n, d, M, Ap = map(int, sys.argv[1:5])
opts = dict(lin_basic=True, best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)
P = split.build(n, d, 0, opts=opts, name='t5_%d_%d_%d' % (n, d, M))
nr = sum(add_oa_row(P, n, M, wd, r) for r, wd in block_rows(n, M, rmin=1))
c0, f = shell_form(P, n, 1)
P.lin.append((Fraction(n * (2 * Ap - M) ** 2, M) - c0, {v: -c for v, c in f.items()}))
info = sdpcore.solve_csdp(P, "agents/novaB/runs/t5", threads=1, cpus=os.environ.get("CPU"))
A, b = betas(P, n, M, info['x'])
print(info['bound']); print('A', {a: round(v, 3) for a, v in A.items() if abs(v) > 1e-6}); print('beta', {k: round(v, 3) for k, v in b.items() if abs(v) > 1e-4})
