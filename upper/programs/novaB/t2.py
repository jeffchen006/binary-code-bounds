import sys, time
from math import comb
from oa import *
import sdpcore
n, d = int(sys.argv[1]), int(sys.argv[2])
for M in map(int, sys.argv[3:]):
    P = split.build(n, d, 0, name='oa_%d_%d_%d' % (n, d, M))
    nr = 0
    for r in range(1, n + 1):
        nr += add_oa_row(P, n, M, {k: comb(r, k) for k in range(1, r + 1)}, r)
    info = sdpcore.solve_csdp(P, "agents/novaB/runs/t2", threads=1, cpus=None)
    A, b = betas(P, n, M, info['x'])
    print(M, M % 64, nr, info['phase.value'], info['bound'], {k: round(v, 2) for k, v in b.items() if 0 < k < n and abs(v) > 1e-4}, flush=True)
