import sys, time, os
from math import comb
from oa import *
from t3 import block_rows
import sdpcore, parity
n, d = int(sys.argv[1]), int(sys.argv[2]); plain = sys.argv[3] == 'plain'
for M in map(int, sys.argv[4:]):
    opts = dict(lin_basic=True, best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False) if plain else None
    P = split.build(n, d, 0, opts=opts, name='oab_%d_%d_%d' % (n, d, M))
    nr = 0
    for r, wd in block_rows(n, M, rmin=1):
        nr += add_oa_row(P, n, M, wd, r)
    info = sdpcore.solve_csdp(P, "agents/novaB/runs/t4", threads=1, cpus=os.environ.get("CPU"))
    print(M, M % 64, nr, info['phase.value'][:12], info['bound'], 'EXCL' if info['bound'] < M else 'feas', flush=True)
