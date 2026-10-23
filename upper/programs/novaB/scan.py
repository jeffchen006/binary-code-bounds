"""scan M downward from M0 until the SDP (plain Schrijver + block OA rows) is feasible. usage: scan.py n d M0 [mode]"""
import sys, os
from oa import *
from t3 import block_rows
import sdpcore
n, d, M = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); mode = sys.argv[4] if len(sys.argv) > 4 else 'plain'
while True:
    opts = dict(lin_basic=True, best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False) if mode == 'plain' else None
    P = split.build(n, d, 0, opts=opts, name='scan_%d_%d_%d_%s' % (n, d, M, mode))
    nr = sum(add_oa_row(P, n, M, wd, r) for r, wd in block_rows(n, M, rmin=1))
    info = sdpcore.solve_csdp(P, "agents/novaB/runs/scan", threads=1, cpus=os.environ.get("CPU"))
    ex = info['bound'] is not None and info['bound'] < M
    print(n, d, M, M % 64, nr, info['phase.value'][:10], info['bound'], 'EXCL' if ex else 'feas', flush=True)
    if not ex: break
    M -= 1
