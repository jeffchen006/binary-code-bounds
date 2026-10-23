#!/usr/bin/env python3
"""Cheap screen: Schrijver SDP (sdp's split.build, n1=0) + all B1 block rows (r>=1, incl. parity) for size exactly M;
max |C| < M excludes M.  CSDP (double) -> numeric only.  usage: schr_b1.py n d M [M ...]   (CPUS env)"""
import sys, os, time
sys.path.insert(0, 'agents/novaB')
sys.path.insert(0, 'agents/sdp')
CPUS = os.environ.get('CPUS')
os.environ['SDP_CPUS'] = CPUS
import split, sdpcore
from oa import add_oa_row, betas
from t3 import block_rows
n, d = int(sys.argv[1]), int(sys.argv[2])
for M in map(int, sys.argv[3:]):
    t = time.time()
    P = split.build(n, d, 0, name='sb1_%d_%d_%d' % (n, d, M))
    nr = sum(add_oa_row(P, n, M, wd, r) for r, wd in block_rows(n, M, rmin=1))
    info = sdpcore.solve_csdp(P, 'agents/exact/data/sb1', threads=1, cpus=CPUS)
    A, b = betas(P, n, M, info['x'])
    print(n, d, M, 'rows', nr, info['phase.value'], 'max|C| =', info['bound'], 'EXCL' if info['bound'] < M else 'feasible',
          '%.0fs' % (time.time() - t), {k: round(v, 2) for k, v in b.items() if 0 < k < n and abs(v) > 1e-4}, flush=True)
