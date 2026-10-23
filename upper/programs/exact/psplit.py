#!/usr/bin/env python3
"""Split (u-split) Schrijver SDP for M = 2 mod 4 with sum-weight w (sdp's parity.build_split2mod4, read-only import).
usage: psplit.py n d M w   (SDP_CPUS env = CPU list)"""
import sys, os, time
sys.path.insert(0, 'agents/sdp')
import parity
n, d, M, w = map(int, sys.argv[1:5])
t = time.time()
P = parity.build_split2mod4(n, d, M, w, name='psplit_%d_%d_%d_w%d' % (n, d, M, w))
print('vars', len(P.vars), 'blocks', len(P.blocks), 'maxblk', max(b['size'] for b in P.blocks), 'lin', len(P.lin), '%.0fs' % (time.time() - t), flush=True)
info = P.solve('agents/exact/data/psplit', threads=int(os.environ.get('TH', '3')))
print(n, d, M, w, info['phase.value'], info['bound'], 'EXCL' if info['bound'] < M else 'feasible', '%.0fs' % (time.time() - t), flush=True)
