#!/usr/bin/env python3
"""sdp's homogenized GMS (gmspar.py pipeline: gms0.build_full + parity rows + homogenize + vscale + equilibrate)
plus B1 block rows (r <= rmax) for size exactly T.  Objective max |C| = 2^n p;  max < T excludes T.
usage: gp_b1.py n d T case w rmax [lin]   (case odd|2mod4|base; SDP_CPUS, TH env)"""
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b1
sys.path.insert(0, b1.NOVAB)
from t3 import block_rows
import sdpcore
n, d, T = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
case = sys.argv[4]; w = None if sys.argv[5] == 'None' else int(sys.argv[5]); rmax = int(sys.argv[6])
lin = len(sys.argv) > 7 and sys.argv[7] == '1'
name = 'gpb1_%d_%d_%d_%s_%s_r%d_l%d' % (n, d, T, case, w, rmax, lin)
t = time.time()
P = b1.build_gp(n, d, T, case, w, lin=lin, light=True, name=name)
nb = 0
for r, wd in block_rows(n, T, rmin=2):
    if r <= rmax:
        nb += b1.add_b1(P, n, T, wd, r)
P = b1.homog(P, n, T)
print(name, 'vars', len(P.vars), 'blocks', len(P.blocks), 'lin', len(P.lin), 'B1 rows', nb, '%.0fs' % (time.time() - t), flush=True)
work = 'agents/exact/data/gpb1'
info = P.solve(work, threads=int(os.environ.get('TH', '3')))
res = '%s %s num=%s pnum=%s' % (name, info['phase.value'], info['bound'], info.get('pbound'))
try:
    cert = P.make_certificate(info, K=80); bnd = P.check_certificate(cert)
    pickle.dump({'cert': cert, 'bound': bnd}, open(os.path.join(work, name + '.cert.pkl'), 'wb'))
    res += ' rig=%.6f %s' % (float(bnd), 'EXCL' if bnd < T else '')
except Exception as e:
    res += ' certfail %r' % e
print(res, '%.0fs' % (time.time() - t), flush=True)
