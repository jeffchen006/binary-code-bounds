"""Base cert: plain Schrijver SDP (Schrijver's own constraints) => |C| <= bound < target. usage: make_base.py n d target outdir"""
import sys, os
from oa import *
import sdpcore, certio
n, d, T, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
opts = dict(lin_basic=True, best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)
P = split.build(n, d, 0, opts=opts, name='base_%d_%d' % (n, d))
info = sdpcore.solve_csdp(P, 'agents/novaB/runs/base', threads=1, cpus=os.environ.get('CPU'))
cert = P.make_certificate(info, K=60); b = P.check_certificate(cert)
print('num', info['bound'], 'rigorous', float(b), 'OK' if b < T else 'FAIL')
os.makedirs(out, exist_ok=True)
certio.export(P, cert, 60, os.path.join(out, 'base_n%d_d%d.json' % (n, d)), {'n': n, 'd': d, 'M': None, 'target': T})
