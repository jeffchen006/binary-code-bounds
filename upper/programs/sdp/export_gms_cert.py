"""Rebuild the gmspar problem (same args/env) and export its pickled certificate to JSON.
usage: export_gms_cert.py n d T case w lin light TAG PURE out.json"""
import sys, os, pickle, mpmath
import gms0, parity, sdpcore, certio
n, d, T = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); case = sys.argv[4]; w = None if sys.argv[5] == 'None' else int(sys.argv[5])
lin = sys.argv[6] == '1'; light = sys.argv[7] != '0'; tag = sys.argv[8]; pure = sys.argv[9] == '1'; out = sys.argv[10]
off = dict(best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)
opts = {} if lin else off
if pure:
    opts = dict(opts, schrijver=False, lin_basic=False)
name = 'gp_%d_%d_%d_%s_%s_l%d_g%d%s' % (n, d, T, case, w, lin, light, tag)
P = gms0.build_full(n, d, 1, opts=opts, name=name, light=light, quad_lin=not pure)
if case != 'base':
    Q = parity.build(n, d, 0, T, case, w, opts=dict(off, lin_basic=False))
    P.lin = P.lin + Q.lin
P = gms0.homogenize(P, n, xle_p=not pure)
pk = 'agents/sdp/runs/gmspar/%s.cert.pkl' % name
if os.path.exists(pk):
    cert = pickle.load(open(pk, 'rb'))['cert']
else:
    # rebuild scaling exactly as gmspar did, then read the solver output
    P.vscale = {v: mpmath.mpf(T) / 2 ** n for v in P.vars}
    if os.environ.get('EQ', '1') == '1':
        sdpcore.equilibrate(P)
    P.write_sdpa('/tmp/work')
    info = P.read_result('agents/sdp/runs/gmspar/%s.out' % name)
    cert = P.make_certificate(info, K=80)
from fractions import Fraction
for v in P.vars:
    P.ub[v] = Fraction(T, 2 ** n)
b = P.check_certificate(cert)
print('internal bound', float(b))
prob = {'n': n, 'd': d, 'M': (None if case == 'base' else T), 'case': (None if case == 'base' else case), 'w': w}
if case == 'base':
    prob['target'] = T
certio.export_gms(P, cert, 80, out, prob)
print('exported', out)
