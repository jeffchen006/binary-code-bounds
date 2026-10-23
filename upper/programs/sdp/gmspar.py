"""GMS 4-point SDP (homogenized, max 2^n p) + optional Best/MEL/KT rows + optional parity rows for exact size T.
usage: gmspar.py n d T case w [lin=0/1] [light=0/1]   (case in base|odd|2mod4; w ignored unless 2mod4)"""
import sys, time, os, pickle
from fractions import Fraction
import gms0, parity, sdpcore
n, d, T = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
case = sys.argv[4]; w = None if sys.argv[5] == 'None' else int(sys.argv[5])
lin = len(sys.argv) > 6 and sys.argv[6] == '1'
light = not (len(sys.argv) > 7 and sys.argv[7] == '0')
off = dict(best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)
opts = {} if lin else off
t = time.time()
name = 'gp_%d_%d_%d_%s_%s_l%d_g%d%s' % (n, d, T, case, w, lin, light, os.environ.get('TAG', ''))
pure = os.environ.get('PURE', '0') == '1'
if pure:
    opts = dict(opts, schrijver=False, lin_basic=False)
P = gms0.build_full(n, d, 1, opts=opts, name=name, light=light, quad_lin=not pure)
if case != 'base':
    # parity rows (normalized form, homogenized below): sum_i A_i K_k(i) >= c_k / T
    Q = parity.build(n, d, 0, T, case, w, opts=dict(off, lin_basic=False))
    P.lin = P.lin + Q.lin      # Q.lin = only parity rows (other rows disabled)
P = gms0.homogenize(P, n, xle_p=not pure)
import mpmath
P.vscale = {v: mpmath.mpf(T) / 2 ** n for v in P.vars}
if os.environ.get('EQ', '1') == '1':
    sdpcore.equilibrate(P)
print(name, 'vars', len(P.vars), 'blocks', len(P.blocks), 'lin', len(P.lin), '%.0fs' % (time.time() - t), flush=True)
work = 'agents/sdp/runs/gmspar'
info = P.solve(work, threads=int(os.environ.get('TH', '4')))
res = '%s %s num=%s pnum=%s' % (name, info['phase.value'], info['bound'], info.get('pbound'))
try:
    cert = P.make_certificate(info, K=80); b = P.check_certificate(cert)
    pickle.dump({'cert': cert, 'bound': b}, open(os.path.join(work, name + '.cert.pkl'), 'wb'))
    res += ' rig=%.6f %s' % (float(b), 'EXCL' if (case != 'base' and b < T) else '')
except Exception as e:
    res += ' certfail %r' % e
res += ' %.0fs' % (time.time() - t)
print(res, flush=True)
open('agents/sdp/gms_results.txt', 'a').write(time.strftime('%m-%d %H:%M ') + res + '\n')
