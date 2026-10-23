"""Export a gmsb1.py certificate (pickle) to certio-GMS JSON with row tags, rebuilding the problem in the same row order.
Row tags: 'gms' (rows of gms0.build_full: Schrijver basic, quadruple rows; checked by check_gms_cert.py), 'parity',
'oa_blocks' (block subspace: list of block sizes), 'oa_cols' (general V: r, column multiplicities), 'box' (x <= p).
Classical rows (--lin runs) are not supported (export only no-lin runs).
usage: export_gmsb1.py n d M case w tag round out.json"""
import sys, os, pickle, json
sys.dont_write_bytecode = True
sys.path.insert(1, 'agents/novaB')
from fractions import Fraction
import mpmath
import gms0, parity, split, sdpcore, certio
from oa import add_oa_row
from t3 import partitions, block_wd
from cp import wd_of

n, d, M = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
case = sys.argv[4]; w = None if sys.argv[5] == 'None' else int(sys.argv[5])
tag, rnd, outp = sys.argv[6], int(sys.argv[7]), sys.argv[8]
rmax = int(os.environ.get('RMAX_EXPORT', '99'))
off = dict(best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)
base = 'gb_%d_%d_%d_%s_%s_l0%s' % (n, d, M, case, w, tag)
work = 'agents/sdp/runs/gmsb1'
pk = pickle.load(open(os.path.join(work, base + '_r%d.cert.pkl' % rnd), 'rb'))
cert, pool = pk['cert'], pk['pool']

P = gms0.build_full(n, d, 1, opts=dict(off), name=base, light=True, quad_lin=True)
meta = [('gms',)] * len(P.lin)
Q = split.build(n, d, 0, opts=dict(off, lin_basic=False), name='q')
if case in ('odd', '2mod4'):
    Q2 = parity.build(n, d, 0, M, case, w, opts=dict(off, lin_basic=False))
    Q.lin += Q2.lin
    meta += [('parity',)] * len(Q2.lin)
# block rows in the order of t3.block_rows (rmin = 2 for odd, else 1)
seen = set()
rmin = 2 if case == 'odd' else 1
for tot in range(1, n + 1):
    for part in partitions(tot):
        r = len(part)
        if r < rmin or M % 2 ** r == 0 or r > rmax:
            continue
        wd = block_wd(part); key = (r, tuple(sorted(wd.items())))
        if key in seen:
            continue
        seen.add(key)
        if add_oa_row(Q, n, M, wd, r):
            meta.append(('oa_blocks', list(part)))
for r, mult in pool:
    wd = wd_of(r, mult)
    if wd is not None and add_oa_row(Q, n, M, wd, r):
        meta.append(('oa_cols', r, list(mult)))
P.lin = P.lin + Q.lin
nl = len(P.lin)
P = gms0.homogenize(P, n, xle_p=True)
meta += [('box',)] * (len(P.lin) - nl)
assert len(meta) == len(P.lin) == len(cert['z']), (len(meta), len(P.lin), len(cert['z']))
for v in P.vars:
    P.ub[v] = Fraction(M, 2 ** n)
b = P.check_certificate(cert)
print('internal exact bound %.9f (M = %d)' % (float(b), M), flush=True)
out = certio.export_gms(P, cert, 80, outp, {'n': n, 'd': d, 'M': M, 'case': case, 'w': w, 'b1': True})
j = 0
for z, m in zip(cert['z'], meta):
    if z == 0:
        continue
    if m[0] == 'oa_blocks':
        out['lin'][j]['oa_blocks'] = m[1]
    elif m[0] == 'oa_cols':
        out['lin'][j]['oa_cols'] = {'r': m[1], 'mult': m[2]}
    out['lin'][j]['cls'] = m[0]
    j += 1
json.dump(out, open(outp, 'w'))
print('wrote', outp, 'rows', len(out['lin']))
