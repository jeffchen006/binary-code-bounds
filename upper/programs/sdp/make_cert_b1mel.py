"""Exact dual certificate: Schrijver SDP + MEL row (A_{n-h} + q sum_{i<h} A_{n-i} <= q, q = floor(2n/d))
+ B1 coset-count rows (block subspaces + general V from a pool) excludes |C| = M.
usage: make_cert_b1mel.py n d M outdir [pool.json] [--nomel] [--sdpa]"""
import sys, os, json
pass
sys.dont_write_bytecode = True
sys.path.insert(1, 'agents/novaB')
from fractions import Fraction
import split, sdpcore, certio
from sdpcore import binom
from oa import add_oa_row, X0forms
from t3 import partitions, block_wd
from cp import wd_of

n, d, M, outdir = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
poolf = sys.argv[5] if len(sys.argv) > 5 and not sys.argv[5].startswith('--') else None
usemel = '--nomel' not in sys.argv
K = 60
os.makedirs(outdir, exist_ok=True)
opts = dict(lin_basic=True, best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)
P = split.build(n, d, 0, opts=opts, name='certm_%d_%d_%d' % (n, d, M))
meta = [None] * len(P.lin)
if usemel:
    h = d // 2; q = (2 * n) // d
    X0 = X0forms(P, n)
    co = {}
    for a, c in [(n - h, 1)] + [(n - i, q) for i in range(h)]:
        for v, f in X0.get(a, {}).items():
            co[v] = co.get(v, 0) - c * f
    mx = max([abs(c) for c in co.values()] + [Fraction(q)])
    P.lin.append((Fraction(q) / mx, {v: c / mx for v, c in co.items()}))
    meta.append(('mel',))
seen = set()
TGT = int(os.environ.get('TARGET', '0'))
for tot in (range(0) if M == 0 else range(1, n + 1)):
    for part in partitions(tot):
        r = len(part)
        if M % 2 ** r == 0: continue
        wd = block_wd(part); key = (r, tuple(sorted(wd.items())))
        if key in seen: continue
        seen.add(key)
        if add_oa_row(P, n, M, wd, r): meta.append(('blocks', list(part)))
for r, mult in (json.load(open(poolf)) if poolf else []):
    wd = wd_of(r, mult)
    if wd is not None and add_oa_row(P, n, M, wd, r): meta.append(('cols', r, list(mult)))
work = 'agents/sdp/runs/cert'
if '--sdpa' in sys.argv:
    K = 80
    if os.environ.get('EQ', '1') == '1':
        sdpcore.equilibrate(P)
    info = P.solve(work, threads=1)
else:
    info = sdpcore.solve_csdp(P, work, threads=1)
cert = P.make_certificate(info, K=K)
b = P.check_certificate(cert)
print('num', info['bound'], 'rigorous', float(b), 'EXCLUDED' if b < (M or TGT) else 'NOT EXCLUDED', flush=True)
if M == 0:
    path = os.path.join(outdir, 'base_n%d_d%d_mel.json' % (n, d))
    out = certio.export(P, cert, K, path, {'n': n, 'd': d, 'target': TGT, 'mel_row': usemel})
else:
    path = os.path.join(outdir, 'excl_n%d_d%d_M%d_b1mel.json' % (n, d, M))
    out = certio.export(P, cert, K, path, {'n': n, 'd': d, 'M': M, 'case': 'oa', 'mel_row': usemel})
j = 0
for z, m in zip(cert['z'], meta):
    if z == 0: continue
    if m:
        if m[0] == 'mel': out['lin'][j]['mel'] = True
        elif m[0] == 'blocks': out['lin'][j]['oa_blocks'] = m[1]
        else: out['lin'][j]['oa_cols'] = {'r': m[1], 'mult': m[2]}
    j += 1
json.dump(out, open(path, 'w'))
print('wrote', path, 'rows used', len(out['lin']), 'mel used', any(r.get('mel') for r in out['lin']),
      'oa rows used', sum(1 for r in out['lin'] if r.get('oa_blocks') or r.get('oa_cols')))
