"""Exact certificate: Schrijver SDP + Best rows (premises A(n-k,d) <= P_k, given) + B1 rows (blocks + pool) excludes |C| = M.
Best row (k = 1..n-1): sum_a C(n-a,k) A_a <= C(n,k) A(n-k,d). Proof: count triples (c, c', K), K a k-set of coordinates,
c, c' in C equal on K: for fixed (c, K) the codewords equal to c on K, with K deleted, form a code of length n-k and
distance >= d, so the count is <= M C(n,k) A(n-k,d); for fixed (c, c') there are C(n - d(c,c'), k) sets K, so the count is
M sum_a A_a C(n-a,k).
usage: make_cert_best.py n d M outdir pool.json "m1,P1 m2,P2 ..."  (premises A(m,d) <= P; other k use A(m,d) <= 2^(m-1))"""
import sys, os, json
pass
sys.dont_write_bytecode = True
sys.path.insert(1, 'agents/novaB')
from fractions import Fraction
import bounds
n, d, M, outdir, poolf = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
prem = {int(a): int(b) for a, b in (kv.split(',') for kv in sys.argv[6].split())}
bounds.OVERRIDE.clear()
for m in range(1, n):
    bounds.OVERRIDE[(m, d)] = prem.get(m, 2 ** (m - 1))
import split, sdpcore, certio
from sdpcore import binom
from oa import add_oa_row, X0forms
from t3 import partitions, block_wd
from cp import wd_of

K = 80
os.makedirs(outdir, exist_ok=True)
opts = dict(lin_basic=True, best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)
P = split.build(n, d, 0, opts=opts, name='certb_%d_%d_%d' % (n, d, M))
meta = [None] * len(P.lin)
X0 = X0forms(P, n)
for k in range(1, n):
    if os.environ.get('ONLYPREM') == '1' and (n - k) not in prem:
        continue
    ub = bounds.A_upper(n - k, d)
    co = {}
    for a, f in X0.items():
        for v, c in f.items():
            co[v] = co.get(v, 0) - binom(n - a, k) * c
    g0 = Fraction(binom(n, k) * ub - binom(n, k))          # a = 0 term: C(n,k) * A_0
    co = {v: c for v, c in co.items() if c}
    if not co:
        continue
    mx = max([abs(c) for c in co.values()] + [abs(g0)])
    P.lin.append((g0 / mx, {v: c / mx for v, c in co.items()}))
    meta.append(('best', k, ub))
seen = set()
for tot in range(1, n + 1):
    for part in partitions(tot):
        r = len(part)
        if M % 2 ** r == 0: continue
        wd = block_wd(part); key = (r, tuple(sorted(wd.items())))
        if key in seen: continue
        seen.add(key)
        if add_oa_row(P, n, M, wd, r): meta.append(('blocks', list(part)))
for r, mult in json.load(open(poolf)):
    wd = wd_of(r, mult)
    if wd is not None and add_oa_row(P, n, M, wd, r): meta.append(('cols', r, list(mult)))
if os.environ.get('NOEQ') != '1':
    sdpcore.equilibrate(P)
info = P.solve('agents/sdp/runs/cert', threads=1)
cert = P.make_certificate(info, K=K)
b = P.check_certificate(cert)
print('num', info['bound'], 'exact', float(b), 'margin %.9f' % (M - float(b)), 'EXCLUDED' if b < M else 'NOT EXCLUDED', flush=True)
path = os.path.join(outdir, 'excl_n%d_d%d_M%d_best_b1.json' % (n, d, M))
used = {}
out = certio.export(P, cert, K, path, {'n': n, 'd': d, 'M': M, 'case': 'oa', 'best_premises': {}})
j = 0
for z, m in zip(cert['z'], meta):
    if z == 0: continue
    if m:
        if m[0] == 'best':
            out['lin'][j]['best'] = {'k': m[1], 'premise': m[2]}; used[n - m[1]] = m[2]
        elif m[0] == 'blocks': out['lin'][j]['oa_blocks'] = m[1]
        else: out['lin'][j]['oa_cols'] = {'r': m[1], 'mult': m[2]}
    j += 1
out['problem']['best_premises'] = {str(k): v for k, v in sorted(used.items())}
json.dump(out, open(path, 'w'))
print('wrote', path, 'best rows used (A(m,d) <= P):', out['problem']['best_premises'])
