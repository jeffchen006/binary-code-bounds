#!/usr/bin/env python3
"""Rebuild a gms.py exclusion SDP (parity or B1 rows, max s), read the SDPA result, round the dual to a rational
certificate, check it exactly (Problem.check_certificate: bound on s), and write a JSON certificate for the
stdlib checker check_cert2.py.
usage: export_cert.py n d M ub sdpa_out json_out (--parity | --b1 RMAX) [--K 200]

JSON (unscaled block coordinates):
  blocks: {'type': 'psi'|'MS', 'm','k','l', 'rows': [[i,i'],...], 'V': int matrix}  or  {'type': 'Z', 'rows': ['e', a...], 'V'}
          Y_b = V V^T / 4^K.  S2-reduced M_S blocks of the solver are mapped back to the full (m,k,l) block:
          the reduced block is T^T B T for an integer matrix T (columns = sigma-orbit sums), so <Y_red, T^T B T> = <T Y_red T^T, B>.
  lin:    {'row': 'size_ge'|'size_le'|'parity' (k)|'b1' (r, blocks), 'z': 'p/q'}  multiplier of the UNNORMALIZED row, s = 1:
          size_ge: y1 - M >= 0;  size_le: M - y1 >= 0;
          parity k: M K_k(0) + sum_i C(n,i) K_k(i) y2(i) - C(n,k) >= 0                 (M odd)
          b1 (r, blocks): sum_k m_k/C(n,k) [M K_k(0) + sum_i C(n,i) K_k(i) y2(i)] - rho(2^r - rho) >= 0,
                          V = span of indicator vectors of disjoint blocks, m_k = #nonzero words of weight k, rho = M mod 2^r.
"""
import sys, os, json, time, argparse
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gms

ap = argparse.ArgumentParser()
ap.add_argument('n', type=int); ap.add_argument('d', type=int); ap.add_argument('M', type=int); ap.add_argument('ub', type=int)
ap.add_argument('out'); ap.add_argument('json')
ap.add_argument('--parity', action='store_true'); ap.add_argument('--b1', type=int, default=None)
ap.add_argument('--mod4min', action='store_true')
ap.add_argument('--cw', action='store_true'); ap.add_argument('--mel', action='store_true'); ap.add_argument('--best', action='store_true')
ap.add_argument('--K', type=int, default=200)
a = ap.parse_args()
n, d, M, K = a.n, a.d, a.M, a.K
t0 = time.time()
gms.META.clear(); gms.LINMETA.clear()
P = gms.build(n, d, Fraction(a.ub), verbose=False)
if a.parity:
    gms.add_parity(P, n, d, M)
if a.b1:
    gms.add_b1(P, n, d, M, a.b1)
if a.mod4min:
    gms.add_mod4min(P, n, d, M)
if a.cw:
    gms.add_cw(P, n, d, M)
if a.mel:
    gms.add_mel(P, n, d, M)
if a.best:
    gms.add_best(P, n, d, M)
if not (a.parity or a.b1 or a.mod4min):
    raise SystemExit('need --parity, --b1 or --mod4min')
assert len(gms.META) == len(P.blocks), (len(gms.META), len(P.blocks))
assert len(gms.LINMETA) == len(P.lin), (len(gms.LINMETA), len(P.lin))
print('built', len(P.vars), len(P.blocks), len(P.lin), '%.0fs' % (time.time() - t0), flush=True)
info = P.read_result(a.out)
print('SDPA', info['phase.value'], info['objValPrimal'], info['objValDual'], flush=True)
cert = P.make_certificate(info, K=K)
b = P.check_certificate(cert)
print('EXACT bound on s (solver model): %.12f' % float(b), 'EXCLUDED' if b < 1 else 'NOT excluded', '%.0fs' % (time.time() - t0), flush=True)

S = 2 ** K
def toint(x):
    y = x * S
    assert y.denominator == 1
    return int(y)

blocks = []
ms2 = {}
for meta, V, blk in zip(gms.META, cert['V'], P.blocks):
    tag = meta[0]
    if tag == 'Z':
        blocks.append({'type': 'Z', 'rows': meta[4], 'V': [[toint(x) for x in row] for row in V]})
    elif tag in ('BLK', 'MS'):
        _, m, k, l, rows = meta
        assert len(rows) == blk['size'] == len(V)
        blocks.append({'type': 'psi' if tag == 'BLK' else 'MS', 'm': m, 'k': k, 'l': l, 'rows': [list(r) for r in rows],
                       'V': [[toint(x) for x in row] for row in V]})
    elif tag == 'MS2':
        _, m, k, l, sgn, cl = meta
        assert len(cl) == blk['size'] == len(V)
        ms2.setdefault((m, k, l), []).append((cl, V))
    else:
        raise ValueError(tag)
for (m, k, l), parts in sorted(ms2.items()):
    rows = sorted(set(tuple(r) for cl, V in parts for c in cl for (r, coef) in c))
    pos = {r: p for p, r in enumerate(rows)}
    Vf = [[] for _ in rows]
    for cl, V in parts:
        ncol = len(V[0]) if V else 0
        T = [[0] * len(cl) for _ in rows]
        for ci, c in enumerate(cl):
            for (r, coef) in c:
                T[pos[tuple(r)]][ci] += coef
        for i in range(len(rows)):
            for j in range(ncol):
                Vf[i].append(toint(sum((T[i][ci] * V[ci][j] for ci in range(len(cl))), Fraction(0))))
    blocks.append({'type': 'MS', 'm': m, 'k': k, 'l': l, 'rows': [list(r) for r in rows], 'V': Vf})

lin = []
for (meta, mx), z in zip(gms.LINMETA, cert['z']):
    assert z >= 0
    zz = z / mx      # multiplier of the unnormalized row
    if meta[0] in ('size_ge', 'size_le'):
        lin.append({'row': meta[0], 'z': str(zz)})
    elif meta[0] == 'parity':
        lin.append({'row': 'parity', 'k': meta[1], 'z': str(zz)})
    elif meta[0] == 'best':
        lin.append({'row': 'best', 'k': meta[1], 'Ap': meta[2], 'z': str(zz)})
    elif meta[0] == 'mel19':
        lin.append({'row': 'mel19', 'q': meta[1], 'z': str(zz)})
    elif meta[0] == 'mel21':
        lin.append({'row': 'mel21', 'i': meta[1], 'S': meta[2], 'T': meta[3], 'z': str(zz)})
    elif meta[0] == 'cw':
        lin.append({'row': 'cw', 'a': meta[1], 'U': meta[2], 'z': str(zz)})
    elif meta[0] == 'mod4min':
        lin.append({'row': 'mod4min', 'k': meta[1], 'z': str(zz)})
    elif meta[0] == 'b1':
        lin.append({'row': 'b1', 'r': meta[1], 'blocks': meta[2], 'z': str(zz)})
    else:
        raise ValueError(meta)
out = {'problem': {'n': n, 'd': d, 'M': M, 'ub': a.ub, 'rows': '+'.join(x for x, f in (('parity', a.parity), ('b1', a.b1), ('mod4min', a.mod4min), ('cw', a.cw), ('mel', a.mel), ('best', a.best)) if f)},
       'K': K, 'note': 'Y_b = V V^T / 4^K (unscaled GMS block coordinates); linear multipliers z (exact rationals) >= 0',
       'blocks': blocks, 'lin': lin}
json.dump(out, open(a.json, 'w'))
print('JSON written:', a.json, 'blocks', len(blocks), 'lin', len(lin), '%.0fs' % (time.time() - t0))
