"""Export tagged exact certificates (certio format, like UB_A_17_6_339) for the referee's spot-check of the P1 audit.
Rows carry 'cls' = basic | best18 | cw22 | mel19 | mel21 | best20 | kt23 | parity; the problem dict lists every constant
used: premises A(m,d) (Best rows), cw values A(n,d,w) (cw22, mel21, best20) and doubly-constant-weight T values (kt23).
usage: export_for_referee.py n d target [M case w]   (no M: base certificate, bound < target)"""
import sys, os, json
sys.dont_write_bytecode = True
pass
import bounds, split, parity, sdpcore, certio
from fractions import Fraction

OUTD = 'agents/sdp/certs_for_referee'
os.makedirs(OUTD, exist_ok=True)
n, d, target = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
M = int(sys.argv[4]) if len(sys.argv) > 4 else None
case = sys.argv[5] if len(sys.argv) > 5 else None
w = (None if sys.argv[6] == 'None' else int(sys.argv[6])) if len(sys.argv) > 6 else None

# self-consistent premises from the P1b pass
prem = {}
for l in open('agents/sdp/runs/audit/self.jsonl'):
    r = json.loads(l)
    if r['n'] == n and r['d'] == d:
        prem = {int(k): v for k, v in r['premises'].items()}
bounds.OVERRIDE.clear()
for m, v in prem.items():
    bounds.OVERRIDE[(m, d)] = v

# record constants
used = {'A': {}, 'cw': {}, 'T': {}}
_A, _cw, _T = bounds.A_upper, bounds.cw_upper, bounds.T_upper
def A_rec(m, dd):
    v = _A(m, dd); used['A']['%d,%d' % (m, dd)] = v; return v
def cw_rec(m, dd, ww):
    v = _cw(m, dd, ww); used['cw']['%d,%d,%d' % (m, dd, ww)] = v; return v
def T_rec(*a):
    v = _T(*a); used['T'][','.join(map(str, a))] = v; return v
bounds.A_upper, bounds.cw_upper, bounds.T_upper = A_rec, cw_rec, T_rec

CLS = ['best18', 'cw22', 'mel19', 'mel21', 'best20', 'kt23']
OFF = {c: False for c in CLS}
def key(g0, co):
    return (g0, frozenset(co.items()))
tag = {}
P0 = split.build(n, d, 0, opts=dict(OFF), name='b')
for g0, co in P0.lin:
    tag[key(g0, co)] = 'basic'
for c in CLS:
    Pc = split.build(n, d, 0, opts=dict(OFF, **{c: True}), name='c')
    for g0, co in Pc.lin:
        tag.setdefault(key(g0, co), c)
P = split.build(n, d, 0, opts={}, name='ref_%d_%d' % (n, d))
if M is not None:
    Q = parity.build(n, d, 0, M, case, w, opts={})
    for g0, co in Q.lin:
        tag.setdefault(key(g0, co), 'parity')
    P = Q
import mpmath
if os.environ.get('VS') == '1':
    P.vscale = {v: mpmath.mpf(target) / 2 ** n for v in P.vars}
if os.environ.get('EQ', '1') == '1':
    sdpcore.equilibrate(P)
info = P.solve('agents/sdp/runs/referee', threads=1)
cert = P.make_certificate(info, K=80)
b = P.check_certificate(cert)
prob = {'n': n, 'd': d, 'premises_A_m_d_le': {k: v for k, v in used['A'].items()}, 'cw_values': used['cw'],
        'dcw_T_values': used['T'], 'note': 'premises = self-consistent values (runs/audit/self.jsonl); classes per split.py'}
if M is None:
    prob['target'] = target
    path = os.path.join(OUTD, 'base_n%d_d%d.json' % (n, d))
else:
    prob.update({'M': M, 'case': case, 'w': w})
    path = os.path.join(OUTD, 'excl_n%d_d%d_M%d_%s_w%s.json' % (n, d, M, case, w))
out = certio.export(P, cert, 80, path, prob)
miss = 0; j = 0
for z, (g0, co) in zip(cert['z'], P.lin):
    if z == 0:
        continue
    t = tag.get(key(g0, co))
    if t is None:
        miss += 1; t = 'UNTAGGED'
    out['lin'][j]['cls'] = t; j += 1
json.dump(out, open(path, 'w'))
print('%s: exact bound %.6f target %s -> %s; rows %d untagged %d' % (os.path.basename(path), float(b), M or target,
      'OK' if b < (M or target) else 'NOT BELOW', len(out['lin']), miss), flush=True)
