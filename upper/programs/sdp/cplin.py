"""Schrijver 3-point SDP + classical linear rows (Best/MEL/KT/cw, split.py) + B1 coset-count rows
(novaB: block subspaces + general-V separation) + optional cap row (premise A(n-1,d) <= Ap) + parity rows.
Team-proven upper bounds (audited) are used as premises in the Best rows via bounds.OVERRIDE.
usage: cplin.py n d M Ap rounds [case w]"""
import sys, os, random, json, time
pass
sys.dont_write_bytecode = True
sys.path.insert(1, 'agents/novaB')
from fractions import Fraction
import bounds
bounds.OVERRIDE.update({(17, 4): 3272, (18, 4): 6543, (19, 4): 13085, (21, 4): 43681, (25, 4): 599174,
                        (17, 6): 339, (18, 6): 672, (22, 6): 6939, (24, 6): 24102, (21, 10): 46,
                        (26, 4): 1198348, (27, 4): 2396696, (26, 6): 84243, (27, 6): 157281, (28, 6): 291265,
                        (28, 8): 32148})
for kv in os.environ.get('OVR', '').split():
    a_, b_, c_ = kv.split(',')
    bounds.OVERRIDE[(int(a_), int(b_))] = Fraction(c_)
if os.environ.get('CWLIVE') == '1':
    import cwlive
    cwlive.apply()
import split, parity, sdpcore
from oa import add_oa_row, betas, shell_form
from t3 import block_rows
from cp import separate, wd_of

n, d, M, Ap, rounds = map(int, sys.argv[1:6])
case = sys.argv[6] if len(sys.argv) > 6 else 'none'
w = int(sys.argv[7]) if len(sys.argv) > 7 and sys.argv[7] != 'None' else None
random.seed(12345)
work = 'agents/sdp/runs/cplin'
os.makedirs(work, exist_ok=True)
tag = 'cl_%d_%d_%d_%d_%s_%s' % (n, d, M, Ap, case, w)
poolf0 = os.environ.get('POOL')
poolf = os.path.join(work, tag + '.pool.json')
pool = json.load(open(poolf0 or poolf)) if os.path.exists(poolf0 or poolf) else []
if os.environ.get('TAG'):
    tag += os.environ['TAG']; poolf = os.path.join(work, tag + '.pool.json')


ALL = ['best18', 'cw22', 'mel19', 'mel21', 'best20', 'kt23']
ON = os.environ.get('LINOPTS', ','.join(ALL)).split(',')
LOPTS = {k: (k in ON) for k in ALL}


def build():
    P = split.build(n, d, 0, opts=LOPTS, name=tag)
    if case in ('odd', '2mod4', '0mod4'):
        off = dict(best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False, lin_basic=False)
        P.lin += parity.build(n, d, 0, M, case, w, opts=off).lin
    nb = 0
    for r, wd in block_rows(n, M, rmin=1):
        nb += add_oa_row(P, n, M, wd, r)
    for r, mult in pool:
        wd = wd_of(r, mult)
        if wd is not None:
            nb += add_oa_row(P, n, M, wd, r)
    if Ap:
        c0, f = shell_form(P, n, 1)
        P.lin.append((Fraction(n * (2 * Ap - M) ** 2, M) - c0, {v: -c for v, c in f.items()}))
    return P, nb


for it in range(rounds):
    t = time.time()
    P, nb = build()
    P.name = tag + '_r%d' % it
    if os.environ.get('VSCALE') == '1':
        import mpmath
        P.vscale = {v: mpmath.mpf(M) / 2 ** n for v in P.vars}
    if os.environ.get('USE_SDPA') == '1':
        if os.environ.get('NOEQ') != '1':
            sdpcore.equilibrate(P)
        info = P.solve(work, threads=1)
        try:
            cert = P.make_certificate(info, K=80); rig = P.check_certificate(cert)
            print('  rigorous bound %.6f' % float(rig), flush=True)
            import pickle; pickle.dump({'cert': cert, 'pool': pool}, open(os.path.join(work, P.name + '.cert.pkl'), 'wb'))
        except Exception as e:
            print('  certfail', repr(e), flush=True)
    else:
        info = sdpcore.solve_csdp(P, work, threads=1)
    xs = info['x']
    if getattr(P, 'vscale', None):
        xs = [P.vscale[v] * y for v, y in zip(P.vars, xs)]
    A, beta = betas(P, n, M, xs)
    print('round', it, 'lin', len(P.lin), 'B1', nb, 'pool', len(pool), info['phase.value'], 'bound', info['bound'],
          'EXCL' if info['bound'] < M else 'feas', '%.0fs' % (time.time() - t), flush=True)
    if info['bound'] < M:
        break
    cuts = separate(beta, n, M)
    keys = {(r, tuple(m)) for r, m in pool}; new = 0
    for v, r, mult, wd in cuts:
        if (r, tuple(mult)) not in keys:
            pool.append((r, list(mult))); keys.add((r, tuple(mult))); new += 1
    json.dump(pool, open(poolf, 'w'))
    print('  cuts', len(cuts), 'new', new, 'worst', cuts[0][:2] if cuts else None, flush=True)
    if new == 0:
        break
