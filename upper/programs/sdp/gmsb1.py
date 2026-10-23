"""GMS 4-point SDP (homogenized, max 2^n p) + parity rows + novaB's B1 coset-count rows, cutting-plane loop.
B1 row (V linear, dim r, m_k nonzero words of weight k, rho = M mod 2^r):
    sum_k m_k/C(n,k) * sum_i A_i K_k(i) >= rho (2^r - rho) / M
Rounds: 0 = block subspaces (all partitions, rho != 0); later rounds add general V found by novaB's
separation (cp.separate) at the current optimum. Excluded when max 2^n p < M (rigorous: certificate bound).
usage: gmsb1.py n d M case w [--rounds R] [--lin] [--rmax R] [--tag T]   (case: odd|2mod4|0mod4|none)"""
import sys, os, time, json, pickle, argparse, random
pass
sys.dont_write_bytecode = True
sys.path.insert(1, 'agents/novaB')
from fractions import Fraction
import mpmath
import gms0, parity, sdpcore, split
from oa import add_oa_row, betas
from t3 import block_rows
from cp import separate, wd_of

ap = argparse.ArgumentParser()
ap.add_argument('n', type=int); ap.add_argument('d', type=int); ap.add_argument('M', type=int)
ap.add_argument('case'); ap.add_argument('w')
ap.add_argument('--rounds', type=int, default=3); ap.add_argument('--lin', action='store_true')
ap.add_argument('--rmax', type=int, default=10); ap.add_argument('--tag', default='')
ap.add_argument('--threads', type=int, default=4); ap.add_argument('--noblock', action='store_true')
a = ap.parse_args()
n, d, M = a.n, a.d, a.M
w = None if a.w == 'None' else int(a.w)
off = dict(best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)
OPTS = {} if a.lin else off
base = 'gb_%d_%d_%d_%s_%s_l%d%s' % (n, d, M, a.case, w, int(a.lin), a.tag)
work = 'agents/sdp/runs/gmsb1'
os.makedirs(work, exist_ok=True)
poolf = os.path.join(work, base + '.pool.json')
pool = json.load(open(poolf)) if os.path.exists(poolf) else []
os.environ['RMAX'] = str(min(a.rmax, 10))
random.seed(12345)


Qref = {}


def build(pool):
    t = time.time()
    P = gms0.build_full(n, d, 1, opts=OPTS, name=base, light=True, quad_lin=True)
    Q = split.build(n, d, 0, opts=dict(off, lin_basic=False), name='q')
    nb = 0
    if a.case in ('odd', '2mod4', '0mod4'):
        Q2 = parity.build(n, d, 0, M, a.case, w, opts=dict(off, lin_basic=False))
        Q.lin += Q2.lin
    if not a.noblock:
        for r, wd in block_rows(n, M, rmin=1 if a.case != 'odd' else 2):
            if r <= a.rmax:
                nb += add_oa_row(Q, n, M, wd, r)
    for r, mult in pool:
        wd = wd_of(r, mult)
        if wd is not None:
            nb += add_oa_row(Q, n, M, wd, r)
    Qref['Q'] = Q
    for g0, co in Q.lin:
        assert all(v in P.vidx for v in co), 'var mismatch'
    P.lin = P.lin + Q.lin
    P = gms0.homogenize(P, n, xle_p=True)
    P.vscale = {v: mpmath.mpf(M) / 2 ** n for v in P.vars}
    sdpcore.equilibrate(P)
    print('build vars %d blocks %d lin %d B1rows %d pool %d %.0fs' % (len(P.vars), len(P.blocks), len(P.lin), nb, len(pool), time.time() - t), flush=True)
    return P


for it in range(a.rounds):
    t = time.time()
    P = build(pool)
    P.name = base + '_r%d' % it
    info = P.solve(work, threads=a.threads)
    res = '%s %s num=%s pnum=%s' % (P.name, info['phase.value'], info['bound'], info.get('pbound'))
    rig = None
    try:
        cert = P.make_certificate(info, K=80); rig = P.check_certificate(cert)
        pickle.dump({'cert': cert, 'bound': rig, 'pool': pool}, open(os.path.join(work, P.name + '.cert.pkl'), 'wb'))
        res += ' rig=%.6f %s' % (float(rig), 'EXCL' if rig < M else '')
    except Exception as e:
        res += ' certfail %r' % e
    res += ' %.0fs' % (time.time() - t)
    print(res, flush=True)
    open('agents/sdp/gms_results.txt', 'a').write(time.strftime('%m-%d %H:%M ') + res + '\n')
    if rig is not None and rig < M:
        break
    # primal point -> distance distribution (X = x/p; SDPA variable y with x = vscale*y)
    y = dict(zip(P.vars, info['x']))
    p = y[('p',)]
    X = {v: float(y[v] / p) for v in P.vars if v != ('p',)}
    A, beta = betas(Qref['Q'], n, M, X)
    print('A', {i: round(v, 4) for i, v in sorted(A.items()) if abs(v) > 1e-6}, flush=True)
    print('beta', {k: round(v, 4) for k, v in sorted(beta.items())}, flush=True)
    cuts = separate(beta, n, M, rs=range(2, min(a.rmax, 10) + 1))
    keys = {(r, tuple(m)) for r, m in pool}; new = 0
    for v, r, mult, wd in cuts:
        if (r, tuple(mult)) in keys:
            continue
        pool.append((r, list(mult))); keys.add((r, tuple(mult))); new += 1
    json.dump(pool, open(poolf, 'w'))
    print('separation: violated %d new %d pool %d worst %s' % (len(cuts), new, len(pool), cuts[0][:2] if cuts else None), flush=True)
    if new == 0:
        break
