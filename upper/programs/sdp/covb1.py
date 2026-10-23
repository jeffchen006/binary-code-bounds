"""P4: odd-code Schrijver SDP for A(N, 3) = A(N+1, 4) + extended-code rows (Best/cw/MEL) + Johnson-type covering rows
(covering.py, parameter T = M, valid for every code with |C| >= M) + B1 coset-count rows at |C| = M exactly + cut loop.
Excluded when the exact certificate bound is < M. Premises as cplin.py (team-audited values; TLY 6551 replaced).
usage: covb1.py N M rounds"""
import sys, os, json, random, time, pickle
pass
sys.dont_write_bytecode = True
sys.path.insert(1, 'agents/novaB')
from fractions import Fraction
import bounds
bounds.OVERRIDE.update({(17, 4): 3272, (18, 4): 6543, (19, 4): 13085, (21, 4): 43681, (25, 4): 599174,
                        (26, 4): 1198348, (27, 4): 2396695})
if os.environ.get('CWLIVE') == '1':
    import cwlive
    cwlive.apply()
import covering, sdpcore
from oa import add_oa_row, betas
from t3 import block_rows
from cp import separate, wd_of

N, M, rounds = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
random.seed(1)
work = 'agents/sdp/runs/covb1'
os.makedirs(work, exist_ok=True)
tag = 'cb_%d_%d' % (N, M)
poolf = os.path.join(work, tag + '.pool.json')
pool = json.load(open(poolf)) if os.path.exists(poolf) else []
for it in range(rounds):
    t = time.time()
    P = covering.build(N, 3, T=M, name=tag + '_r%d' % it)
    nb = 0
    if os.environ.get('MEL07') == '1':
        # Mounits-Etzion-Litsyn 2007, Corollary 2 (n even, d = 3): |C| (2n + 3 - S1/n) <= 2^n,
        # S1 = 6(A_{n-3} + A_{n-2}) + 3n(A_{n-1} + A_n). At |C| = M: S1 >= n(2n + 3) - n 2^n / M.
        assert N % 2 == 0
        from oa import X0forms
        X0 = X0forms(P, N)
        co = {}
        for a, c in ((N - 3, 6), (N - 2, 6), (N - 1, 3 * N), (N, 3 * N)):
            for v, f in X0.get(a, {}).items():
                co[v] = co.get(v, 0) + c * f
        g0 = -Fraction(N * (2 * N + 3)) + Fraction(N * 2 ** N, M)
        mx = max([abs(c) for c in co.values()] + [abs(g0)])
        P.lin.append((g0 / mx, {v: c / mx for v, c in co.items()}))
    for r, wd in block_rows(N, M, rmin=1):
        nb += add_oa_row(P, N, M, wd, r)
    for r, mult in pool:
        wd = wd_of(r, mult)
        if wd is not None:
            nb += add_oa_row(P, N, M, wd, r)
    if os.environ.get('NOEQ') != '1':
        sdpcore.equilibrate(P)
    info = P.solve(work, threads=1)
    try:
        cert = P.make_certificate(info, K=80); rig = P.check_certificate(cert)
        pickle.dump({'cert': cert, 'pool': pool}, open(os.path.join(work, P.name + '.cert.pkl'), 'wb'))
    except Exception as e:
        rig = None; print('certfail', repr(e)[:80])
    print('round', it, 'lin', len(P.lin), 'B1', nb, info['phase.value'], 'num', info['bound'],
          'exact %.6f' % float(rig) if rig is not None else 'exact NA', 'EXCL' if (rig is not None and rig < M) else 'feas',
          '%.0fs' % (time.time() - t), flush=True)
    if rig is not None and rig < M:
        break
    A, beta = betas(P, N, M, [float(v) for v in info['x']])
    cuts = separate(beta, N, M)
    keys = {(r, tuple(m)) for r, m in pool}; new = 0
    for v, r, mult, wd in cuts:
        if (r, tuple(mult)) not in keys:
            pool.append((r, list(mult))); keys.add((r, tuple(mult))); new += 1
    json.dump(pool, open(poolf, 'w'))
    print('  cuts', len(cuts), 'new', new, flush=True)
    if new == 0:
        break
