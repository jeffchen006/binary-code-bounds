"""P1 audit (round 5): high-precision re-solve of A(n,d), d even, with exact certificates.
V1 Delsarte LP (even weights, HiGHS + exact rational dual check, box A_i <= C(n,i) for residuals)
V2 plain Schrijver 3-point SDP
V3 Schrijver + classical rows (Best18, cw22, MEL19/21, Best20, KT23; premises = Brouwer table 2026-10-21 snapshot)
V4 parity / 2-mod-4 exclusion descent from floor(V3) (odd: S_u odd; 2 mod 4: Kim-Toan/AVZ split over w)
Every SDP value is an exact certificate bound (min over solver attempts: SDPA-DD equilibrated, plain, scaled, GMP).
usage: audit.py n d   -> appends one JSON line to runs/audit/results.jsonl"""
import sys, os, json, time, math
from fractions import Fraction
sys.dont_write_bytecode = True
import numpy as np
import mpmath
import split, parity, sdpcore

BASE = 'agents/sdp'
WORK = os.path.join(BASE, 'runs', 'audit')
os.makedirs(WORK, exist_ok=True)
OFF = dict(best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)


def kraw(n, k, x):
    return sum((-1) ** j * math.comb(x, j) * math.comb(n - x, k - j) for j in range(k + 1))


def delsarte(n, d):
    from scipy.optimize import linprog
    D = [i for i in range(d, n + 1) if i % 2 == 0]
    K = [[kraw(n, k, i) for i in D] for k in range(1, n + 1)]
    # max sum A  s.t.  -sum_i A_i K_k(i) <= C(n,k)
    res = linprog(-np.ones(len(D)), A_ub=-np.array(K, dtype=float), b_ub=[math.comb(n, k) for k in range(1, n + 1)],
                  bounds=[(0, None)] * len(D), method='highs')
    y = [max(Fraction(-float(v)).limit_denominator(10 ** 12), Fraction(0)) for v in res.ineqlin.marginals]
    # |C| <= 1 + sum_k y_k C(n,k) + sum_i C(n,i) max(0, 1 + sum_k y_k K_k(i))
    b = 1 + sum(y[k - 1] * math.comb(n, k) for k in range(1, n + 1))
    for c, i in enumerate(D):
        r = 1 + sum(y[k - 1] * K[k - 1][c] for k in range(1, n + 1))
        if r > 0:
            b += r * math.comb(n, i)
    return 1 - res.fun, b


def solve_exact(P, tag, target):
    """min exact certificate bound over solver attempts; returns (num, rig, how)."""
    best = (None, None, None)
    for how in ('eq', 'plain', 'vs', 'gmp'):
        Q = P
        Q.name = '%s_%s' % (tag, how)
        Q.vscale = {}
        for b in Q.blocks:
            b['scale'] = b.get('scale0', b['scale'])
            b['scale0'] = b['scale']
        old = sdpcore.SDPA
        try:
            if how == 'vs':
                Q.vscale = {v: mpmath.mpf(target) / 2 ** NN[0] for v in Q.vars}
                sdpcore.equilibrate(Q)
            elif how == 'eq':
                sdpcore.equilibrate(Q)
            elif how == 'gmp':
                sdpcore.SDPA = os.path.join(BASE, '..', '..', 'opt', 'sdpa-gmp', 'sdpa_gmp')
                sdpcore.equilibrate(Q)
            info = Q.solve(WORK, param=os.path.join(BASE, 'param.sdpa') if how == 'gmp' else None, threads=1)
            cert = Q.make_certificate(info, K=80)
            rig = Q.check_certificate(cert)
        except Exception as e:
            print(tag, how, 'fail', repr(e)[:80], flush=True)
            continue
        finally:
            sdpcore.SDPA = old
        num = info.get('bound')
        print(tag, how, info['phase.value'], num, float(rig), flush=True)
        if best[1] is None or rig < best[1]:
            best = (num, rig, how)
        # good enough: certificate within 1e-4 of the solver value
        if num is not None and abs(float(rig) - num) < 1e-4 * max(1.0, abs(num)) and info['phase.value'] in ('pdOPT', 'pFEAS', 'pdFEAS'):
            break
        if how == 'plain' and best[1] is not None and abs(float(best[1]) - (best[0] or 0)) < 1e-3:
            break
    return best


NN = [0]


def main(n, d):
    NN[0] = n
    import bounds
    t0 = time.time()
    out = {'n': n, 'd': d}
    try:
        out['table_upper'] = bounds.A_upper(n, d)
    except Exception:
        out['table_upper'] = None
    T = out['table_upper'] or 2 ** n
    lp_num, lp_rig = delsarte(n, d)
    out['delsarte_num'] = lp_num; out['delsarte'] = float(lp_rig)
    P = split.build(n, d, 0, opts=dict(OFF), name='au')
    num, rig, how = solve_exact(P, 'au_%d_%d_schr' % (n, d), T)
    out['schrijver_num'] = num; out['schrijver'] = float(rig) if rig is not None else None; out['schrijver_how'] = how
    P = split.build(n, d, 0, opts={}, name='au')
    num, rig, how = solve_exact(P, 'au_%d_%d_lin' % (n, d), T)
    out['lin_num'] = num; out['lin'] = float(rig) if rig is not None else None; out['lin_how'] = how
    # V4 parity descent from floor(min(V2, V3))
    vals = [v for v in (out['schrijver'], out['lin']) if v is not None]
    F = math.floor(min(vals)) if vals else None
    M = F; steps = []
    while M is not None and len(steps) < 8 and M >= 1:
        cs = parity.cases(M, n)
        if cs[0][0] == '0mod4':
            break
        ok = True
        for case, w in cs:
            Q = parity.build(n, d, 0, M, case, w, opts={})
            _, r, _ = solve_exact(Q, 'au_%d_%d_par%d_%s%s' % (n, d, M, case, w), M)
            if r is None or r >= M:
                ok = False; break
        steps.append([M, ok])
        if not ok:
            break
        M -= 1
    out['floor_sdp'] = F; out['parity_steps'] = steps; out['parity_bound'] = M
    out['time'] = round(time.time() - t0)
    with open(os.path.join(WORK, 'results.jsonl'), 'a') as f:
        f.write(json.dumps(out) + '\n')
    print(json.dumps(out), flush=True)


if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
