"""d8 request: is the forced DSL distance distribution (n=28, d=4, M=4718592) compatible with (a) B1 rows (exact check over
all block subspaces with rho != 0), (b) the Schrijver 3-point SDP (caps x(i,0,0) <= A_i/C(n,i), max |C|; exact certificate)?"""
import sys, os, math
sys.dont_write_bytecode = True
from fractions import Fraction as F
sys.path.insert(1, 'agents/novaB')
n, d = 28, 4
A = {0: 1, 4: 805, 6: 12908, 8: 110187, 10: 459536, 12: 1072239, 14: 1407240, 16: 1072239, 18: 459536, 20: 110187,
     22: 12908, 24: 805, 28: 1}
M = sum(A.values()); print('M', M, M == 4718592)
def K(k, x): return sum((-1) ** j * math.comb(x, j) * math.comb(n - x, k - j) for j in range(k + 1))
T = {k: sum(A[i] * K(k, i) for i in A) for k in range(n + 1)}        # sum_{|u|=k} S_u^2 / M
print('nonzero shells', {k: T[k] for k in T if T[k]})
assert all(T[k] >= 0 for k in T)
from t3 import partitions, block_wd
worst = None; nrows = 0
for tot in range(1, n + 1):
    for part in partitions(tot):
        r = len(part); rho = M % 2 ** r
        if rho == 0: continue
        wd = block_wd(part); nrows += 1
        lhs = sum(F(m, math.comb(n, k)) * T[k] for k, m in wd.items())
        rhs = F(rho * (2 ** r - rho), M)
        sl = lhs - rhs
        if worst is None or sl < worst[0]: worst = (sl, part)
print('B1 block rows', nrows, 'min slack', float(worst[0]), 'at', worst[1], 'VIOLATED' if worst[0] < 0 else 'all satisfied')
if '--sdp' in sys.argv:
    pass
    import split, sdpcore, mpmath
    from oa import X0forms
    P = split.build(n, d, 0, opts=dict(best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False), name='dsl28')
    X0 = X0forms(P, n)
    for a, f in X0.items():
        for v, c in f.items():       # c = C(n,a): A_a = c * x  <= A[a]
            P.lin.append((F(A.get(a, 0)), {v: -c}))
    best = None
    for how in ('eq', 'vs', 'plain'):
        P.vscale = {}
        for b in P.blocks:
            b['scale'] = b.get('scale0', b['scale']); b['scale0'] = b['scale']
        if how == 'vs':
            P.vscale = {v: mpmath.mpf(M) / 2 ** n for v in P.vars}
        if how != 'plain':
            sdpcore.equilibrate(P)
        P.name = 'dsl28_' + how
        out = 'agents/sdp/runs/dsl28/%s.out' % P.name
        if '--read' in sys.argv and os.path.exists(out) and os.path.getsize(out) > 0:
            info = P.read_result(out)
            info['bound'] = float(P.c0) - float(info['objValDual'])
        else:
            info = P.solve('agents/sdp/runs/dsl28', threads=1)
        try:
            cert = P.make_certificate(info, K=80); rig = P.check_certificate(cert)
        except Exception as e:
            print(how, 'cert fail', e); continue
        print(how, info['phase.value'], 'num', info['bound'], 'exact', float(rig), 'M - exact = %.6f' % (M - float(rig)), flush=True)
        if best is None or rig < best: best = rig
    print('RESULT: exact max |C| =', float(best), 'EXCLUDED (infeasible)' if best < M else 'NOT excluded (Schrijver-feasible)')
