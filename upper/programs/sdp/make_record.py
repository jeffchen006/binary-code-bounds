"""Produce JSON dual certificates proving A(n,d) <= target_final via Schrijver SDP + parity exclusions.
usage: make_record.py n d outdir M_start M_end   (excludes all M in [M_end+1 .. M_start], base shows A <= M_start)"""
import sys, os, json, time
import split, parity, certio
n, d, outdir, Ms, Me = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
K = 60
opts = dict(lin_basic=True, best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)
os.makedirs(outdir, exist_ok=True)
work = 'agents/sdp/runs/record_%d_%d' % (n, d)
def do(P, problem, fname):
    info = P.solve(work, threads=4)
    cert = P.make_certificate(info, K=K)
    b = P.check_certificate(cert)
    certio.export(P, cert, K, os.path.join(outdir, fname), problem)
    print(fname, info['phase.value'], 'num', info['bound'], 'rig', float(b), flush=True)
    return b
P = split.build(n, d, 0, opts=opts, name='rec_base_%d_%d' % (n, d))
b = do(P, {'n': n, 'd': d, 'M': None, 'target': Ms + 1}, 'base_n%d_d%d.json' % (n, d))
assert b < Ms + 1
for M in range(Ms, Me, -1):
    for case, w in parity.cases(M, n):
        assert case in ('odd', '2mod4'), 'M = 0 mod 4 cannot be excluded'
        Q = parity.build(n, d, 0, M, case, w, opts=opts, name='rec_%d_%d_%d_%s_%s' % (n, d, M, case, w))
        b = do(Q, {'n': n, 'd': d, 'M': M, 'case': case, 'w': w}, 'excl_n%d_d%d_M%d_%s_w%s.json' % (n, d, M, case, w))
        assert b < M, ('NOT excluded', M, case, w, float(b))
print('DONE: A(%d,%d) <= %d' % (n, d, Me))
