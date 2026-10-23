"""Schrijver+parity certificates for given (n,d,M,case,w-list). usage: make_cases.py n d M case outdir w1 w2 ..."""
import sys, os
import parity, certio
n, d, M, case, outdir = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
ws = [None if x == 'None' else int(x) for x in sys.argv[6:]]
K = 60
opts = dict(lin_basic=True, best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)
os.makedirs(outdir, exist_ok=True)
for w in ws:
    Q = parity.build(n, d, 0, M, case, w, opts=opts, name='rc_%d_%d_%d_%s_%s' % (n, d, M, case, w))
    info = Q.solve('agents/sdp/runs/rc', threads=1)
    cert = Q.make_certificate(info, K=K); b = Q.check_certificate(cert)
    fn = 'excl_n%d_d%d_M%d_%s_w%s.json' % (n, d, M, case, w)
    certio.export(Q, cert, K, os.path.join(outdir, fn), {'n': n, 'd': d, 'M': M, 'case': case, 'w': w})
    print(fn, info['phase.value'], float(b), 'EXCL' if b < M else 'FAIL', flush=True)
