"""Exact dual certificate: Schrijver SDP (plain, Schrijver's own constraints) + B1 block-type OA rows excludes |C| = M.
usage: make_cert_oa.py n d M outdir [solver=csdp|dd]"""
import sys, os, json, time
from math import comb
from fractions import Fraction
from oa import *
from t3 import partitions, block_wd
import sdpcore, certio

def build(n, d, M):
    opts = dict(lin_basic=True, best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)
    P = split.build(n, d, 0, opts=opts, name='cert_%d_%d_%d' % (n, d, M))
    nbase = len(P.lin)
    meta = [None] * nbase
    seen = set()
    for tot in range(1, n + 1):
        for part in partitions(tot):
            r = len(part)
            if M % 2 ** r == 0: continue
            wd = block_wd(part)
            key = (r, tuple(sorted(wd.items())))
            if key in seen: continue
            seen.add(key)
            if add_oa_row(P, n, M, wd, r):
                meta.append(list(part))
    return P, meta

if __name__ == '__main__':
    n, d, M, outdir = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    solver = sys.argv[5] if len(sys.argv) > 5 else 'csdp'
    K = 60
    os.makedirs(outdir, exist_ok=True)
    P, meta = build(n, d, M)
    work = 'agents/novaB/runs/cert'
    if solver == 'csdp':
        info = sdpcore.solve_csdp(P, work, threads=1, cpus=None)
    else:
        info = P.solve(work, threads=1)
    cert = P.make_certificate(info, K=K)
    b = P.check_certificate(cert)
    print('num', info['bound'], 'rigorous', float(b), 'EXCLUDED' if b < M else 'NOT EXCLUDED', flush=True)
    path = os.path.join(outdir, 'excl_n%d_d%d_M%d_oa.json' % (n, d, M))
    out = certio.export(P, cert, K, path, {'n': n, 'd': d, 'M': M, 'case': 'oa'})
    # attach block lists to the exported lin rows (same order: nonzero z only)
    j = 0
    for z, m in zip(cert['z'], meta):
        if z == 0: continue
        if m: out['lin'][j]['oa_blocks'] = m
        j += 1
    json.dump(out, open(path, 'w'))
    print('wrote', path, 'rows used', len(out['lin']), 'oa rows used', sum(1 for r in out['lin'] if r.get('oa_blocks')))
