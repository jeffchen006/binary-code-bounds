"""usage: certs_from_scan.py scan.txt pool.json outdir  -> make_cert2 for every EXCL line, with the pool prefix used there"""
import sys, subprocess, os
scan, pool, out = sys.argv[1:4]
for ln in open(scan):
    p = ln.split()
    if len(p) < 10 or p[-1] != 'EXCL': continue
    n, d, M, Ap = p[0], p[1], p[2], p[4]
    P = p[p.index('pool') + 1]
    r = subprocess.run(['python3', 'make_cert2.py', n, d, M, Ap, out, pool, P], capture_output=True, text=True)
    print(M, r.stdout.strip().replace('\n', ' | '), r.stderr.strip()[-200:], flush=True)
