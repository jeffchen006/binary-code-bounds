#!/usr/bin/env python3
"""Rebuild a GMS SDP exactly as gms.py did, read the SDPA result, make a rational certificate, check it exactly.
usage: gms_cert.py n d ub outfile [--parity M | --mod4 M w] [--extra list]
Plain: prints exact upper bound on |C|.  Parity/mod4: exact upper bound on s (s < 1 excludes the case)."""
import sys, pickle, time, argparse
from fractions import Fraction
import gms
ap = argparse.ArgumentParser()
ap.add_argument('n', type=int); ap.add_argument('d', type=int); ap.add_argument('ub', type=int); ap.add_argument('out')
ap.add_argument('--parity', type=int); ap.add_argument('--mod4', type=int, nargs=2); ap.add_argument('--extra', default='')
ap.add_argument('--K', type=int, default=200)
a = ap.parse_args()
t = time.time()
P = gms.build(a.n, a.d, Fraction(a.ub), verbose=False)
if a.extra:
    gms.add_extra(P, a.n, a.d, {k: True for k in a.extra.split(',')}, verbose=False)
if a.parity:
    gms.add_parity(P, a.n, a.d, a.parity)
if a.mod4:
    gms.add_mod4(P, a.n, a.d, a.mod4[0], a.mod4[1])
print('built', len(P.vars), len(P.blocks), len(P.lin), '%.0fs' % (time.time() - t), flush=True)
info = P.read_result(a.out)
cert = P.make_certificate(info, K=a.K)
pickle.dump(cert, open(a.out + '.cert.pkl', 'wb'))
b = P.check_certificate(cert)
what = 's' if (a.parity or a.mod4) else '|C|'
print('EXACT bound on %s: %.9f' % (what, float(b)), ('EXCLUDED' if (what == 's' and b < 1) else ''), '%.0fs' % (time.time() - t), flush=True)
