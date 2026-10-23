#!/usr/bin/env python3
"""Read a finished gmspar-type SDPA solution, print distance distribution / shell averages, and search
B1 rows violated at that point (block V's exhaustively + novaB's local search over general V).
usage: sep_out.py n d T case w outfile [lin] [light]"""
import sys, os, random, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b1
sys.path.insert(0, b1.NOVAB)
import cp, t3
n, d, T = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
case = sys.argv[4]; w = None if sys.argv[5] == 'None' else int(sys.argv[5])
out = sys.argv[6]
lin = len(sys.argv) <= 7 or sys.argv[7] == '1'
light = len(sys.argv) <= 8 or sys.argv[8] == '1'
t = time.time()
P = b1.build_gp(n, d, T, case, w, lin, light)
P = b1.homog(P, n, T)
print('built', len(P.vars), '%.0fs' % (time.time() - t), flush=True)
info, p, A, beta = b1.betas_from_out(P, n, T, out)
print('phase', info['phase.value'], 'p*2^n =', float(p * 2 ** n), 'sum A =', float(sum(A.values())))
print('A:', {a: round(float(v), 5) for a, v in sorted(A.items())})
bf = {k: float(v) for k, v in beta.items()}
print('beta:', {k: round(v, 4) for k, v in bf.items()})
viol = []
for r, wd in t3.block_rows(n, T, rmin=1):
    rho = T % 2 ** r
    lhs = sum(m * bf[k] for k, m in wd.items())
    viol.append((lhs - rho * (2 ** r - rho), r, wd))
viol.sort(key=lambda z: z[0])
print('block rows: %d, most violated:' % len(viol))
for z in viol[:10]:
    print('  %.4f r=%d %s' % z)
random.seed(1)
cuts = cp.separate(bf, n, T, rs=range(2, 9), tries=60, steps=400)
print('general V cuts:', len(cuts))
for c in cuts[:10]:
    print('  %.4f r=%d mult=%s wd=%s' % (c[0], c[1], c[2], c[3]))
json.dump({'A': {a: float(v) for a, v in A.items()}, 'beta': bf, 'cuts': [(c[0], c[1], c[2]) for c in cuts[:200]],
           'blockviol': [(z[0], z[1], sorted(z[2].items())) for z in viol[:200]]},
          open(out + '.sep.json', 'w') if os.access(os.path.dirname(out), os.W_OK) and 'exact' in out else
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', os.path.basename(out) + '.sep.json'), 'w'))
print('done %.0fs' % (time.time() - t))
