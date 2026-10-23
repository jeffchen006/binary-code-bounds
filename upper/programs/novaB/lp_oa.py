import sys
from math import comb
from lp import *
def coord_rows(n, M, rmax=None):
    rows = []
    for r in range(1, (rmax or n) + 1):
        rho = M % 2 ** r
        if rho:
            rows.append(({k: comb(r, k) / comb(n, k) for k in range(1, r + 1)}, rho * (2 ** r - rho) / M))
    return rows
n, d = int(sys.argv[1]), int(sys.argv[2])
for M in map(int, sys.argv[3:]):
    v0, _ = delsarte(n, d)
    v, A = delsarte(n, d, coord_rows(n, M))
    b = betas(n, M, A) if A else {}
    print(M, 'M mod 64 =', M % 64, 'LP', round(v0, 3), 'LP+OAcoord', v and round(v, 3), {k: round(x, 2) for k, x in b.items() if 0 < k < n and abs(x) > 1e-6})
