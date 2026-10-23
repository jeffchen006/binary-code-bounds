#!/usr/bin/env python3
"""Delsarte LP for an even (n,M,d) code: range of sum_j j*A_j (= (2/M) * sum over coords of w(M-w)).
If max over the LP region < 2*n*m*(M-m)/M, then a code with all coordinate parts <= m cannot exist."""
import sys
from math import comb
from scipy.optimize import linprog
from fractions import Fraction

def K(k, x, n):
    return sum((-1)**j * comb(x, j) * comb(n - x, k - j) for j in range(k + 1))

def lp(n, d, M, sense):
    js = [j for j in range(d, n + 1) if j % 2 == 0]
    # variables A_j, j in js
    A_eq = [[1] * len(js)]; b_eq = [M - 1]
    A_ub = []; b_ub = []
    for k in range(1, n + 1):
        A_ub.append([-K(k, j, n) for j in js]); b_ub.append(K(k, 0, n))
    c = [sense * j for j in js]
    r = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * len(js), method='highs')
    return r, js

n, d, M = map(int, sys.argv[1:4])
for sense, name in ((1, 'min'), (-1, 'max')):
    r, js = lp(n, d, M, sense)
    if r.status != 0:
        print(name, 'infeasible', r.message); continue
    val = sense * r.fun
    print(name, 'sum j A_j =', round(val, 4), ' dist:', {j: round(a, 3) for j, a in zip(js, r.x) if a > 1e-9})
    if name == 'max':
        for m in range((M + 1) // 2, M):
            f = 2 * n * m * (M - m) / M
            if f > val + 1e-9:
                print(f'  m={m}: need sum jA_j >= {f:.3f} > LPmax -> impossible')
