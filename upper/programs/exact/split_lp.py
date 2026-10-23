#!/usr/bin/env python3
"""Split Delsarte bound.  Even (N,M,d) code C, shortened at the coordinate/value with the largest
count m: C0 (m words, even weight, length N-1) and C1 (M-m words, odd weight).
Variables: ordered-pair distance distributions a (C0xC0), b (C1xC1), c (C0xC1) on the N-1 coordinates.
Constraints: for every k, [[sum a K_k, sum c K_k],[sum c K_k, sum b K_k]] PSD (Krawtchouk, length N-1);
length-N Delsarte on the whole code; every coordinate has max part <= m (average-distance bound).
Reports, for each m, whether M is feasible."""
import sys
from math import comb
import cvxpy as cp
import numpy as np


def K(k, x, n):
    return sum((-1) ** j * comb(x, j) * comb(n - x, k - j) for j in range(k + 1))


def feasible(N, d, M, m, balance=True, solver='CLARABEL', extra=None):
    n = N - 1
    M1 = M - m
    je = [j for j in range(d, n + 1) if j % 2 == 0]
    jo = [j for j in range(d - 1, n + 1) if j % 2 == 1]
    a = cp.Variable(len(je), nonneg=True)
    b = cp.Variable(len(je), nonneg=True)
    c = cp.Variable(len(jo), nonneg=True)
    cons = [cp.sum(a) == m * m - m, cp.sum(b) == M1 * M1 - M1, cp.sum(c) == m * M1]
    for k in range(0, n + 1):
        Sa = m * K(k, 0, n) + sum(K(k, j, n) * a[i] for i, j in enumerate(je))
        Sb = M1 * K(k, 0, n) + sum(K(k, j, n) * b[i] for i, j in enumerate(je))
        Sc = sum(K(k, j, n) * c[i] for i, j in enumerate(jo))
        cons.append(cp.bmat([[Sa, Sc], [Sc, Sb]]) >> 0)
    # full code, length N: distance D = j (within halves) or j+1 (cross)
    for k in range(1, N + 1):
        tot = M * K(k, 0, N) + sum(K(k, j, N) * (a[i] + b[i]) for i, j in enumerate(je)) \
            + 2 * sum(K(k, j + 1, N) * c[i] for i, j in enumerate(jo))
        cons.append(tot >= 0)
    if balance:
        # sum over the n other coordinates of w(M-w) >= n*m*(M-m); ordered-pair distance sum = 2*sum w(M-w)
        dsum = sum(j * (a[i] + b[i]) for i, j in enumerate(je)) + 2 * sum(j * c[i] for i, j in enumerate(jo))
        cons.append(dsum >= 2 * n * m * (M - m))
    if extra:
        cons += extra(a, b, c, je, jo)
    prob = cp.Problem(cp.Minimize(0), cons)
    try:
        prob.solve(solver=solver)
    except Exception as e:
        return 'error ' + str(e)[:60]
    return prob.status


if __name__ == '__main__':
    N, d, M = map(int, sys.argv[1:4])
    for m in range((M + 1) // 2, M + 1):
        st = feasible(N, d, M, m)
        st2 = feasible(N, d, M, m, balance=False)
        print(f'M={M} m={m}: with balance: {st}; without: {st2}', flush=True)
