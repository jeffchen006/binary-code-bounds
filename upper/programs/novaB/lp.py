"""Delsarte LP (even-weight code, d even) + exact-size rows, float (scipy HiGHS). For exploration only."""
import sys
from math import comb
import numpy as np
from scipy.optimize import linprog

def kraw(n, k, x):
    return sum((-1) ** j * comb(x, j) * comb(n - x, k - j) for j in range(k + 1))

def delsarte(n, d, rows=(), Mfix=None, extra=None):
    """vars A_i, i in D={d..n even}. max 1+sum A_i. rows: list of (coef dict k->w, rhs) meaning
    sum_k w_k * sum_i A_i K_k(i) >= rhs (A_0=1 included)."""
    D = [i for i in range(d, n + 1) if i % 2 == 0]
    K = np.array([[kraw(n, k, i) for i in D] for k in range(n + 1)], float)
    K0 = np.array([kraw(n, k, 0) for k in range(n + 1)], float)
    Aub, bub = [], []
    for k in range(1, n + 1):          # sum A_i K_k(i) >= -K_k(0)
        Aub.append(-K[k]); bub.append(K0[k])
    for w, rhs in rows:
        r = sum(wk * K[k] for k, wk in w.items()); c = sum(wk * K0[k] for k, wk in w.items())
        Aub.append(-r); bub.append(c - rhs)
    if extra:
        for a, b in extra(D):
            Aub.append(a); bub.append(b)
    res = linprog(-np.ones(len(D)), A_ub=np.array(Aub), b_ub=np.array(bub), bounds=[(0, None)] * len(D), method='highs')
    if res.status != 0:
        return None, None
    A = dict(zip(D, res.x))
    return 1 + res.x.sum(), A

def betas(n, M, A):
    return {k: M * (kraw(n, k, 0) + sum(a * kraw(n, k, i) for i, a in A.items())) / comb(n, k) for k in range(n + 1)}

if __name__ == '__main__':
    n, d = int(sys.argv[1]), int(sys.argv[2])
    v, A = delsarte(n, d)
    print('LP', v)
    b = betas(n, v, A)
    print({k: round(x, 3) for k, x in b.items()})
