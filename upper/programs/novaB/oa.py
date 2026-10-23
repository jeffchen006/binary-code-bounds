"""B1: orthogonal-array integrality rows on top of the sdp agent's Schrijver SDP (split.py, n1=0).
Row for a linear code V <= F2^n of dim r with weight counts m_k (nonzero words):
   sum_k m_k/C(n,k) * sum_i A_i K_k(i) >= rho(2^r - rho)/M ,  rho = M mod 2^r      (A_0 = 1, sum A_i = M)
Valid for every code of size exactly M (coset counts of V^perp are integers)."""
import sys, os
sys.path.insert(0, 'agents/sdp')
pass
from fractions import Fraction
import split
from sdpcore import binom

def kraw(n, k, x):
    return sum((-1) ** j * binom(x, j) * binom(n - x, k - j) for j in range(k + 1))

def X0forms(P, n):
    X0 = {}
    for v in P.vars:
        if all(s == 0 for s in v[0]) and v[1] == v[2]:
            a = sum(v[1])
            X0.setdefault(a, {})[v] = Fraction(binom(n, a))
    return X0

def shell_form(P, n, k, X0=None):
    """returns (const, {var:coef}) of  sum_i A_i K_k(i)"""
    X0 = X0 or X0forms(P, n)
    co = {}
    for a, f in X0.items():
        K = kraw(n, k, a)
        for v, c in f.items():
            co[v] = co.get(v, 0) + c * K
    return Fraction(kraw(n, k, 0)), {v: c for v, c in co.items() if c}

_SF = {}
def add_oa_row(P, n, M, mk, r):
    """mk: dict k->count of nonzero words of weight k in V (dim r)."""
    rho = M % (2 ** r)
    rhs = Fraction(rho * (2 ** r - rho), M)
    if rhs == 0:
        return False
    if _SF.get('P') is not P:
        X0 = X0forms(P, n)
        _SF.clear(); _SF['P'] = P
        for k in range(n + 1):
            _SF[k] = shell_form(P, n, k, X0)
    g0 = -rhs; co = {}
    for k, m in mk.items():
        c0, f = _SF[k]
        w = Fraction(m, binom(n, k))
        g0 += w * c0
        for v, c in f.items():
            co[v] = co.get(v, 0) + w * c
    co = {v: c for v, c in co.items() if c}
    mx = max([abs(c) for c in co.values()] + [abs(g0)])
    P.lin.append((g0 / mx, {v: c / mx for v, c in co.items()}))
    return True

def betas(P, n, M, x):
    """avg S_u^2 over shell k at a primal point x (list aligned with P.vars order)"""
    X0 = X0forms(P, n)
    vals = dict(zip(P.vars, x)) if not isinstance(x, dict) else x
    A = {0: 1.0}
    for a, f in X0.items():
        A[a] = sum(float(c) * float(vals[v]) for v, c in f.items())
    out = {}
    for k in range(n + 1):
        out[k] = M * sum(A[a] * kraw(n, k, a) for a in A) / binom(n, k)
    return A, out
