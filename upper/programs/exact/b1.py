#!/usr/bin/env python3
"""B1 coset-count rows (novaB, METHOD_B1.md) on sdp's GMS pipeline (gms0.build_full + homogenize), read-only import.
Row for a linear V <= F2^n, dim r, nonzero weight counts m_k:
    sum_k m_k/C(n,k) * sum_i A_i K_k(i) >= rho(2^r - rho)/T,   rho = T mod 2^r,  A_0 = 1, sum A_i = |C|.
Valid for every code of size exactly T.  Objective of the homogenized SDP = max |C|;  max < T excludes T.

Functions:
  build_gp(n, d, T, case, w, lin, light)  -> P exactly as agents/sdp/gmspar.py builds it (before homogenize)
  X0forms(P, n)                           -> {a: {var: coef}} with A_a = sum coef*X(var)
  add_b1(P, n, T, mk, r)                  -> append one row (normalized form, before homogenize)
  betas_from_out(P, n, T, outfile)        -> distance distribution A and shell averages beta_k at an SDPA solution
"""
import sys, os
SDP = 'agents/sdp'
NOVAB = 'agents/novaB'
sys.path.insert(0, SDP)
pass
from fractions import Fraction
import mpmath
import gms0, parity, sdpcore
from sdpcore import binom

OFF = dict(best18=False, cw22=False, mel19=False, mel21=False, best20=False, kt23=False)


def kraw(n, k, x):
    return sum((-1) ** j * binom(x, j) * binom(n - x, k - j) for j in range(k + 1))


def build_gp(n, d, T, case, w, lin=True, light=True, name='b1'):
    opts = {} if lin else OFF
    P = gms0.build_full(n, d, 1, opts=opts, name=name, light=light, quad_lin=True)
    if case != 'base':
        Q = parity.build(n, d, 0, T, case, w, opts=dict(OFF, lin_basic=False))
        P.lin = P.lin + Q.lin
    return P


def X0forms(P, n):
    X0 = {}
    for v in P.vars:
        try:
            ok = isinstance(v, tuple) and len(v) == 3 and isinstance(v[0], tuple) and all(s == 0 for s in v[0]) and v[1] == v[2]
        except TypeError:
            ok = False
        if ok:
            a = sum(v[1])
            X0.setdefault(a, {})[v] = Fraction(binom(n, a))
    return X0


_SF = {}


def shell_forms(P, n):
    if _SF.get('P') is not P:
        X0 = X0forms(P, n)
        _SF.clear(); _SF['P'] = P
        for k in range(n + 1):
            co = {}
            for a, f in X0.items():
                K = kraw(n, k, a)
                for v, c in f.items():
                    co[v] = co.get(v, 0) + c * K
            _SF[k] = (Fraction(kraw(n, k, 0)), {v: c for v, c in co.items() if c})
    return _SF


def add_b1(P, n, T, mk, r):
    rho = T % (2 ** r)
    rhs = Fraction(rho * (2 ** r - rho), T)
    if rhs == 0:
        return False
    SF = shell_forms(P, n)
    g0 = -rhs; co = {}
    for k, m in mk.items():
        c0, f = SF[k]
        wgt = Fraction(m, binom(n, k))
        g0 += wgt * c0
        for v, c in f.items():
            co[v] = co.get(v, 0) + wgt * c
    co = {v: c for v, c in co.items() if c}
    mx = max([abs(c) for c in co.values()] + [abs(g0)])
    P.lin.append((g0 / mx, {v: c / mx for v, c in co.items()}))
    return True


def betas_from_out(P, n, T, outfile):
    """P must be the homogenized problem with vscale set exactly as in the run. X = x/p."""
    info = P.read_result(outfile)
    vs = getattr(P, 'vscale', {})
    val = {v: info['x'][i] * vs.get(v, 1) for i, v in enumerate(P.vars)}
    p = val[('p',)]
    X0 = X0forms(P, n)
    A = {0: mpmath.mpf(1)}
    for a, f in X0.items():
        A[a] = sum(mpmath.mpf(c.numerator) / c.denominator * val[v] for v, c in f.items()) / p
    beta = {k: T * sum(A[a] * kraw(n, k, a) for a in A) / binom(n, k) for k in range(n + 1)}
    return info, p, A, beta


def homog(P, n, T):
    P = gms0.homogenize(P, n, xle_p=True)
    P.vscale = {v: mpmath.mpf(T) / 2 ** n for v in P.vars}
    sdpcore.equilibrate(P)
    return P
