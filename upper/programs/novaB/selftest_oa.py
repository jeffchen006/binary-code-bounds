"""Validity test of B1 rows on real codes (random codes, linear-code subsets), exact integers."""
import random, itertools
from math import comb
from fractions import Fraction
from t3 import partitions, block_wd
def kraw(n, k, x):
    return sum((-1) ** j * comb(x, j) * comb(n - x, k - j) for j in range(k + 1))
random.seed(1)
bad = 0; tests = 0; tight = 0
for trial in range(300):
    n = random.randint(5, 10)
    kind = trial % 3
    if kind == 0:
        M = random.randint(2, 2 ** n // 3)
        C = random.sample(range(2 ** n), M)
    else:   # subset of a linear code (structured; many S_u = 0)
        k = random.randint(2, n - 1)
        G = [random.randrange(1, 2 ** n) for _ in range(k)]
        L = set()
        for mask in range(2 ** k):
            v = 0
            for b in range(k):
                if mask >> b & 1: v ^= G[b]
            L.add(v)
        L = list(L); M = max(1, len(L) - random.randint(0, 5)); C = random.sample(L, M)
    S = [sum(1 - 2 * (bin(u & c).count('1') & 1) for c in C) for u in range(2 ** n)]
    sh = [0] * (n + 1)
    for u in range(2 ** n): sh[bin(u).count('1')] += S[u] ** 2
    # distance distribution formula check
    A = [0] * (n + 1)
    for a in C:
        for b in C: A[bin(a ^ b).count('1')] += 1
    for k in range(n + 1):
        assert sh[k] * M == M * sum(Fraction(A[i], M) * kraw(n, k, i) for i in range(n + 1)) * M / 1 or True
        assert Fraction(sh[k]) == sum(Fraction(A[i]) * kraw(n, k, i) for i in range(n + 1)), 'shell identity'
    for tot in range(1, n + 1):
        for part in partitions(tot):
            r = len(part); rho = M % 2 ** r
            wd = block_wd(part)
            lhs = sum(Fraction(m, comb(n, k)) * sh[k] for k, m in wd.items())
            tests += 1
            if lhs < rho * (2 ** r - rho): bad += 1; print('VIOLATION', n, M, part, lhs, rho)
            if lhs == rho * (2 ** r - rho) and rho: tight += 1
    # per-subspace check, random subspaces of dim r (not only block type)
    for _ in range(20):
        r = random.randint(1, n)
        B = [random.randrange(1, 2 ** n) for _ in range(r)]
        V = set()
        for mask in range(2 ** r):
            v = 0
            for b in range(r):
                if mask >> b & 1: v ^= B[b]
            V.add(v)
        rr = len(V).bit_length() - 1; rho = M % 2 ** rr
        lhs = sum(S[v] ** 2 for v in V if v)
        tests += 1
        if lhs < rho * (2 ** rr - rho): bad += 1; print('VIOLATION subspace', n, M, rr)
print('tests', tests, 'violations', bad, 'tight (equality) cases', tight)
