"""Validity of general-V B1 rows (aggregated over S_n, via shells) and of the cap row, on real codes. Exact integers."""
import random
from math import comb
from fractions import Fraction
random.seed(2)
def wd_of(r, mult):
    wd = {}
    for x in range(1, 2 ** r):
        w = sum(m for c, m in enumerate(mult, 1) if m and bin(x & c).count('1') & 1)
        if w == 0: return None
        wd[w] = wd.get(w, 0) + 1
    return wd
def best_code(n, d, tries=30):   # random greedy code with min distance d (to test the cap row with Ap = true max of shortened sizes)
    best = []
    for _ in range(tries):
        pts = list(range(2 ** n)); random.shuffle(pts); C = []
        for p in pts:
            if all(bin(p ^ c).count('1') >= d for c in C): C.append(p)
        if len(C) > len(best): best = C
    return best
tests = bad = tight = 0
for trial in range(400):
    n = random.randint(5, 10)
    if trial % 2:
        M = random.randint(2, 2 ** n // 3); C = random.sample(range(2 ** n), M)
    else:
        C = best_code(n, random.choice([2, 3, 4]), 3); C = random.sample(C, max(1, len(C) - random.randint(0, 3))); M = len(C)
    S = [sum(1 - 2 * (bin(u & c).count('1') & 1) for c in C) for u in range(2 ** n)]
    sh = [0] * (n + 1)
    for u in range(2 ** n): sh[bin(u).count('1')] += S[u] ** 2
    for _ in range(60):
        r = random.randint(1, min(n, 6)); T = 2 ** r - 1
        mult = [0] * T
        for b in range(r): mult[(1 << b) - 1] = 1
        for _ in range(random.randint(0, n - r)): mult[random.randrange(T)] += 1
        wd = wd_of(r, mult)
        if wd is None: continue
        rho = M % 2 ** r
        lhs = sum(Fraction(m, comb(n, k)) * sh[k] for k, m in wd.items())
        tests += 1
        if lhs < rho * (2 ** r - rho): bad += 1; print('VIOL', n, M, r, mult)
        if rho and lhs == rho * (2 ** r - rho): tight += 1
    # cap row with Ap = max shortened size (true value for this code)
    Ap = max(max(sum(1 for c in C if (c >> i) & 1 == b) for b in (0, 1)) for i in range(n))
    tests += 1
    if sh[1] > n * (2 * Ap - M) ** 2: bad += 1; print('CAP VIOL', n, M, Ap)
print('tests', tests, 'violations', bad, 'tight', tight)
