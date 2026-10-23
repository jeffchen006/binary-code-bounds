"""Valid upper bounds: A(n,d) (Brouwer table), A(n,d,w) (Agrell cw table + trivial rules),
T(w1,n1,w2,n2,d) (Agrell dcw tables + Johnson-type recursion)."""
import re, os, math
from functools import lru_cache
import htmltab
H = os.path.dirname(os.path.abspath(__file__))

def _num(s):
    s = s.replace('\n', ' ').strip()
    m = re.findall(r'\d+', s)
    return [int(x) for x in m]

# ---- A(n,d) upper bounds from Brouwer table (ref/binary-1.txt) ----
def _load_brouwer():
    A = {}
    lines = open(os.path.join(H, '..', '..', 'ref', 'binary-1.txt'), encoding='latin-1').read().splitlines()
    for ln in lines:
        p = ln.split()
        if len(p) == 8 and p[0].isdigit() and 6 <= int(p[0]) <= 28:
            n = int(p[0])
            for k, tok in enumerate(p[1:]):
                d = 4 + 2 * k
                v = tok.split('-')[-1]
                A[(n, d)] = int(v)
    return A
_AB = _load_brouwer()
OVERRIDE = {}   # proven improvements found by us: (n,d)->ub

def A_upper(n, d):
    if d % 2:
        return A_upper(n + 1, d + 1)
    if n <= 0: return 1
    if d <= 0: return 2 ** n
    if d == 2: return 2 ** (n - 1)
    if d > n: return 1
    if 3 * d > 2 * n: return 2
    if (n, d) in OVERRIDE: return OVERRIDE[(n, d)]
    if (n, d) in _AB: return _AB[(n, d)]
    if n < 6: raise ValueError
    raise KeyError((n, d))

# ---- constant weight ----
def _load_cw():
    cw = {}
    for ctx, t in htmltab.tables(os.path.join(H, 'ref', 'agrell_cw.html')):
        m = re.search(r'A\(n,(\d+),w\)\s*$', ctx.strip())
        if not m: continue
        d = int(m.group(1))
        hdr = t[0]
        ws = []
        for c in hdr[1:]:
            v = _num(c[2])
            ws.append(v[0] if v else None)
        for row in t[1:]:
            if not row or not _num(row[0][2]): continue
            n = _num(row[0][2])[0]
            col = 0
            for c in row[1:]:
                if c[0] == 'th': break
                if col < len(ws) and ws[col] is not None and _num(c[2]):
                    v = _num(c[2])
                    cw[(n, d, ws[col])] = v[-1]
                col += c[1]
    return cw
_CW = _load_cw()

@lru_cache(None)
def cw_upper(n, d, w):
    if d % 2: d += 1
    if w < 0 or w > n: return 0
    w = min(w, n - w)
    if d <= 0: return math.comb(n, w)
    if d <= 2: return math.comb(n, w)
    if 2 * w < d: return 1
    if 2 * w == d: return n // w if w > 0 else 1
    best = None
    if (n, d, w) in _CW: best = _CW[(n, d, w)]
    # Johnson bound recursion (always valid)
    j = (n * cw_upper(n - 1, d, w - 1)) // w
    j2 = (n * cw_upper(n - 1, d, w)) // (n - w) if n - w > 0 else j
    j = min(j, j2)
    return min(best, j) if best is not None else j

# ---- doubly constant weight T(w1,n1,w2,n2,d) ----
def _load_dcw():
    T = {}
    for d in (4, 6, 8, 10, 12, 14):
        p = os.path.join(H, 'ref', 'agrell_dcw%d.html' % d)
        for ctx, t in htmltab.tables(p):
            ms = re.findall(r'T\((\d+),n1,(\d+),n2,(\d+)\)', ctx)
            if not ms: continue
            w1, w2, dd = map(int, ms[-1]); assert dd == d
            hdr = t[0]
            n2s = []
            for c in hdr[1:]:
                v = _num(c[2]); n2s.append(v[0] if v else None)
            for row in t[1:]:
                if not row or not _num(row[0][2]): continue
                n1 = _num(row[0][2])[0]
                col = 0
                for c in row[1:]:
                    if c[0] == 'th': break
                    if col < len(n2s) and n2s[col] is not None and _num(c[2]):
                        T[(w1, n1, w2, n2s[col], d)] = _num(c[2])[-1]
                    col += c[1]
    return T
_T = _load_dcw()

@lru_cache(None)
def T_upper(w1, n1, w2, n2, d):
    """max size of code of length n1+n2, weight w1 on first n1 coords, w2 on last n2, min dist >= d."""
    if d % 2: d += 1
    if w1 < 0 or w2 < 0 or w1 > n1 or w2 > n2: return 0
    w1 = min(w1, n1 - w1); w2 = min(w2, n2 - w2)
    if d <= 2: return math.comb(n1, w1) * math.comb(n2, w2)
    if w1 == 0: return cw_upper(n2, d, w2)
    if w2 == 0: return cw_upper(n1, d, w1)
    if 2 * (w1 + w2) < d: return 1
    cands = [cw_upper(n1 + n2, d, w1 + w2)]
    # each part: fixing the n2-part, first parts form cw code with distance d - 0 ... use projection:
    # words sharing the same second part must differ by >= d in first part
    cands.append(cw_upper(n1, d, w1) * math.comb(n2, w2))
    cands.append(cw_upper(n2, d, w2) * math.comb(n1, w1))
    # Johnson-type: count incidences of first-part coordinates
    cands.append((n1 * T_upper(w1 - 1, n1 - 1, w2, n2, d)) // w1)
    cands.append((n2 * T_upper(w1, n1, w2 - 1, n2 - 1, d)) // w2)
    for key in [(w1, n1, w2, n2, d), (w2, n2, w1, n1, d)]:
        if key in _T: cands.append(_T[key])
    return min(cands)

if __name__ == '__main__':
    print(len(_CW), len(_T))
    print([cw_upper(17, 4, w) for w in range(3, 9)])
    print([cw_upper(19, 4, w) for w in range(3, 10)])
    print([cw_upper(20, 6, w) for w in range(4, 11)])
    print(T_upper(1, 3, 2, 25, 4), T_upper(1,2,2,4,4), T_upper(2,4,2,8,4))
    print(A_upper(19,4), A_upper(18,5), A_upper(25,8))
