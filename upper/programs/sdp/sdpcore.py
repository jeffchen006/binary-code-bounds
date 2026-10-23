"""Generic exact-data SDP container + SDPA-GMP driver + rigorous rational certificate.

Primal (the code relaxation):   max  c0 + sum_v c_v x_v
  s.t.  A_k(x) = A_k0 + sum_v x_v A_kv  PSD    (blocks k, integer/rational data)
        g_l(x) = g_l0 + sum_v g_lv x_v >= 0    (linear)
        lb_v <= x_v <= ub_v                    (box, valid for every code)

Certificate: rational Y_k = V_k V_k^T (PSD by construction), rational z_l >= 0.
For every feasible x:
  c(x) <= c(x) + sum <Y_k, A_k(x)> + sum z_l g_l(x) = K + sum_v r_v x_v
        <= K + sum_v max(r_v*lb_v, r_v*ub_v).
"""
from fractions import Fraction
import math, os, subprocess, re, json, sys
import mpmath

mpmath.mp.dps = 60
HERE = os.path.dirname(os.path.abspath(__file__))
SOLVER = os.environ.get('SDP_SOLVER', 'dd')
SDPA = os.path.join(HERE, '..', '..', 'opt', 'sdpa-gmp', {'gmp': 'sdpa_gmp', 'qd': 'sdpa_qd4'}.get(SOLVER, 'sdpa_dd2'))
CPUS = os.environ.get('SDP_CPUS')


class Problem:
    def __init__(self, name):
        self.name = name
        self.vars = []          # list of keys
        self.vidx = {}
        self.c0 = Fraction(0)
        self.c = {}             # v -> Fraction
        self.blocks = []        # list of dict: {'size':m, 'scale':[floats or None], 'ent':{v_or_None: {(i,j):Fraction}}}
        self.lin = []           # list of (g0, {v:coef})
        self.lb = {}
        self.ub = {}

    def var(self, key, lb=0, ub=1):
        if key not in self.vidx:
            self.vidx[key] = len(self.vars)
            self.vars.append(key)
            self.lb[key] = Fraction(lb)
            self.ub[key] = Fraction(ub)
        return key

    def add_block(self, size, ent, scale=None):
        """ent: {v or None(const): {(i,j): value}} with i<=j (upper triangle, symmetric)."""
        self.blocks.append({'size': size, 'ent': ent, 'scale': scale})

    def add_lin(self, g0, coefs):
        self.lin.append((Fraction(g0), {v: Fraction(a) for v, a in coefs.items() if a != 0}))

    # ---------------- SDPA output ----------------
    def write_sdpa(self, path):
        m = len(self.vars)
        # linear constraints + box as one diagonal block
        diag = []
        for g0, co in self.lin:
            diag.append((g0, co))
        for v in self.vars:
            diag.append((-self.lb[v], {v: Fraction(1)}))
            diag.append((self.ub[v], {v: Fraction(-1)}))
        self._diag = diag
        nb = len(self.blocks) + 1
        bs = [b['size'] for b in self.blocks] + [-len(diag)]
        lines = ['"%s"' % self.name, str(m), str(nb), ' '.join(map(str, bs))]
        # SDPA: min sum cS_i x_i ; cS = -c
        vs = getattr(self, 'vscale', {})
        lines.append(' '.join(mpmath.nstr(-mpmath.mpf(self.c.get(v, 0).numerator) / self.c.get(v, Fraction(0)).denominator * vs.get(v, 1), 40) if self.c.get(v, 0) else '0' for v in self.vars))
        out = []
        for bi, b in enumerate(self.blocks):
            sc = b['scale']
            for v, mat in b['ent'].items():
                matno = 0 if v is None else self.vidx[v] + 1
                sgn = -1 if v is None else 1   # F0 = -A0
                for (i, j), val in mat.items():
                    if val == 0:
                        continue
                    x = mpmath.mpf(val.numerator) / val.denominator
                    if v is not None:
                        x = x * vs.get(v, 1)
                    if sc is not None:
                        x = x * sc[i] * sc[j]
                    out.append('%d %d %d %d %s' % (matno, bi + 1, i + 1, j + 1, mpmath.nstr(sgn * x, 45)))
        bi = len(self.blocks)
        for r, (g0, co) in enumerate(diag):
            if g0 != 0:
                out.append('0 %d %d %d %s' % (bi + 1, r + 1, r + 1, mpmath.nstr(-mpmath.mpf(g0.numerator) / g0.denominator, 45)))
            for v, a in co.items():
                out.append('%d %d %d %d %s' % (self.vidx[v] + 1, bi + 1, r + 1, r + 1, mpmath.nstr(mpmath.mpf(a.numerator) / a.denominator * vs.get(v, 1), 45)))
        with open(path, 'w') as f:
            f.write('\n'.join(lines + out) + '\n')

    def solve(self, workdir, param=None, threads=8):
        os.makedirs(workdir, exist_ok=True)
        dat = os.path.join(workdir, self.name + '.dat-s')
        res = os.path.join(workdir, self.name + '.out')
        self.write_sdpa(dat)
        param = param or os.path.join(HERE, {'gmp': 'param.sdpa', 'qd': 'param_qd.sdpa'}.get(SOLVER, 'param_dd.sdpa'))
        env = dict(os.environ, OMP_NUM_THREADS=str(threads))
        cmd = [SDPA, '-ds', dat, '-o', res, '-p', param]
        subprocess.run(cmd, env=env, stdout=open(res + '.log', 'w'), stderr=subprocess.STDOUT)
        info = self.read_result(res)
        try:
            info['bound'] = float(self.c0) - float(info['objValDual'])
            info['pbound'] = float(self.c0) - float(info['objValPrimal'])
        except Exception:
            info['bound'] = None
        return info

    def read_result(self, res):
        txt = open(res).read()
        if '_' in txt:
            def _comb(m):
                return mpmath.nstr(sum(mpmath.mpf(p) for p in m.group(0).split('_')), 60)
            txt = re.sub(r'[-+]?\d\.\d+[eE][-+]\d+(?:_[-+]?\d\.\d+[eE][-+]\d+)+', _comb, txt)
        def grab(name):
            i = txt.index(name)
            j = txt.index('{', i)
            depth = 0
            for k in range(j, len(txt)):
                if txt[k] == '{':
                    depth += 1
                elif txt[k] == '}':
                    depth -= 1
                    if depth == 0:
                        return txt[j:k + 1]
        info = {}
        for key in ['phase.value', 'objValPrimal', 'objValDual']:
            mm = re.search(key + r'\s*=\s*(\S+)', txt)
            info[key] = mm.group(1) if mm else None
        xv = [mpmath.mpf(s) for s in re.findall(r'[-+]?\d\.\d+[eE][-+]\d+|[-+]?\d+\.\d*', grab('xVec'))]
        ym = grab('yMat')
        # parse nested blocks
        blocks = []
        body = ym.strip()[1:-1]
        # split top-level groups
        depth = 0; cur = ''
        groups = []
        for ch in body:
            if ch == '{':
                depth += 1
            if depth > 0:
                cur += ch
            if ch == '}':
                depth -= 1
                if depth == 0:
                    groups.append(cur); cur = ''
        for g in groups:
            g = g.strip()
            if g.startswith('{{') or g.count('{') > 1:
                rows = re.findall(r'\{([^{}]*)\}', g)
                M = [[mpmath.mpf(s) for s in r.replace(',', ' ').split()] for r in rows]
                blocks.append(('m', M))
            else:
                vals = [mpmath.mpf(s) for s in g.strip('{}').replace(',', ' ').split()]
                blocks.append(('d', vals))
        info['x'] = xv
        info['Y'] = blocks
        return info

    # ---------------- certificate ----------------
    def make_certificate(self, info, K=200):
        """Round dual to rationals. Y_k = V V^T with V = Q sqrt(lam+) (unscaled), z>=0."""
        Ys = info['Y']
        cert = {'V': [], 'z': []}
        for bi, b in enumerate(self.blocks):
            typ, M = Ys[bi]
            m = b['size']
            if typ == 'd':
                assert len(M) == m == 1, (bi, m, len(M))
                M = [M]
            A = mpmath.matrix(m, m)
            for i in range(m):
                for j in range(m):
                    A[i, j] = M[i][j]
            sc = b['scale']
            if sc is not None:
                for i in range(m):
                    for j in range(m):
                        A[i, j] = A[i, j] * sc[i] * sc[j]
            A = (A + A.T) / 2
            E, Q = mpmath.eigsy(A)
            V = []
            for i in range(m):
                row = []
                for r in range(m):
                    if E[r] > 0:
                        row.append(Q[i, r] * mpmath.sqrt(E[r]))
                    else:
                        row.append(mpmath.mpf(0))
                V.append(row)
            # drop zero columns
            cols = [r for r in range(m) if E[r] > 0]
            Vr = [[_toq(V[i][r], K) for r in cols] for i in range(m)]
            cert['V'].append(Vr)
        typ, d = Ys[len(self.blocks)]
        nl = len(self.lin)
        cert['z'] = [max(_toq(d[r], K), Fraction(0)) for r in range(nl)]
        return cert

    def check_certificate(self, cert):
        """Exact: returns Fraction upper bound valid for all feasible x."""
        K = self.c0
        r = {v: self.c.get(v, Fraction(0)) for v in self.vars}
        for bi, b in enumerate(self.blocks):
            V = cert['V'][bi]
            m = b['size']
            Y = [[sum((V[i][s] * V[j][s] for s in range(len(V[i]))), Fraction(0)) for j in range(m)] for i in range(m)]
            for v, mat in b['ent'].items():
                tot = Fraction(0)
                for (i, j), val in mat.items():
                    tot += val * (Y[i][j] if i == j else 2 * Y[i][j])
                if v is None:
                    K += tot
                else:
                    r[v] += tot
        for z, (g0, co) in zip(cert['z'], self.lin):
            assert z >= 0
            if z == 0:
                continue
            K += z * g0
            for v, a in co.items():
                r[v] += z * a
        bound = K
        for v in self.vars:
            bound += max(r[v] * self.lb[v], r[v] * self.ub[v])
        return bound


def _toq(x, K):
    """mpf -> exact Fraction with denominator 2^K (rounded)."""
    return Fraction(int(mpmath.nint(x * mpmath.mpf(2) ** K)), 2 ** K)


def binom(n, k):
    return math.comb(n, k) if 0 <= k <= n else 0


CSDP = os.path.join(HERE, '..', '..', 'opt', 'sdpa-gmp', 'csdp')


def solve_csdp(P, workdir, threads=4, cpus=None):
    """Solve with CSDP (double). Returns info dict like Problem.solve (x, Y in same layout)."""
    import numpy as np
    os.makedirs(workdir, exist_ok=True)
    dat = os.path.join(workdir, P.name + '.dat-s')
    sol = os.path.join(workdir, P.name + '.csdp.sol')
    P.write_sdpa(dat)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    cmd = [CSDP, dat, sol]
    out = subprocess.run(cmd, env=env, capture_output=True, text=True).stdout
    open(sol + '.log', 'w').write(out)
    info = {}
    mm = re.search(r'Primal objective value:\s*(\S+)', out); info['objValPrimal'] = mm.group(1) if mm else None
    mm = re.search(r'Dual objective value:\s*(\S+)', out); info['objValDual'] = mm.group(1) if mm else None
    mm = re.search(r'Success: (.*)', out); info['phase.value'] = ('OK ' + mm.group(1)[:30]) if mm else ('FAIL ' + out.strip().splitlines()[-1][:60] if out.strip() else 'FAIL')
    lines = open(sol).read().split('\n')
    xv = [mpmath.mpf(s) for s in lines[0].split()]
    sizes = [b['size'] for b in P.blocks] + [len(P._diag)]
    Ys = [np.zeros((s, s)) for s in sizes]
    for ln in lines[1:]:
        p = ln.split()
        if len(p) != 5 or p[0] != '2':
            continue
        b, i, j, v = int(p[1]) - 1, int(p[2]) - 1, int(p[3]) - 1, float(p[4])
        Ys[b][i, j] = v; Ys[b][j, i] = v
    Y = []
    for b, s in enumerate(sizes):
        if b < len(P.blocks):
            Y.append(('m', [[mpmath.mpf(Ys[b][i, j]) for j in range(s)] for i in range(s)]))
        else:
            Y.append(('d', [mpmath.mpf(Ys[b][i, i]) for i in range(s)]))
    info['x'] = xv; info['Y'] = Y
    # CSDP: primal = max tr(C X) (= our dual), objective values reported with SDPA sign convention
    try:
        info['bound'] = float(P.c0) - float(info['objValDual'])
        info['pbound'] = float(P.c0) - float(info['objValPrimal'])
    except Exception:
        info['bound'] = None
    return info


def equilibrate(P):
    """Fold a diagonal row scaling into each block's 'scale' so that diagonal magnitudes are ~1 (congruence)."""
    vs = getattr(P, 'vscale', {})
    for b in P.blocks:
        m = b['size']
        sc = b['scale'] if b['scale'] is not None else [mpmath.mpf(1)] * m
        diag = [mpmath.mpf(0)] * m
        for v, mat in b['ent'].items():
            w = 1 if v is None else vs.get(v, 1)
            for (i, j), val in mat.items():
                if i == j:
                    diag[i] += abs(mpmath.mpf(val.numerator) / val.denominator) * sc[i] ** 2 * w
        b['scale'] = [sc[i] / mpmath.sqrt(diag[i]) if diag[i] > 0 else sc[i] for i in range(m)]
