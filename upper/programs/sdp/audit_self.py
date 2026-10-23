"""P1b/P2: self-consistent premises. V3 (Schrijver + Best/cw/MEL/KT rows) where every premise A(m,d), m < n, is OUR value:
min(our certified V2/V3/V4 from this pass, team-audited bound), except entries closed in the table (lower = upper), which keep
the table value. Removes the dependence on disputed/unverified published values (e.g. TLY A(18,4) <= 6551).
Then a parity descent (V4) from the result. Runs n = 17..28 in increasing order (premises need smaller n first).
usage: audit_self.py d1,d2,...   (one process per d; appends to runs/audit/self.jsonl)"""
import sys, os, json, math, time
sys.dont_write_bytecode = True
import bounds, split, parity
from audit import solve_exact, NN

BASE = '.'
B = json.load(open(os.path.join(BASE, 'tools/bounds.json')))['bounds']
prev = {}
for l in open(os.path.join(BASE, 'agents/sdp/runs/audit/results.jsonl')):
    r = json.loads(l)
    prev[(r['n'], r['d'])] = r
OUT = os.path.join(BASE, 'agents/sdp/runs/audit/self.jsonl')


def premise(m, d):
    k = '%d,%d' % (m, d)
    e = B.get(k)
    if e is None or m < 17:
        return None                   # n <= 16: table values are exact and classical
    if e['lower'] == e['upper'] and (m, d) not in ((18, 4),):
        return e['upper']             # closed entries keep the table value
    return SELF.get((m, d))


SELF = {}
for d in [int(x) for x in sys.argv[1].split(',')]:
    for n in range(17, 29):
        if d > n:
            continue
        NN[0] = n
        t0 = time.time()
        bounds.OVERRIDE.clear()
        used = {}
        for m in range(max(d, 6), n):
            p = premise(m, d)
            if p is not None:
                bounds.OVERRIDE[(m, d)] = p; used[m] = p
        e = B['%d,%d' % (n, d)]
        team = e.get('team_upper')
        r = prev.get((n, d), {})
        T = e['upper']
        P = split.build(n, d, 0, opts={}, name='as')
        num, rig, how = solve_exact(P, 'as_%d_%d_lin' % (n, d), T)
        v3 = float(rig) if rig is not None else None
        vals = [v for v in (r.get('schrijver'), v3) if v is not None]
        F = math.floor(min(vals) + 1e-9)
        M = F; steps = []
        while len(steps) < 8 and M >= 1:
            cs = parity.cases(M, n)
            if cs[0][0] == '0mod4':
                break
            ok = True
            for case, w in cs:
                Q = parity.build(n, d, 0, M, case, w, opts={})
                _, rr, _ = solve_exact(Q, 'as_%d_%d_par%d_%s%s' % (n, d, M, case, w), M)
                if rr is None or rr >= M:
                    ok = False; break
            steps.append([M, ok])
            if not ok:
                break
            M -= 1
        ours = min([x for x in (M, team) if x is not None])
        SELF[(n, d)] = ours
        out = {'n': n, 'd': d, 'live': e['upper'], 'team': team, 'premises': used, 'v3_self': v3, 'v4_self': M,
               'v3_table': r.get('lin'), 'v4_table': r.get('parity_bound'), 'self_value': ours, 'time': round(time.time() - t0)}
        open(OUT, 'a').write(json.dumps(out) + '\n')
        print(json.dumps(out), flush=True)
