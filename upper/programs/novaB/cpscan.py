"""usage: cpscan.py n d M0 Ap rounds  -- cutting-plane exclusion scan downward; keeps all cuts found (they are valid for any M,
re-instantiated with the new rho); stops at the first M not excluded."""
import sys, os, random, json
from cp import build, separate, wd_of
from oa import betas
import sdpcore
n, d, M, Ap, rounds = map(int, sys.argv[1:6])
random.seed(int(os.environ.get("SEED", "7")))
pool = []; keys = set()
work = "agents/novaB/runs/cps_%d" % n
res = {}
while True:
    ok = False
    for it in range(rounds):
        P, meta = build(n, d, M, Ap, pool)
        info = sdpcore.solve_csdp(P, work, threads=1, cpus=os.environ.get("CPU"))
        if info['bound'] is not None and info['bound'] < M and info['phase.value'].startswith('OK'):
            ok = True; break
        A, beta = betas(P, n, M, info['x'])
        cuts = separate(beta, n, M)
        new = 0
        for v, r, mult, wd in cuts:
            key = (r, tuple(mult))
            if key in keys: continue
            keys.add(key); pool.append((r, mult)); new += 1
        if new == 0:
            cuts = separate(beta, n, M, tries=400, steps=800)
            for v, r, mult, wd in cuts:
                key = (r, tuple(mult))
                if key in keys: continue
                keys.add(key); pool.append((r, mult)); new += 1
            if new == 0: break
    print(n, d, M, 'Ap', Ap, 'rounds', it + 1, 'pool', len(pool), 'bound', info['bound'], 'EXCL' if ok else 'feas', flush=True)
    json.dump({'n': n, 'd': d, 'Ap': Ap, 'extra': pool}, open('runs/cpscan_%d_%d_c%s.json' % (n, d, os.environ.get('CPU')), 'w'))
    if not ok: break
    M -= 1
    if M < int(os.environ.get('MSTOP', '0')): break
