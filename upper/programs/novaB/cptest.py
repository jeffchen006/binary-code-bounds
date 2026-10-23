"""usage: cptest.py n d rounds M1 M2 ...  -- cutting-plane exclusion test at fixed M (no cap)."""
import sys, os, random
from cp import build, separate
from oa import betas
import sdpcore
n, d, R = map(int, sys.argv[1:4])
for M in map(int, sys.argv[4:]):
    random.seed(1); pool = []; keys = set()
    for it in range(R):
        P, meta = build(n, d, M, 0, pool)
        info = sdpcore.solve_csdp(P, "agents/novaB/runs/cpt_%d" % n, threads=1, cpus=os.environ.get("CPU"))
        if info['bound'] < M: break
        A, beta = betas(P, n, M, info['x'])
        new = 0
        for v, r, mult, wd in separate(beta, n, M, rs=range(2, min(n, 9))):
            if (r, tuple(mult)) in keys: continue
            keys.add((r, tuple(mult))); pool.append((r, mult)); new += 1
        if new == 0: break
    print('n', n, 'M', M, 'rounds', it + 1, 'pool', len(pool), 'bound', info['bound'], 'EXCL' if info['bound'] < M else 'feas', flush=True)
