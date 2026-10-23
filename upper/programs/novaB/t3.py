import sys, time
from math import comb
from oa import *
import sdpcore

def partitions(total, maxpart=None):
    if maxpart is None: maxpart = total
    if total == 0:
        yield (); return
    for p in range(min(total, maxpart), 0, -1):
        for rest in partitions(total - p, p):
            yield (p,) + rest

def block_wd(blocks):
    wd = {0: 1}
    for a in blocks:
        nw = dict(wd)
        for w, c in wd.items():
            nw[w + a] = nw.get(w + a, 0) + c
        wd = nw
    wd.pop(0); return wd

def block_rows(n, M, rmin=2):
    seen = set()
    for tot in range(1, n + 1):
        for part in partitions(tot):
            r = len(part)
            if r < rmin or M % 2 ** r == 0: continue
            wd = block_wd(part)
            key = (r, tuple(sorted(wd.items())))
            if key in seen: continue
            seen.add(key); yield r, wd

if __name__ == '__main__':
    n, d = int(sys.argv[1]), int(sys.argv[2])
    for M in map(int, sys.argv[3:]):
        P = split.build(n, d, 0, name='oab_%d_%d_%d' % (n, d, M))
        nr = 0
        for r, wd in block_rows(n, M):
            nr += add_oa_row(P, n, M, wd, r)
        t = time.time()
        info = sdpcore.solve_csdp(P, "agents/novaB/runs/t3", threads=1, cpus=None)
        A, b = betas(P, n, M, info['x'])
        print(M, M % 64, nr, info['phase.value'], info['bound'], '%.0fs' % (time.time() - t), {k: round(v, 2) for k, v in b.items() if 0 < k < n and abs(v) > 1e-4}, flush=True)
