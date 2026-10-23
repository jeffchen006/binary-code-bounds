import sys, os, json
from cp import build
from oa import betas
import sdpcore
n, d, M, Ap = map(int, sys.argv[1:5]); pool = [tuple(e) for e in json.load(open(sys.argv[5]))['extra']] if len(sys.argv) > 5 else []
P, meta = build(n, d, M, Ap, pool)
info = sdpcore.solve_csdp(P, "agents/novaB/runs/t6", threads=1, cpus=os.environ.get("CPU"))
A, b = betas(P, n, M, info['x'])
print(info['bound']); print('A', {a: round(v, 3) for a, v in A.items() if abs(v) > 1e-6}); print('beta', {k: round(v, 3) for k, v in b.items() if abs(v) > 1e-4})
