import sys, time
from oa import *
import sdpcore
n, d, M = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
P = split.build(n, d, 0, name='plain_%d_%d' % (n, d))
t = time.time()
info = sdpcore.solve_csdp(P, "agents/novaB/runs/t1", threads=1, cpus=None)
print(info['phase.value'], info['bound'], info['pbound'], '%.0fs' % (time.time() - t))
A, b = betas(P, n, M, info['x'])
print('A', {a: round(v, 4) for a, v in A.items() if abs(v) > 1e-9})
print('beta', {k: round(v, 4) for k, v in b.items()})
