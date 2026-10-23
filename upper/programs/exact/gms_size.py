import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gms
from fractions import Fraction
n, d = int(sys.argv[1]), int(sys.argv[2])
t = time.time()
P = gms.build(n, d, Fraction(2 ** n))
nnz = sum(len(m) for b in P.blocks for m in b['ent'].values())
print(n, d, 'vars', len(P.vars), 'blocks', len(P.blocks), 'maxblk', max(b['size'] for b in P.blocks),
      'sum blk^2', sum(b['size'] ** 2 for b in P.blocks), 'nnz', nnz, '%.0fs' % (time.time() - t), flush=True)
