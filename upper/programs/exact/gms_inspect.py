import sys, pickle, os
import gms
from gms import pair_key, triple_key, C
from fractions import Fraction
n, d = int(sys.argv[1]), int(sys.argv[2]); out = sys.argv[3]
P = gms.build(n, d, Fraction(80), verbose=False)
info = P.read_result(out)
x = {v: float(info['x'][i]) for i, v in enumerate(P.vars)}
pickle.dump((P.vars, x), open(out + '.x.pkl', 'wb'))
y1 = x[pair_key(n, 0)]
print('y1 =', y1)
D = [j for j in range(d, n + 1) if j % 2 == 0]
for j in D:
    k = pair_key(n, j)
    print('A_%d = %.6f' % (j, C(n, j) * x.get(k, 0) / y1))
import bounds
for j in D:
    print('cw bound A(%d,%d,%d) = %d' % (n, d, j, bounds.cw_upper(n, d, j)))
