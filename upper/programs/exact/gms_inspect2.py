import sys, pickle
import gms
from gms import pair_key, triple_key, C
sys.path.insert(0, '../../sdp')
import bounds
n, d = int(sys.argv[1]), int(sys.argv[2]); out = sys.argv[3]
vars_, x = pickle.load(open(out + '.x.pkl', 'rb'))
y1 = x[pair_key(n, 0)]
D = [j for j in range(d, n + 1) if j % 2 == 0]
for i in D:
    y2 = x.get(pair_key(n, i), 0)
    if y2 < 1e-12: continue
    print('--- pair distance', i, 'A_i=%.4f' % (C(n, i) * y2 / y1))
    tot = 0
    for j in D:
        for t in range(max(0, i + j - n), min(i, j) + 1):
            c = i + j - 2 * t
            if c not in D: continue
            k3 = triple_key(n, i, j, t)
            y3 = x.get(k3, 0)
            N = C(i, t) * C(n - i, j - t) * y3 / y2
            if N > 1e-9:
                Tv = bounds.T_upper(t, i, j - t, n - i, d)
                tot += N
                print('  w: d(0,w)=%d d(u,w)=%d (t=%d): N=%.4f  T_upper=%d' % (j, c, t, N, Tv))
    print('  total other codewords: %.4f (should be M-2 = %.4f)' % (tot, y1 - 2))
