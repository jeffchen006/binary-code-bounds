"""min-beta1 exclusion: SDP + B1 rows + row |C| >= M, objective max -g(x), g = K_1(0) + sum_i A_i K_1(i) = sum_i S_{e_i}^2 / M.
Premise A(n-1,d) <= Ap gives g <= capg = n (2Ap - M)^2 / M. If rigorous max(-g) < -capg, no code of size M exists.
usage: mb1.py n d M Ap rounds [pool.json] ; env CERT=outdir to write a certificate"""
import sys, os, json, random
from fractions import Fraction
from cp import build, separate
from oa import betas, shell_form
import sdpcore, certio
n, d, M, Ap, R = map(int, sys.argv[1:6])
pool = [tuple(e) for e in json.load(open(sys.argv[6]))['extra']] if len(sys.argv) > 6 else []
random.seed(3); keys = set((r, tuple(m)) for r, m in pool)
capg = Fraction(n * (2 * Ap - M) ** 2, M)
for it in range(R):
    P, meta = build(n, d, M, 0, pool)
    # size >= M row
    P.lin.append((P.c0 - M, dict(P.c))); meta.append(('size', M))
    c0, f = shell_form(P, n, 1)
    P.c0 = -c0; P.c = {v: -c for v, c in f.items()}
    P.name = 'mb1_%d_%d_%d' % (n, d, M)
    info = sdpcore.solve_csdp(P, 'agents/novaB/runs/mb1', threads=1, cpus=os.environ.get('CPU'))
    ok = info['phase.value'].startswith('OK')
    gmin = -info['bound'] if info['bound'] is not None else None
    print('round', it, info['phase.value'][:22], 'min g', gmin, 'capg', float(capg), 'EXCL' if (ok and gmin > capg) else '', flush=True)
    if ok and gmin > capg:
        cert = P.make_certificate(info, K=60); B = P.check_certificate(cert)
        print('rigorous min g >=', float(-B), 'EXCLUDED' if -B > capg else 'not excluded', flush=True)
        if -B > capg and os.environ.get('CERT'):
            out = os.environ['CERT']; os.makedirs(out, exist_ok=True)
            path = os.path.join(out, 'excl_n%d_d%d_M%d_minbeta1_cap%d.json' % (n, d, M, Ap))
            o = certio.export(P, cert, 60, path, {'n': n, 'd': d, 'M': M, 'case': 'oa', 'mode': 'minbeta1', 'premise_A_n1_d_le': Ap})
            j = 0
            for z, m in zip(cert['z'], meta):
                if z == 0: continue
                if m and m[0] == 'size': o['lin'][j]['size_ge_M'] = True
                elif m: o['lin'][j]['oa_cols'] = {'r': m[1], 'mult': list(m[2])}
                j += 1
            json.dump(o, open(path, 'w')); print('wrote', path)
        break
    A, beta = betas(P, n, M, info['x'])
    new = 0
    for v, r, mult, wd in separate(beta, n, M):
        if (r, tuple(mult)) in keys: continue
        keys.add((r, tuple(mult))); pool.append((r, mult)); new += 1
    print('  new cuts', new, flush=True)
    if new == 0: break
json.dump({'extra': pool}, open('runs/mb1_%d_%d_%d.json' % (n, d, M), 'w'))
