"""Exact cert: plain Schrijver + B1 rows (block + general V from a cp json) [+ column cap row with premise A(n-1,d) <= Ap].
usage: make_cert2.py n d M Ap outdir [cp.json]"""
import sys, os, json
from cp import build
import sdpcore, certio
n, d, M, Ap, outdir = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
extra = json.load(open(sys.argv[6]))['extra'] if len(sys.argv) > 6 else []
if len(sys.argv) > 7: extra = extra[:int(sys.argv[7])]
K = 60
os.makedirs(outdir, exist_ok=True)
P, meta = build(n, d, M, Ap, [tuple(e) for e in extra])
if os.environ.get('SOLVER') == 'dd':
    sdpcore.CPUS = os.environ.get('CPU')
    info = P.solve('agents/novaB/runs/cert2dd', threads=1)
else:
    info = sdpcore.solve_csdp(P, 'agents/novaB/runs/cert2', threads=1, cpus=os.environ.get('CPU'))
cert = P.make_certificate(info, K=K)
b = P.check_certificate(cert)
print(info['phase.value'][:20], 'num', info['bound'], 'rigorous', float(b), 'EXCLUDED' if b < M else 'NOT EXCLUDED', flush=True)
prob = {'n': n, 'd': d, 'M': M, 'case': 'oa'}
if Ap: prob['premise_A_n1_d_le'] = Ap
path = os.path.join(outdir, 'excl_n%d_d%d_M%d_oa%s.json' % (n, d, M, '_cap%d' % Ap if Ap else ''))
out = certio.export(P, cert, K, path, prob)
j = 0
for z, m in zip(cert['z'], meta):
    if z == 0: continue
    if m and m[0] == 'cap':
        out['lin'][j]['cap'] = m[1]
        if m[2]: out['lin'][j]['cap_eps'] = str(m[2])
    elif m: out['lin'][j]['oa_cols'] = {'r': m[1], 'mult': list(m[2])}
    j += 1
json.dump(out, open(path, 'w'))
print('wrote', path, 'rows', len(out['lin']))
