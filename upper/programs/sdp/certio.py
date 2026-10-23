"""Export a rational dual certificate (from sdpcore.make_certificate) for a Schrijver (n1=0) problem to JSON."""
import json
from fractions import Fraction


def key2str(v):
    return ','.join(str(p[0]) for p in v)


def export(P, cert, K, path, problem):
    D = 2 ** K
    out = {'problem': problem, 'K': K, 'blocks': [], 'lin': []}
    for b, V in zip(P.blocks, cert['V']):
        meta = b['meta']
        assert meta['kind'] == 'schr'
        Vi = []
        for row in V:
            r = []
            for x in row:
                y = x * D
                assert y.denominator == 1
                r.append(int(y))
            Vi.append(r)
        out['blocks'].append({'k': meta['k'], 'typ': meta['typ'], 'rows': meta['rows'], 'V': Vi})
    for z, (g0, co) in zip(cert['z'], P.lin):
        if z == 0:
            continue
        y = z * D
        assert y.denominator == 1
        out['lin'].append({'z': int(y), 'g0': str(g0), 'coef': {key2str(v): str(c) for v, c in co.items()}})
    with open(path, 'w') as f:
        json.dump(out, f)
    return out


def gkey(v):
    if v == ('p',):
        return 'p'
    if v[0] == 'Q':
        return 'Q,' + ','.join(str(x) for x in v[1:])
    return 'T,' + ','.join(str(p[0]) for p in v)


def export_gms(P, cert, K, path, problem):
    D = 2 ** K
    out = {'problem': problem, 'K': K, 'blocks': [], 'lin': []}
    for b, V in zip(P.blocks, cert['V']):
        meta = dict(b['meta'])
        meta['rows'] = [list(r) if isinstance(r, tuple) else r for r in meta['rows']]
        Vi = [[int(x * D) for x in row] for row in V]
        assert all((x * D).denominator == 1 for row in V for x in row)
        out['blocks'].append({'meta': meta, 'V': Vi})
    for z, (g0, co) in zip(cert['z'], P.lin):
        if z == 0:
            continue
        y = z * D
        assert y.denominator == 1
        out['lin'].append({'z': int(y), 'g0': str(g0), 'coef': {gkey(v): str(c) for v, c in co.items()}})
    with open(path, 'w') as f:
        json.dump(out, f)
    return out
