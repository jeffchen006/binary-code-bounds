#!/usr/bin/env python3
"""Driver: classify even (n,M,d) codes level by level with ext, check by double counting.

usage: classify.py d n M [--cpus 62-71] [--data DIR]
  classifies even (n,M',d) codes for all M' >= M (up to Amax), recursively
  classifying the needed seeds of length n-1 first.
"""
import sys, os, subprocess, math, argparse
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.join(HERE, os.environ.get('EXTBIN', 'ext'))


def plotkin_even(n, d):
    # upper bound on A(n,d), d even, used only to cap the classification range
    if 2 * d > n:
        return max(1, 2 * (d // (2 * d - n)))
    if 2 * d == n:
        return 4 * d
    return 1 << n


def parse_cpus(s):
    out = []
    for part in s.split(','):
        if '-' in part:
            a, b = part.split('-'); out += list(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


def fname(data, n, M):
    return os.path.join(data, f'n{n}_M{M}.txt')


def count_lines(f):
    with open(f) as fh:
        return sum(1 for l in fh if l.strip())


def run_level(args, d, n, M, amax_prev):
    out = fname(args.data, n, M)
    if os.path.exists(out) and os.path.exists(out + '.dc'):
        return
    lo = (M + 1) // 2
    seeds = []
    for Mp in range(lo, min(M, amax_prev) + 1):
        f = fname(args.data, n - 1, Mp)
        with open(f) as fh:
            seeds += [l for l in fh if l.strip()]
    cpus = parse_cpus(args.cpus)
    nj = min(len(seeds), max(1, len(cpus) * 30)) if seeds else 1
    chunks = [seeds[i::nj] for i in range(nj)]
    tmpdir = os.path.join(args.data, 'tmp'); os.makedirs(tmpdir, exist_ok=True)
    import queue, threading
    from concurrent.futures import ThreadPoolExecutor
    cpuq = queue.Queue()
    for c in cpus: cpuq.put(c)
    def work(j):
        ch = chunks[j]
        inf = os.path.join(tmpdir, f'n{n}_M{M}_in{j}')
        with open(inf, 'w') as fh:
            fh.writelines(ch)
        env = dict(os.environ, STATS=os.path.join(tmpdir, f'n{n}_M{M}_st{j}'))
        if args.final is not None and n == args.top:
            env['FINAL'] = str(args.final)
        c = cpuq.get()
        try:
            with open(inf) as fi, open(os.path.join(tmpdir, f'n{n}_M{M}_out{j}'), 'w') as fo:
                r = subprocess.run([EXT, str(n), str(d), str(M)],
                                   stdin=fi, stdout=fo, env=env)
        finally:
            cpuq.put(c)
        assert r.returncode == 0
    with ThreadPoolExecutor(max_workers=len(cpus)) as ex:
        list(ex.map(work, range(nj)))
    # merge and double count
    lhs = Fraction(0); rhs = Fraction(0); ncodes = 0; hist = {}
    Gn = 2 ** n * math.factorial(n)
    with open(out + '.part', 'w') as fo:
        for j in range(nj):
            with open(os.path.join(tmpdir, f'n{n}_M{M}_out{j}')) as fh:
                for l in fh:
                    fo.write(l)
            with open(os.path.join(tmpdir, f'n{n}_M{M}_st{j}')) as fh:
                for l in fh:
                    t = l.split()
                    rhs += Fraction(t[2]) / int(t[3])
                    ncodes += int(t[5])
                    lhs += Fraction(int(t[7]), Gn)
                    if len(t) > 8 and t[8] == 'H':
                        for h in t[9:]:
                            a, b = h.split(':'); hist[int(a)] = hist.get(int(a), 0) + int(b)
    os.rename(out + '.part', out)
    ok = (lhs == rhs)
    with open(out + '.dc', 'w') as fh:
        fh.write(f'n={n} M={M} d={d} codes={ncodes} lhs={lhs} rhs={rhs} ok={ok}\n')
        if hist:
            fh.write('final M+omega histogram: ' + ' '.join(f'{a}:{hist[a]}' for a in sorted(hist)) + '\n')
    print(f'(n={n}, M={M}, d={d}): {ncodes} classes; double count {"OK" if ok else "FAIL"} ({lhs} vs {rhs})' + (f' final M+omega hist {dict(sorted(hist.items()))}' if hist else ''), flush=True)
    for j in range(nj):
        for s in ('in', 'out', 'st'):
            os.remove(os.path.join(tmpdir, f'n{n}_M{M}_{s}{j}'))
    if not ok:
        sys.exit('double counting failed')


def classify(args, d, n, Mmin):
    amax = plotkin_even(n, d)
    if Mmin <= 1:
        # base: single codeword
        f = fname(args.data, n, 1)
        if not os.path.exists(f):
            with open(f, 'w') as fh:
                fh.write('1 0 aut=%d e=%d\n' % (math.factorial(n), n))
            open(f + '.dc', 'w').write('base\n')
        Mmin = 2
    if Mmin > amax:
        return
    amax_prev = plotkin_even(n - 1, d)
    classify(args, d, n - 1, (Mmin + 1) // 2)
    for M in range(Mmin, amax + 1):
        if (M + 1) // 2 > amax_prev:
            break
        run_level(args, d, n, M, amax_prev)
        if open(fname(args.data, n, M) + '.dc').read().split('codes=')[1].split()[0] == '0':
            # no codes of size M: larger sizes do not exist either
            for M2 in range(M + 1, amax + 1):
                open(fname(args.data, n, M2), 'w').close()
                open(fname(args.data, n, M2) + '.dc', 'w').write('empty (subset of empty)\n')
            break


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('d', type=int); ap.add_argument('n', type=int); ap.add_argument('M', type=int)
    ap.add_argument('--cpus', default=None)
    ap.add_argument('--data', default=None)
    ap.add_argument('--final', type=int, default=None, help='at length --top: compute omega of lengthening graph, keep codes with M+omega >= FINAL')
    ap.add_argument('--top', type=int, default=None)
    ap.add_argument('--only', type=int, default=None, help='classify only this size at length n (seeds must exist)')
    a = ap.parse_args()
    if a.data is None:
        a.data = os.path.join(HERE, '..', 'data', f'd{a.d}')
    os.makedirs(a.data, exist_ok=True)
    if a.only is not None:
        run_level(a, a.d, a.n, a.only, plotkin_even(a.n - 1, a.d))
    else:
        classify(a, a.d, a.n, a.M)
