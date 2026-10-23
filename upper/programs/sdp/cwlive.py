"""L2: live constant-weight UPPER bounds A(n,d,w) from Brouwer's table (aeb.win.tue.nl/codes/Andw.html, fetched 2026-10-23,
copy in ext/andw_20261023/; upper values incl. Polak quadruple SDP 'Po', KKT, Ostergard etc.). A cell 'a-b' means lower a,
upper b; a single number is exact. apply() lowers bounds._CW to min(Agrell 2000 table, live upper) and clears the cache."""
import re, os, html
import bounds

PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ext', 'andw_20261023', 'Andw.html')


def parse(path=PATH):
    s = open(path, encoding='utf-8').read()
    out = {}
    for m in re.finditer(r'name="d(\d+)">Bounds on A\(n,\d+,w\)</a></h1>\s*<table[^>]*>(.*?)</table>', s, re.S):
        d = int(m.group(1)); tab = m.group(2)
        rows = re.findall(r'<tr>(.*?)</tr>', tab, re.S)
        ws = [int(x) for x in re.findall(r'<th>(\d+)</th>', rows[0])]
        for r in rows[1:]:
            mm = re.match(r'\s*<th>(\d+)</th>(.*)', r, re.S)
            if not mm:
                continue
            n = int(mm.group(1))
            cells = re.findall(r'<td[^>]*>(.*?)</td>', mm.group(2), re.S)
            for w, c in zip(ws, cells):
                c = re.sub(r'<sup>.*?</sup>', '', c, flags=re.S)
                c = html.unescape(re.sub(r'<[^>]+>', '', c)).strip()
                nums = re.findall(r'\d+', c)
                if not nums:
                    continue
                out[(n, d, w)] = int(nums[-1])      # 'a-b' -> b ; 'a' -> a (exact)
    return out


def apply():
    live = parse()
    changed = {}
    for (n, d, w), u in live.items():
        k = (n, d, min(w, n - w))
        old = bounds._CW.get(k)
        if old is None or u < old:
            bounds._CW[k] = u
            changed[k] = (old, u)
    bounds.cw_upper.cache_clear()
    return changed


if __name__ == '__main__':
    ch = apply()
    print(len(parse()), 'live cells;', len(ch), 'lowered vs Agrell table')
    for k in sorted(ch):
        if k[0] <= 28:
            print(k, ch[k])
