"""Build AUDIT_PUBLISHED.md from runs/audit/results.jsonl + live table (tools/bounds.json) + sources (live table text)."""
import json, math, os
BASE = '.'
B = json.load(open(os.path.join(BASE, 'tools/bounds.json')))['bounds']
SRC = {}
def put(tag, items):
    for n, d, v in items:
        SRC.setdefault((n, d), []).append((v, tag))
put('Schrijver 2004 (3-pt SDP)', [(19,6,1280),(19,8,142),(20,8,274),(22,10,87),(23,6,13766),(25,6,48008),(25,8,5477),(25,10,503),
    (26,10,859),(27,8,17768),(28,8,32151)])
put('Schrijver Terwilliger preprint', [(25,6,47998),(26,10,886)])
put('Laurent 2007 (SDP)', [(20,8,273),(25,6,47997)])
put('MEL 2001', [(23,4,173015),(26,6,84260),(27,6,157286),(25,8,5557),(26,8,9673),(26,10,990)])
put('MEL 2002 (LP+rows)', [(21,4,43688),(25,4,599184),(27,6,157285),(26,8,9672),(28,8,32204),(26,10,989)])
put('MEL 2007 (assoc. scheme LP)', [(22,4,87333),(23,4,172361),(25,6,47538)])
put('Haas 2008', [(20,4,26168),(28,4,4792950)])
put('Johnson bound', [(24,4,344308)])
put('AVZ 2001 (2 mod 4)', [(21,4,43689),(22,6,6941)])
put('GMS 2012 (4-pt SDP)', [(18,6,673),(19,6,1237),(20,6,2279),(23,6,13674),(19,8,135),(20,8,256),(25,8,5421),(26,8,9275),
    (27,8,17099),(21,10,47),(22,10,84),(24,10,268),(25,10,466),(26,10,836),(27,10,1585),(28,10,2817),(25,12,55),(26,12,96)])
put('Honkala 1987', [(25,12,56),(26,12,98)])
put('van Pul 1982', [(18,8,72),(21,10,48),(22,10,88)])
put('Kim-Toan 2012 (SDP)', [(18,8,71),(19,8,131)])
put('Ostergard (exact)', [(18,8,64),(19,8,128),(17,8,36)])
put('TLY 2023 (split SDP)', [(18,4,6551)])

rows = [json.loads(l) for l in open(os.path.join(BASE, 'agents/sdp/runs/audit/results.jsonl'))]
rows = {(r['n'], r['d']): r for r in rows}
SELF = {}
for l in open(os.path.join(BASE, 'agents/sdp/runs/audit/self.jsonl')):
    r = json.loads(l); SELF[(r['n'], r['d'])] = r
def fl(x):
    return None if x is None else math.floor(x + 1e-9)
out = ['# Audit of published upper bounds (sdp, round 5, P1)', '',
       'Every value below is an EXACT certificate bound (rational dual, PSD by construction, checked in fractions), not a',
       'solver number. V1 Delsarte LP (even weights); V2 plain Schrijver 3-point; V3 Schrijver + Best/cw/MEL/Best20/Kim-Toan',
       'rows (premises: Brouwer table 2026-10-21); V4 = V3 + parity / 2-mod-4 (AVZ, Kim-Toan) exclusion descent.',
       '"ours" = min(V1, V2, V4 self) (no unverified premise); V3/V4 use table premises and are for comparison only. live = live Brouwer table upper bound; team = audited team bound not yet in the table.',
       'Verdict: OK = ours <= live (published value reproduced or beaten); GAP = ours > live (published value not reproduced',
       'by these relaxations; the source column says whether that is expected).',
       'V4 self = V3+V4 with self-consistent premises: A(m,d), m < n, = our own certified value or a team-audited bound',
       '(closed table entries keep their value); removes dependence on unverified published values (e.g. TLY 6551).', '',
       '| n | d | live | team | source of live value | V1 LP | V2 Schr | V3 +rows | V4 parity | V4 self | ours | verdict |',
       '|---|---|---|---|---|---|---|---|---|---|---|---|']
gaps = []
for (n, d) in sorted(rows):
    r = rows[(n, d)]
    k = '%d,%d' % (n, d)
    live = B.get(k, {}).get('upper'); low = B.get(k, {}).get('lower'); team = B.get(k, {}).get('team_upper')
    if live is None or low is None or (live == low and n < 17):
        continue
    v = [fl(r.get('delsarte')), fl(r.get('schrijver')), fl(r.get('lin')), r.get('parity_bound')]
    sv0 = SELF.get((n, d), {}).get('v4_self')
    # ours = premise-free (V1, V2) or self-consistent (V4 self); table-premised V3/V4 are shown for information only
    ours = min(x for x in (v[0], v[1], sv0) if x is not None) if sv0 is not None else min(x for x in v if x is not None)
    src = '; '.join('%s %s' % (t, val) for val, t in SRC.get((n, d), []) if val == live) or \
          ('(classical: BBMOS 1978 / AVZ 2001)' if live != low else 'exact (lower = upper)')
    if ours < live:
        verdict = 'BEATS live (by %d)' % (live - ours)
    elif ours == live:
        verdict = 'OK (reproduced)'
    else:
        verdict = 'GAP +%d' % (ours - live)
        gaps.append((n, d, live, ours, src))
    sv = SELF.get((n, d), {}).get('v4_self', '')
    out.append('| %d | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (n, d, live, team or '', src, v[0], v[1], v[2], v[3], sv, ours, verdict))
out += ['', '## Entries where the published value is NOT reproduced by LP / Schrijver / classical rows / parity', '']
for n, d, live, ours, src in gaps:
    out.append('- A(%d,%d): live %d, ours %d (+%d). Source: %s' % (n, d, live, ours, ours - live, src))
open(os.path.join(BASE, 'agents/sdp/AUDIT_PUBLISHED.md'), 'w').write('\n'.join(out) + '\n')
print('\n'.join(out[-len(gaps) - 1:]))
