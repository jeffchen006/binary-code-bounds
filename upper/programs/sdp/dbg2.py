import split, check_sdp_cert as C
from fractions import Fraction
n,d=24,6
opts=dict(lin_basic=True,best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False)
P=split.build(n,d,0,opts=opts)
def allowed(a): return a==0 or (a>=d and a%2==0)
def x(i,j,t):
    tri=tuple(sorted((i,j,i+j-2*t)))
    if not all(allowed(a) for a in tri): return {}
    if tri==(0,0,0): return {None:1}
    return {tri:1}
for b in P.blocks:
    k,typ,rows=b['meta']['k'],b['meta']['typ'],b['meta']['rows']
    for p,i in enumerate(rows):
        for q,j in enumerate(rows):
            if q<p: continue
            ent={}
            for t in range(max(0,i+j-n),min(i,j)+1):
                be=C.beta(n,i,j,k,t)
                if typ==0:
                    for kk,v in x(i,j,t).items(): ent[kk]=ent.get(kk,0)+be*v
                else:
                    for kk,v in x(i+j-2*t,0,0).items(): ent[kk]=ent.get(kk,0)+be*v
                    for kk,v in x(i,j,t).items(): ent[kk]=ent.get(kk,0)-be*v
            mine={}
            for v,mat in b['ent'].items():
                if (p,q) in mat:
                    mine[None if v is None else tuple(s[0] for s in v)]=mat[(p,q)]
            ent={a:c for a,c in ent.items() if c!=0}
            if mine!=ent: print('DIFF k',k,'typ',typ,'i,j',i,j,mine,ent)
print('done')
