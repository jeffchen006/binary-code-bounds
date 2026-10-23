import numpy as np, sys, itertools, gms0, gms, split
from validate_gms import pats
n=int(sys.argv[1]); d=int(sys.argv[2]); T=int(sys.argv[3])
off=dict(best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False,lin_basic=False)
P=split.build(n,d,0,opts=off,name='bc')
def okd(a): return a==0 or (a>=d and a%2==0)
def key3(a,b,c):
    k=tuple(sorted((a,b,c)))
    if not all(okd(x) for x in k): return None
    if k==(0,0,0): return 'one'
    kk=tuple((x,) for x in k); return kk if kk in P.vidx else P.var(kk)
nb=len(P.blocks)
gms0.add_mempty(P,n,d,T,key3,None)
rng=np.random.default_rng(3)
val={v:rng.normal() for v in P.vars}
pc=lambda z: bin(z).count('1')
def X(ws):
    ws=list(dict.fromkeys(ws))
    for a,b in itertools.combinations(ws,2):
        if not okd(pc(a^b)): return 0.0
    if len(ws)==1: return 1.0
    if len(ws)==2:
        a=pc(ws[0]^ws[1]); return val[tuple((x,) for x in sorted((0,a,a)))]
    if len(ws)==3:
        a,b,c=ws; k=tuple(sorted((pc(a^b),pc(a^c),pc(b^c)))); return val.get(tuple((x,) for x in k),0.0)
    w0=ws[0]; k=('Q',)+gms.canon4(pats([w^w0 for w in ws],n)); return val.get(k,0.0)
N=2**n
rows=[()]+[(v,) for v in range(N)]+[(v,w) for v in range(N) for w in range(v+1,N) if okd(pc(v^w))]
corner=2**n/T
F=np.zeros((len(rows),len(rows)))
for a,r in enumerate(rows):
    for b,s in enumerate(rows):
        S=list(r)+list(s)
        F[a,b]=corner if not S else X(S)
fe=np.linalg.eigvalsh(F)
ev=[]
for b in P.blocks[nb:]:
    A=np.zeros((b['size'],)*2)
    for v,mat in b['ent'].items():
        x=1.0 if v is None else val[v]
        for (i,j),w in mat.items():
            e=float(w)*x*float(b['scale'][i]*b['scale'][j]); A[i,j]+=e
            if i!=j: A[j,i]+=e
    ev.append(np.linalg.eigvalsh(A))
be=np.concatenate(ev)
def cont(a,b): return all(np.min(np.abs(b-x))<1e-7*max(1,abs(x)) for x in a)
print('n',n,'d',d,'rows',len(rows),'full min %.6f block min %.6f'%(fe.min(),be.min()),'blocks->full',cont(be,fe),'full->blocks',cont(fe,be),'#eig',len(fe),len(be))
