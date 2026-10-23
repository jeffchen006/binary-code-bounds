import numpy as np, sys, gms, itertools
from validate_gms import pats
n=int(sys.argv[1]); d=int(sys.argv[2]); m=int(sys.argv[3])
P=gms.build(n,d,ms=[m],quad_lin=False)
rng=np.random.default_rng(1)
val={v:rng.normal() for v in P.vars}
def okd(a): return a==0 or (a>=d and a%2==0)
pc=lambda z: bin(z).count('1')
u=(1<<m)-1
Np=[v for v in range(2**n) if pc(v)%2==0 and pc(v)>=d and okd(pc(u^v))]
def X(ws):
    ws=list(dict.fromkeys(ws))  # distinct preserving order
    for a,b in itertools.combinations(ws,2):
        if not okd(pc(a^b)): return 0.0
    if len(ws)==2:
        k=tuple(sorted((0,m,m))); return val[tuple((x,) for x in k)]
    if len(ws)==3:
        a,b,c=ws; k=tuple(sorted((pc(a^b),pc(a^c),pc(b^c)))); kk=tuple((x,) for x in k)
        return val.get(kk,0.0)
    k=('Q',)+gms.canon4(pats(ws,n)); return val.get(k,0.0)
F=np.array([[X([0,u,v,w]) for w in Np] for v in Np])
fe=np.linalg.eigvalsh(F)
ev=[]
for b in P.blocks:
    if b.get('meta',{}).get('kind')!='gms': continue
    A=np.zeros((b['size'],)*2)
    for v,mat in b['ent'].items():
        x=1.0 if v is None else val[v]
        for (i,j),w in mat.items():
            e=float(w)*x*float(b['scale'][i]*b['scale'][j]); A[i,j]+=e
            if i!=j: A[j,i]+=e
    ev.append(np.linalg.eigvalsh(A))
be=np.concatenate(ev)
ok1=all(np.min(np.abs(fe-x))<1e-8*max(1,abs(x)) for x in be)
ok2=all(np.min(np.abs(be-x))<1e-8*max(1,abs(x)) for x in fe)
print('n',n,'d',d,'m',m,'|N\'|',len(Np),'full min %.6f block min %.6f'%(fe.min(),be.min()),'blocks->full',ok1,'full->blocks',ok2)
