# strength test: full 2^n x 2^n matrices R, R' (split Terwilliger) vs. my blocks; eigenvalues must coincide (as sets)
import numpy as np, random, itertools, split, sys
n=int(sys.argv[1]); n1=int(sys.argv[2]); n2=n-n1
P=split.build(n,1,n1,even=False,opts=dict(lin_basic=False,best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False))
rng=np.random.default_rng(0)
val={v:rng.normal() for v in P.vars}
pc=lambda z: bin(z).count('1')
m1=(1<<n1)-1; m2=((1<<n)-1)^m1
def key(Y,Z):
    a,b=Y,Z
    if n1>0: prs=((pc(a&m1),pc(a&m2)),(pc(b&m1),pc(b&m2)),(pc((a^b)&m1),pc((a^b)&m2)))
    else: prs=((pc(a),),(pc(b),),(pc(a^b),))
    return tuple(sorted(prs))
N=1<<n
R=np.zeros((N,N)); R2=np.zeros((N,N))
for Y in range(N):
    for Z in range(N):
        k=key(Y,Z)
        x=1.0 if all(sum(p)==0 for p in k) else val[k]
        R[Y,Z]=x
        k0=key(Y^Z,0)
        x0=1.0 if all(sum(p)==0 for p in k0) else val[k0]
        R2[Y,Z]=x0-x
for b in P.blocks:
    pass
# blocks: first the Schrijver blocks (2 per k), then split blocks; test split ones (and Schrijver ones vs full-symmetrized?)
ev=[]
for b in P.blocks:
    A=np.zeros((b['size'],)*2)
    for v,mat in b['ent'].items():
        x=1.0 if v is None else val[v]
        for (i,j),w in mat.items():
            e=float(w)*x*float(b['scale'][i]*b['scale'][j]); A[i,j]+=e
            if i!=j: A[j,i]+=e
    ev.append(np.linalg.eigvalsh(A))
nS=2*(n//2+1)
schr=ev[:nS]; spl=ev[nS:]
full=np.sort(np.concatenate([np.linalg.eigvalsh(R),np.linalg.eigvalsh(R2)]))
bl=np.sort(np.concatenate(spl)) if n1>0 else np.sort(np.concatenate(schr))
# every block eigenvalue must occur among full eigenvalues, and min eig must agree
def contained(a,b):
    return all(np.min(np.abs(b-x))<1e-8*max(1,abs(x)) for x in a)
print('full min',full.min(),'block min',bl.min(),'block eigs in full:',contained(bl,full),'full eigs in blocks:',contained(full,bl))
miss=[x for x in full if np.min(np.abs(bl-x))>1e-8*max(1,abs(x))]
print(len(miss), miss[:10], len(full), len(bl))
fR=np.linalg.eigvalsh(R); fR2=np.linalg.eigvalsh(R2)
blo=np.sort(np.concatenate(spl[0::2] if n1>0 else schr[0::2])); bl2=np.sort(np.concatenate(spl[1::2] if n1>0 else schr[1::2]))
print('R ok', contained(fR,blo), contained(blo,fR), "R' ok", contained(fR2,bl2), contained(bl2,fR2))
