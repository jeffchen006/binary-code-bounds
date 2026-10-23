# check real codes satisfy all constraints of split.build (with bounds tables, n small)
import random, math, numpy as np, sys
from fractions import Fraction
import split
def code_x(C, n, n1):
    M=len(C); lam={}
    m1=(1<<n1)-1; m2=((1<<n)-1)^m1
    pc=lambda z: bin(z).count('1')
    for x in C:
        for y in C:
            a=x^y
            for z in C:
                b=x^z
                if n1>0:
                    key=tuple(sorted(((pc(a&m1),pc(a&m2)),(pc(b&m1),pc(b&m2)),(pc((a^b)&m1),pc((a^b)&m2)))))
                    cfg=(pc(a&m1),pc(b&m1),pc(a&b&m1),pc(a&m2),pc(b&m2),pc(a&b&m2))
                else:
                    key=tuple(sorted(((pc(a),),(pc(b),),(pc(a^b),))))
                    cfg=(pc(a),pc(b),pc(a&b),0,0,0)
                lam[(key,cfg)]=lam.get((key,cfg),0)+1
    vals={}
    n2=n-n1
    for (key,(i,j,t,i2,j2,t2)),L in lam.items():
        if n1>0: g=split.multi(n1,i-t,j-t,t)*split.multi(n2,i2-t2,j2-t2,t2)
        else: g=split.multi(n,i-t,j-t,t)
        v=Fraction(L,M*g)
        if key in vals: assert vals[key]==v,(key,vals[key],v)
        vals[key]=v
    return vals,M
def check(P, vals, M):
    val=lambda v: vals.get(v,Fraction(0))
    for v in vals:
        if v[0]!='Q' and not all(sum(p)==0 for p in v): assert v in P.vidx, ('nonzero var not in problem',v)
    obj=P.c0+sum(c*val(v) for v,c in P.c.items()); assert obj==M,(obj,M)
    for g0,co in P.lin:
        s=g0+sum(a*val(v) for v,a in co.items())
        assert s>=0,(g0,co,s)
    mine=1e9
    for b in P.blocks:
        A=np.zeros((b['size'],)*2)
        for v,mat in b['ent'].items():
            x=1 if v is None else float(val(v))
            for (i,j),w in mat.items():
                e=float(w)*x*float(b['scale'][i]*b['scale'][j])
                A[i,j]+=e
                if i!=j: A[j,i]+=e
        mine=min(mine,np.linalg.eigvalsh(A).min()/max(1,abs(A).max()))
    return mine
if __name__=='__main__':
  random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)
  for trial in range(25):
    n=random.randint(7,11); d=random.choice([4,4,6]); n1=random.randint(0,n//2)
    words=list(range(2**n)); random.shuffle(words)
    C=[]
    lim=random.choice([5,10,40,1000])
    for w in words:
        if bin(w).count('1')%2: continue
        if all(bin(w^c).count('1')>=d for c in C): C.append(w)
        if len(C)>=lim: break
    vals,M=code_x(C,n,n1)
    P=split.build(n,d,n1)
    print(n,d,n1,M,len(P.vars),'min rel eig %.2e'%check(P,vals,M),flush=True)
