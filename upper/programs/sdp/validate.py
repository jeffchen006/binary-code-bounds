# check that actual codes satisfy every constraint of the built Problem
import random, math, itertools, numpy as np, sys
from fractions import Fraction
import schrijver
def code_x(C, n):
    M=len(C); lam={}
    for x in C:
        for y in C:
            a=x^y
            for z in C:
                b=x^z
                key=(bin(a).count('1'),bin(b).count('1'),bin(a&b).count('1'))
                lam[key]=lam.get(key,0)+1
    xs={}
    for (i,j,t),L in lam.items():
        g=math.factorial(n)//(math.factorial(t)*math.factorial(i-t)*math.factorial(j-t)*math.factorial(n-i-j+t))
        xs[(i,j,t)]=Fraction(L, M*g)
    vals={}
    for (i,j,t),v in xs.items():
        k=schrijver.tri(i,j,t)
        if k in vals: assert vals[k]==v,(k,vals[k],v)
        vals[k]=v
    return vals, M
def check(P, vals, M):
    val=lambda v: vals.get(v,Fraction(0))
    for v in vals:
        if v!=(0,0,0): assert v in P.vidx, ('nonzero var not in problem',v)
    obj=P.c0+sum(c*val(v) for v,c in P.c.items()); assert obj==M,(obj,M)
    for g0,co in P.lin:
        assert g0+sum(a*val(v) for v,a in co.items())>=0
    mine=1e9
    for b in P.blocks:
        A=np.zeros((b['size'],)*2)
        for v,mat in b['ent'].items():
            x=1 if v is None else float(val(v))
            for (i,j),w in mat.items():
                A[i,j]+=float(w)*x*float(b['scale'][i]*b['scale'][j])
                if i!=j: A[j,i]+=float(w)*x*float(b['scale'][i]*b['scale'][j])
        mine=min(mine,np.linalg.eigvalsh(A).min()/max(1,abs(A).max()))
    return mine
random.seed(1)
for trial in range(40):
    n=random.randint(6,10); d=random.choice([3,4,5,6]); even=(d%2==0) and random.random()<0.7
    words=list(range(2**n)); random.shuffle(words)
    C=[]
    for w in words:
        if even and bin(w).count('1')%2: continue
        if all(bin(w^c).count('1')>=d for c in C): C.append(w)
        if random.random()<0.02: break
    vals,M=code_x(C,n)
    P=schrijver.build(n,d,even=even)
    print(n,d,even,M,'min rel eig %.2e'%check(P,vals,M))
