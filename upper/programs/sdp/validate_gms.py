import random, math, itertools, sys
from fractions import Fraction
from validate2 import code_x, check
import gms
def pats(ws, n):
    c=[0]*8
    w0=ws[0]
    for b in range(n):
        bits=tuple((w>>b)&1 for w in ws)
        c[gms._pat_index(bits)]+=1
    return tuple(c)
def gms_vals(C,n):
    vals,M=code_x(C,n,0)
    cnt={}
    for S in itertools.combinations(C,4):
        T=gms.canon4(pats(S,n))
        cnt[T]=cnt.get(T,0)+1
    for T,k in cnt.items():
        orbitvecs=set()
        for tab in gms._PT:
            v=[0]*8
            for q in range(8): v[tab[q]]+=T[q]
            orbitvecs.add(tuple(v))
        mult=math.factorial(n)
        for x in T: mult//=math.factorial(x)
        orb=Fraction(2**(n-1)*mult*len(orbitvecs),24)
        vals[('Q',)+T]=Fraction(2**(n-1)*k)/(M*orb)
    return vals,M
if __name__=='__main__':
 random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 7)
 for trial in range(12):
    n=random.randint(7,10); d=random.choice([2,4,4,6])
    words=[w for w in range(2**n) if bin(w).count('1')%2==0]; random.shuffle(words)
    C=[]; lim=random.randint(4,40)
    for w in words:
        if all(bin(w^c).count('1')>=d for c in C): C.append(w)
        if len(C)>=lim: break
    vals,M=gms_vals(C,n)
    P=gms.build(n,d)
    print(n,d,M,'min rel eig %.2e'%check(P,vals,M),flush=True)
