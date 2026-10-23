import mpmath, re, sys
mpmath.mp.dps=60
from sdpcore import Problem
dat='runs/gmspar/gp_21_10_47_odd_None_l0_g1qd.dat-s'; out='runs/gmspar/gp_21_10_47_odd_None_l0_g1qd.out'
L=open(dat).read().split('\n')
m=int(L[1]); nb=int(L[2]); bs=[int(x) for x in L[3].split()]
c=[mpmath.mpf(x) for x in L[4].split()]
P=Problem('x'); info=P.read_result(out)
Y=info['Y']
targets=set(range(1,m+1))
res=[mpmath.mpf(0)]*(m+1)
for ln in L[5:]:
    p=ln.split()
    if len(p)!=5: continue
    k,b,i,j=map(int,p[:4]); val=mpmath.mpf(p[4])
    if k==0: continue
    t,M=Y[b-1]
    if t=='m': y=M[i-1][j-1]
    else: y=M[i-1] if i==j else 0
    res[k]+= val*(y if i==j else 2*y)
worst=sorted(((abs(res[k]-c[k-1]),k) for k in range(1,m+1)),reverse=True)[:5]
print('max dual residual |<F_k,Y> - c_k|:',[(mpmath.nstr(a,5),k) for a,k in worst])
