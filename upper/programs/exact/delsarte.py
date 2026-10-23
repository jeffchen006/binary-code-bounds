import sys
from math import comb
from scipy.optimize import linprog
def K(k, x, n): return sum((-1)**j*comb(x,j)*comb(n-x,k-j) for j in range(k+1))
def lp(n,d):
    js=[j for j in range(d,n+1) if j%2==0]
    A_ub=[[-K(k,j,n) for j in js] for k in range(1,n+1)]; b_ub=[K(k,0,n) for k in range(1,n+1)]
    r=linprog([-1]*len(js),A_ub=A_ub,b_ub=b_ub,bounds=[(0,None)]*len(js),method='highs')
    return 1-r.fun, dict((j,round(a,3)) for j,a in zip(js,r.x) if a>1e-9)
for n,d in [(17,8),(21,10),(25,12),(13,6),(9,4),(22,10),(26,12)]: print(n,d,lp(n,d))
