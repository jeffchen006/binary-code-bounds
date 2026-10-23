import gms0, parity, sdpcore, mpmath, pickle, sys
from fractions import Fraction
import numpy as np
n,d,T=21,10,47
off=dict(best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False)
name='gp_21_10_47_odd_None_l0_g1qd'
P=gms0.build_full(n,d,1,opts=off,name=name)
Q=parity.build(n,d,0,T,'odd',None,opts=dict(off,lin_basic=False)); P.lin=P.lin+Q.lin
P=gms0.homogenize(P,n)
P.vscale={v: mpmath.mpf(T)/2**n for v in P.vars}; sdpcore.equilibrate(P)
P.write_sdpa('/tmp/work')
info=P.read_result('runs/gmspar/%s.out'%name)
negs=[]
for bi,b in enumerate(P.blocks):
    t,M=info['Y'][bi]
    A=np.array([[float(x) for x in row] for row in M]) if t=='m' else np.array([[float(M[0])]])
    ev=np.linalg.eigvalsh((A+A.T)/2)
    negs.append((ev.min()/max(1e-300,abs(ev).max()), ev.min(), ev.max(), b['meta'].get('kind'), b['meta'].get('m'),b['meta'].get('k'),b['meta'].get('l')))
negs.sort()
for x in negs[:6]: print('rel min eig %.2e min %.2e max %.2e'%x[:3], x[3:])
c=pickle.load(open('runs/gmspar/%s.cert.pkl'%name,'rb'))['cert']
K0=P.c0; r={v:P.c.get(v,Fraction(0)) for v in P.vars}
for bi,bl in enumerate(P.blocks):
    V=c['V'][bi]; m=bl['size']
    Y=[[sum((V[i][s]*V[j][s] for s in range(len(V[i]))),Fraction(0)) for j in range(m)] for i in range(m)]
    for v,mat in bl['ent'].items():
        tot=sum((val*(Y[i][j] if i==j else 2*Y[i][j]) for (i,j),val in mat.items()),Fraction(0))
        if v is None: K0+=tot
        else: r[v]+=tot
for z,(g0,co) in zip(c['z'],P.lin):
    K0+=z*g0
    for v,a in co.items(): r[v]+=z*a
pos=sorted(((float(x),v) for v,x in r.items() if x>0),reverse=True)
print('K0',float(K0),'sum pos r',sum(x for x,_ in pos),'top',pos[:6])
