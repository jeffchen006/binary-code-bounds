import gms0, sdpcore, mpmath, pickle, sys
from fractions import Fraction
import numpy as np
n,d,T=17,8,39
off=dict(best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False)
P=gms0.build_full(n,d,1,opts=off,name='gp_17_8_39_base_None_l0_g1eq'); P=gms0.homogenize(P,n)
P.vscale={v: mpmath.mpf(T)/2**n for v in P.vars}; sdpcore.equilibrate(P)
P._diag=None
P.write_sdpa('/tmp/work')
info=P.read_result('runs/gmspar/gp_17_8_39_base_None_l0_g1eq.out')
# eigenvalue stats of scaled Y blocks
negs=[]
for bi,b in enumerate(P.blocks):
    M=np.array([[float(x) for x in row] for row in info['Y'][bi][1]]) if info['Y'][bi][0]=='m' else np.array([[float(info['Y'][bi][1][0])]])
    ev=np.linalg.eigvalsh((M+M.T)/2)
    negs.append((ev.min(), ev.max(), b['meta'].get('kind')))
negs.sort()
print('most negative eigs', negs[:5])
cert=P.make_certificate(info,K=80)
b=P.check_certificate(cert); print('rig',float(b))
# loss breakdown: compute r_v
K0=P.c0; r={v:P.c.get(v,Fraction(0)) for v in P.vars}
for bi,bl in enumerate(P.blocks):
    V=cert['V'][bi]; m=bl['size']
    Y=[[sum((V[i][s]*V[j][s] for s in range(len(V[i]))),Fraction(0)) for j in range(m)] for i in range(m)]
    for v,mat in bl['ent'].items():
        tot=sum((val*(Y[i][j] if i==j else 2*Y[i][j]) for (i,j),val in mat.items()),Fraction(0))
        if v is None: K0+=tot
        else: r[v]+=tot
for z,(g0,co) in zip(cert['z'],P.lin):
    K0+=z*g0
    for v,a in co.items(): r[v]+=z*a
pos=sorted(((float(x),v) for v,x in r.items() if x>0),reverse=True)
print('K0',float(K0),'sum pos r',sum(x for x,_ in pos),'top',pos[:5])
