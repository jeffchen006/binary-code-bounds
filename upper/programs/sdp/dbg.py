import json, split, check_sdp_cert as C
from fractions import Fraction
f='../../results/verified/UB_A_24_6_24105/base_n24_d6.json'
J=json.load(open(f))
opts=dict(lin_basic=True,best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False)
P=split.build(24,6,0,opts=opts)
print(len(P.blocks), len(J['blocks']), len(P.lin), len(J['lin']))
for b,jb in zip(P.blocks,J['blocks']):
    if b['meta']['rows']!=jb['rows'] or b['meta']['k']!=jb['k'] or b['meta']['typ']!=jb['typ']: print('meta mismatch',b['meta'],jb['k'],jb['typ'])
# rebuild cert dict
D=2**J['K']
cert={'V':[[[Fraction(x,D) for x in row] for row in jb['V']] for jb in J['blocks']]}
zs=[]
jl={ (Fraction(r['g0']), frozenset((tuple(int(s) for s in k.split(',')),Fraction(v)) for k,v in r['coef'].items())): Fraction(r['z'],D) for r in J['lin']}
for g0,co in P.lin:
    key=(g0, frozenset((tuple(p[0] for p in v),c) for v,c in co.items()))
    zs.append(jl.get(key,Fraction(0)))
print('matched z', sum(1 for z in zs if z), 'of', len(jl))
cert['z']=zs
print(float(P.check_certificate(cert)))
pr,b=C.check(f); print(float(b))
# per-component comparison
from fractions import Fraction
def mine(P,cert):
    K=P.c0; r={v:P.c.get(v,Fraction(0)) for v in P.vars}
    for bi,b in enumerate(P.blocks):
        V=cert['V'][bi]; m=b['size']
        Y=[[sum((V[i][s]*V[j][s] for s in range(len(V[i]))),Fraction(0)) for j in range(m)] for i in range(m)]
        for v,mat in b['ent'].items():
            tot=sum((val*(Y[i][j] if i==j else 2*Y[i][j]) for (i,j),val in mat.items()),Fraction(0))
            if v is None: K+=tot
            else: r[v]+=tot
    for z,(g0,co) in zip(cert['z'],P.lin):
        K+=z*g0
        for v,a in co.items(): r[v]+=z*a
    return K,{tuple(p[0] for p in v):x for v,x in r.items()}
K1,r1=mine(P,cert)
# checker internals: patch to capture L
import check_sdp_cert as CC
src=open('check_sdp_cert.py').read()
ns={}
exec(src.replace("    bound = L.pop(None, Fraction(0))","    globals()['LL']=dict(L)\n    bound = L.pop(None, Fraction(0))").replace("if __name__ == '__main__':","if False:"),ns)
ns['check'](f); L=ns['LL']
print('K', float(K1), float(L[None]))
for v in sorted(set(r1)|set(k for k in L if k is not None)):
    a=r1.get(v,0); b=L.get(v,0)
    if a!=b: print(v, float(a), float(b))
