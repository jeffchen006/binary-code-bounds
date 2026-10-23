exec(open('dbg.py').read().split('# per-component comparison')[0].replace("print(float(P.check_certificate(cert)))","").replace("pr,b=C.check(f); print(float(b))",""))
from fractions import Fraction
# lin-only contributions, mine
r_lin={}
for z,(g0,co) in zip(cert['z'],P.lin):
    for v,a in co.items():
        kk=tuple(p[0] for p in v); r_lin[kk]=r_lin.get(kk,0)+z*a
# checker lin-only
r_lin2={}
D1=Fraction(1,2**J['K'])
for row in J['lin']:
    for key,v in row['coef'].items():
        kk=tuple(int(s) for s in key.split(',')); r_lin2[kk]=r_lin2.get(kk,0)+row['z']*D1*Fraction(v)
bad=[k for k in set(r_lin)|set(r_lin2) if r_lin.get(k,0)!=r_lin2.get(k,0)]
print('lin diffs',len(bad))
# block-only
r_b={}
for bi,b in enumerate(P.blocks):
    V=cert['V'][bi]; m=b['size']
    Y=[[sum((V[i][s]*V[j][s] for s in range(len(V[i]))),Fraction(0)) for j in range(m)] for i in range(m)]
    for v,mat in b['ent'].items():
        if v is None: continue
        tot=sum((val*(Y[i][j] if i==j else 2*Y[i][j]) for (i,j),val in mat.items()),Fraction(0))
        kk=tuple(p[0] for p in v); r_b[kk]=r_b.get(kk,0)+tot
r_b2={}
n,dd=24,6
def allowed(a): return a==0 or (a>=dd and a%2==0)
def x(i,j,t):
    tri=tuple(sorted((i,j,i+j-2*t)))
    if not all(allowed(a) for a in tri): return {}
    if tri==(0,0,0): return {None:1}
    return {tri:1}
for jb in J['blocks']:
    k,typ,rows,VV=jb['k'],jb['typ'],jb['rows'],jb['V']
    m=len(rows)
    Y=[[sum(a*b for a,b in zip(VV[p],VV[q])) for q in range(m)] for p in range(m)]
    for p in range(m):
        for q in range(m):
            i,j=rows[p],rows[q]
            for t in range(max(0,i+j-n),min(i,j)+1):
                be=C.beta(n,i,j,k,t)
                forms=[(x(i,j,t),be)] if typ==0 else [(x(i+j-2*t,0,0),be),(x(i,j,t),-be)]
                for f_,c in forms:
                    for kk,v in f_.items():
                        if kk is None: continue
                        r_b2[kk]=r_b2.get(kk,0)+c*v*Fraction(Y[p][q],4**J['K'])
bad=[k for k in set(r_b)|set(r_b2) if r_b.get(k,0)!=r_b2.get(k,0)]
print('block diffs',len(bad), bad[:3], [ (float(r_b.get(k,0)),float(r_b2.get(k,0))) for k in bad[:3]])
# compare V
for bi,jb in enumerate(J['blocks']):
    V=cert['V'][bi]
    if any(Fraction(jb['V'][i][s],2**J['K'])!=V[i][s] for i in range(len(V)) for s in range(len(V[i]))): print('V mismatch',bi)
