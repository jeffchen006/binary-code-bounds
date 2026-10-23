exec(open('dbg3.py').read().split('# block-only')[0])
target=(10,14,14)
n,dd=24,6
def allowed(a): return a==0 or (a>=dd and a%2==0)
def x(i,j,t):
    tri=tuple(sorted((i,j,i+j-2*t)))
    if not all(allowed(a) for a in tri): return {}
    if tri==(0,0,0): return {None:1}
    return {tri:1}
for bi,b in enumerate(P.blocks):
    V=cert['V'][bi]; m=b['size']
    Y=[[sum((V[i][s]*V[j][s] for s in range(len(V[i]))),Fraction(0)) for j in range(m)] for i in range(m)]
    a=0
    for v,mat in b['ent'].items():
        if v is None or tuple(p[0] for p in v)!=target: continue
        a=sum((val*(Y[i][j] if i==j else 2*Y[i][j]) for (i,j),val in mat.items()),Fraction(0))
    jb=J['blocks'][bi]; k,typ,rows=jb['k'],jb['typ'],jb['rows']
    c=0
    for p in range(m):
        for q in range(m):
            i,j=rows[p],rows[q]
            for t in range(max(0,i+j-n),min(i,j)+1):
                be=C.beta(n,i,j,k,t)
                forms=[(x(i,j,t),be)] if typ==0 else [(x(i+j-2*t,0,0),be),(x(i,j,t),-be)]
                for f_,cc in forms:
                    c+=cc*f_.get(target,0)*Y[p][q]
    if a!=c: print(bi,k,typ,float(a),float(c))
b=P.blocks[12]
for v,mat in b['ent'].items():
    lower=[k for k in mat if k[0]>k[1]]
    if lower: print('lower entries', v, lower[:5])
bi=12; V=cert['V'][bi]; m=P.blocks[bi]['size']
Y=[[sum((V[i][s]*V[j][s] for s in range(len(V[i]))),Fraction(0)) for j in range(m)] for i in range(m)]
asym=[(i,j) for i in range(m) for j in range(m) if Y[i][j]!=Y[j][i]]
print('asym',asym[:5], 'size',m, 'Vrowlens',set(len(r) for r in V))
b=P.blocks[12]; jb=J['blocks'][12]; k,typ,rows=jb['k'],jb['typ'],jb['rows']
for v,mat in b['ent'].items():
    if tuple(p[0] for p in v)!=target: continue
    for (p,q),val in mat.items():
        i,j=rows[p],rows[q]
        cc=0
        for t in range(max(0,i+j-n),min(i,j)+1):
            be=C.beta(n,i,j,k,t)
            cc+=be*x(i,j,t).get(target,0)
        cc2=0
        for t in range(max(0,i+j-n),min(i,j)+1):
            be=C.beta(n,j,i,k,t)
            cc2+=be*x(j,i,t).get(target,0)
        print(p,q,i,j,val,cc,cc2)
