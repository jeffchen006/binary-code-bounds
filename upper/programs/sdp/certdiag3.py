exec(open('certdiag2.py').read().split("c=pickle.load")[0])
import mpmath
mpmath.mp.dps=60
# raw residual for a few top variables using raw Y (scaled) and raw diag z, including box multipliers
tops=[('Q', 0, 1, 2, 4, 3, 1, 6, 4),('Q', 1, 0, 1, 5, 2, 2, 7, 3)]
diag=info['Y'][len(P.blocks)][1]
nl=len(P.lin)
for v in tops:
    s=mpmath.mpf(P.c.get(v,0).numerator)/P.c.get(v,Fraction(1)).denominator*P.vscale[v]
    for bi,b in enumerate(P.blocks):
        t,M=info['Y'][bi]
        mat=b['ent'].get(v)
        if not mat: continue
        sc=b['scale']
        for (i,j),val in mat.items():
            e=mpmath.mpf(val.numerator)/val.denominator*sc[i]*sc[j]*P.vscale[v]
            Yij=M[i][j] if t=='m' else M[0]
            s+= e*(Yij if i==j else 2*Yij)
    # LP rows
    for l,(g0,co) in enumerate(P.lin):
        if v in co: s+= diag[l]*mpmath.mpf(co[v].numerator)/co[v].denominator*P.vscale[v]
    # box rows: lb then ub for each var in order
    k=P.vidx[v]
    s_lb = diag[nl+2*k]*P.vscale[v]; s_ub=-diag[nl+2*k+1]*P.vscale[v]
    print(v,'raw resid (scaled) without box',mpmath.nstr(s,8),' box lb',mpmath.nstr(s_lb,8),' ub',mpmath.nstr(s_ub,8),' total',mpmath.nstr(s+s_lb+s_ub,8))
