import gms0, mpmath, math
n,d,T=21,10,47
off=dict(best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False)
P=gms0.build_full(n,d,1,opts=off,name='ct')
P=gms0.homogenize(P,n)
vs=float(T)/2**n
worst=[]
for bi,b in enumerate(P.blocks):
    m=b['size']; sc=b['scale']
    diag=[0.0]*m
    for v,mat in b['ent'].items():
        for (i,j),val in mat.items():
            if i==j:
                diag[i]+=abs(float(val))*float(sc[i])**2*(1.0 if v is None else vs)
    nz=[x for x in diag if x>0]
    if nz:
        worst.append((max(nz)/min(nz), b['meta'].get('kind'), b['meta'].get('m'), b['meta'].get('k'), b['meta'].get('l'), m, min(nz), max(nz)))
worst.sort(reverse=True)
for w in worst[:12]: print('%.2e'%w[0], w[1:5], 'size',w[5],'min %.2e max %.2e'%(w[6],w[7]))
import statistics
print('overall diag range', min(w[6] for w in worst), max(w[7] for w in worst))
