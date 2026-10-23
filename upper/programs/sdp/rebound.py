import gms0, parity, sdpcore, mpmath, pickle, sys
from fractions import Fraction
n,d,T=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]); case=sys.argv[4]; w=None if sys.argv[5]=='None' else int(sys.argv[5]); name=sys.argv[6]
off=dict(best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False)
P=gms0.build_full(n,d,1,opts=off,name=name)
if case!='base':
    Q=parity.build(n,d,0,T,case,w,opts=dict(off,lin_basic=False)); P.lin=P.lin+Q.lin
P=gms0.homogenize(P,n)
c=pickle.load(open('runs/gmspar/%s.cert.pkl'%name,'rb'))['cert']
print('ub=1    :',float(P.check_certificate(c)))
for v in P.vars: P.ub[v]=Fraction(T,2**n)
print('ub=T/2^n:',float(P.check_certificate(c)))
