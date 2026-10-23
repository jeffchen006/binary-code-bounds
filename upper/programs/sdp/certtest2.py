import split, sys, time
from fractions import Fraction
n,d,n1=map(int,sys.argv[1:4])
P=split.build(n,d,n1)
info=P.read_result(sys.argv[4])
t=time.time()
for K in (30,60):
    cert=P.make_certificate(info,K=K)
    b=P.check_certificate(cert)
    print('K',K,'rigorous bound',float(b),'floor',b.numerator//b.denominator,'%.1fs'%(time.time()-t),flush=True)
