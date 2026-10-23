import sdpap, sys, time
from sdpap import fileio
f=sys.argv[1]
t=time.time()
A,b,c,K,J=fileio.importsdpa(f)
x,y,info,tinfo,sinfo=sdpap.solve(A,b,c,K,J,{'print':'display','epsilonStar':1e-9,'epsilonDash':1e-9,'lambdaStar':1e4,'numThreads':4})
print(info, sinfo.get('primalObj'), sinfo.get('dualObj'), time.time()-t)
