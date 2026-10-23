# fast double-precision scan: build problem, write SDPA, solve with sdpap (SDPA double)
import split, sys, time, os, io, contextlib
import sdpap
from sdpap import fileio
S='agents/sdp/runs/dscan'
os.makedirs(S,exist_ok=True)
def dsolve(P, th=4):
    f=os.path.join(S,P.name+'.dat-s'); P.write_sdpa(f)
    A,b,c,K,J=fileio.importsdpa(f)
    buf=io.StringIO()
    with contextlib.redirect_stdout(buf):
        x,y,info,tinfo,sinfo=sdpap.solve(A,b,c,K,J,{'print':'no','epsilonStar':1e-10,'epsilonDash':1e-10,'lambdaStar':1e4,'numThreads':th})
    return float(P.c0)-sinfo['dualObj'], float(P.c0)-sinfo['primalObj'], info['phasevalue']
if __name__=='__main__':
    for spec in sys.argv[1:]:
        p=spec.split(','); n,d,n1=map(int,p[:3]); opts={}
        for kv in p[3:]:
            k,v=kv.split('='); opts[k]=bool(int(v))
        t=time.time()
        P=split.build(n,d,n1,opts=opts,name='d_%d_%d_%d'%(n,d,n1)+''.join('_%s%d'%(k,int(v)) for k,v in opts.items()))
        db,pb,ph=dsolve(P)
        line='dscan n=%d d=%d n1=%d opts=%s vars=%d dual_bound=%.5f primal=%.5f %s %.0fs'%(n,d,n1,opts,len(P.vars),db,pb,ph,time.time()-t)
        print(line,flush=True)
        open('agents/sdp/results_all.txt','a').write(line+'\n')
