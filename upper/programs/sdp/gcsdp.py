import sys, time, mpmath, gms0, sdpcore
n,d,T=map(int,sys.argv[1:4])
off=dict(best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False)
t=time.time()
P=gms0.build_full(n,d,1,opts=off,name='gc_%d_%d'%(n,d)); P=gms0.homogenize(P,n)
P.vscale={v: mpmath.mpf(T)/2**n for v in P.vars}; sdpcore.equilibrate(P)
info=sdpcore.solve_csdp(P,'agents/sdp/runs/gcsdp',threads=3,cpus=None)
c=P.make_certificate(info,K=60); b=P.check_certificate(c)
print(n,d,info['phase.value'],'num',info['bound'],'pnum',info['pbound'],'rig',float(b),'%.0fs'%(time.time()-t),flush=True)
