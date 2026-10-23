import gms, split, parity, sdpcore, sys, time
n,d=int(sys.argv[1]),int(sys.argv[2]); kind=sys.argv[3]
off=dict(best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False)
t=time.time()
if kind=='schr': P=split.build(n,d,0,opts=off,name='cs_schr_%d_%d'%(n,d))
elif kind=='gms': P=gms.build(n,d,opts=off,name='cs_gms_%d_%d'%(n,d))
elif kind=='gmsall': P=gms.build(n,d,name='cs_gmsall_%d_%d'%(n,d))
info=sdpcore.solve_csdp(P,'agents/sdp/runs/csdp',threads=4,cpus=None)
c=P.make_certificate(info,K=60); b=P.check_certificate(c)
print(n,d,kind,info['phase.value'],'num',info['bound'],'pnum',info['pbound'],'rig',float(b),'%.0fs'%(time.time()-t),flush=True)
