import gms0, sys, time, sdpcore
n,d,T=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
off=dict(best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False)
t=time.time()
P=gms0.build_full(n,d,T,opts=off,light=True,name='gdd_%d_%d_%d'%(n,d,T))
info=P.solve('agents/sdp/runs/gms0dd',threads=4)
print(n,d,T,info['phase.value'],'bound',info['bound'],'pbound',info.get('pbound'),'%.0fs'%(time.time()-t),flush=True)
