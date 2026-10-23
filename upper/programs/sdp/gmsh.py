import gms0, sys, time, sdpcore
n,d=int(sys.argv[1]),int(sys.argv[2]); solver=sys.argv[3] if len(sys.argv)>3 else 'csdp'
off=dict(best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False)
t=time.time()
P=gms0.build_gms(n,d,opts=off,name='gh_%d_%d'%(n,d))
if solver=='csdp':
    info=sdpcore.solve_csdp(P,'agents/sdp/runs/gmsh',threads=4,cpus=None)
else:
    info=P.solve('agents/sdp/runs/gmsh',threads=4)
print(n,d,solver,info['phase.value'],'bound',info['bound'],'pbound',info.get('pbound'),'%.0fs'%(time.time()-t),flush=True)
