import gms0, time
t=time.time()
import sys
n,d=int(sys.argv[1]),int(sys.argv[2])
P=gms0.build_full(n,d,1,opts=dict(best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False),name='sz')
print('vars',len(P.vars),'blocks',len(P.blocks),'maxblk',max(b['size'] for b in P.blocks),'sum blk^2',sum(b['size']**2 for b in P.blocks),'%.0fs'%(time.time()-t))
