import gms0, sys, time
for n,d in [(21,10),(25,12),(18,6)]:
    t=time.time()
    P=gms0.build_full(n,d,1,opts=dict(best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False),name='sz')
    print(n,d,'vars',len(P.vars),'blocks',len(P.blocks),'maxblk',max(b['size'] for b in P.blocks),'nnz',sum(len(m) for b in P.blocks for m in b['ent'].values()),'%.0fs'%(time.time()-t),flush=True)
