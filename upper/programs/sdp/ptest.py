import parity, sys, time
n,d,M,w=map(int,sys.argv[1:5])
t=time.time()
P=parity.build_split2mod4(n,d,M,w)
print('vars',len(P.vars),'blocks',len(P.blocks),'maxblk',max(b['size'] for b in P.blocks),'lin',len(P.lin),'%.0fs'%(time.time()-t),flush=True)
info=P.solve('agents/sdp/runs/psplit',threads=4)
print(n,d,M,w,info['phase.value'],info['bound'],'EXCL' if info['bound']<M else 'feasible','%.0fs'%(time.time()-t))
