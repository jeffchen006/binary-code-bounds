import parity, sys, time
n,d,n1,M=map(int,sys.argv[1:5]); case=sys.argv[5]; w=None if sys.argv[6]=='None' else int(sys.argv[6])
t=time.time()
P=parity.build(n,d,n1,M,case,w)
print('vars',len(P.vars),'blocks',len(P.blocks),'maxblk',max(b['size'] for b in P.blocks),'lin',len(P.lin),'%.0fs'%(time.time()-t),flush=True)
info=P.solve('agents/sdp/runs/ptest2',threads=4)
print(n,d,n1,M,case,w,info['phase.value'],info['bound'],'EXCL' if info['bound']<M else 'feasible','%.0fs'%(time.time()-t),flush=True)
