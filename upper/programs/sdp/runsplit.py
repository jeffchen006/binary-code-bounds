import split, sys, time, json, os
# usage: runsplit.py tag threads "n,d,n1[,opt=val...]" ...
tag=sys.argv[1]; th=int(sys.argv[2])
S='agents/sdp/runs/'+tag
os.makedirs(S,exist_ok=True)
for spec in sys.argv[3:]:
    p=spec.split(','); n,d,n1=map(int,p[:3]); opts={}
    for kv in p[3:]:
        k,v=kv.split('='); opts[k]=bool(int(v))
    t0=time.time()
    P=split.build(n,d,n1,opts=opts,name='s_%d_%d_%d'%(n,d,n1)+''.join('_%s%d'%(k,v) for k,v in opts.items()))
    nv=len(P.vars); mb=max(b['size'] for b in P.blocks)
    info=P.solve(S,threads=th)
    line='%s n=%d d=%d n1=%d opts=%s vars=%d maxblk=%d %s bound=%.6f pbound=%.6f %.0fs'%(tag,n,d,n1,opts,nv,mb,info['phase.value'],info['bound'],info.get('pbound',0),time.time()-t0)
    print(line,flush=True)
    open('agents/sdp/results_all.txt','a').write(line+'\n')
