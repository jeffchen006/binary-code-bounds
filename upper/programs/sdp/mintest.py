import parity, sys
n,d,M=24,6,int(sys.argv[1]); case=sys.argv[2]; w=None if sys.argv[3]=='None' else int(sys.argv[3])
base=dict(lin_basic=True,best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False)
for extra in [{}, {'best18':True}, {'cw22':True}, {'mel19':True,'mel21':True,'best20':True}, {'kt23':True}]:
    o=dict(base); o.update(extra)
    P=parity.build(n,d,0,M,case,w,opts=o,name='min_%d_%s_%s_%s'%(M,case,w,'_'.join(extra)))
    info=P.solve('agents/sdp/runs/min',threads=4)
    print(M,case,w,extra,info['phase.value'],info['bound'],flush=True)
