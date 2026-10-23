import parity
off=dict(best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False)
for w in range(0,23,2):
    for nm,o in [('basic',off),('all',{})]:
        P=parity.build(23,10,0,150,'2mod4',w,opts=o,name='s2310_%d_%s'%(w,nm))
        info=P.solve('agents/sdp/runs/s2310',threads=1)
        print(w,nm,info['phase.value'],round(info['bound'],5),'EXCL' if info['bound']<150 else '',flush=True)
