import split
off=dict(best18=False,cw22=False,mel19=False,mel21=False,best20=False,kt23=False)
for n,d,n1 in [(18,4,2),(17,6,3),(20,8,4)]:
  for nm,o in [('schr only',dict(off,split=False)),('split only',dict(off,schrijver=False)),('both',dict(off)),('split only nolin',dict(off,schrijver=False,lin_basic=False)),('schr only nolin',dict(off,split=False,lin_basic=False))]:
    P=split.build(n,d,n1,opts=o,name='st_%d_%d_%d_%s'%(n,d,n1,nm.replace(' ','_')))
    info=P.solve('agents/sdp/runs/st',threads=4)
    print(n,d,n1,nm,info['phase.value'],round(info['bound'],6),flush=True)
