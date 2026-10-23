import schrijver, sys
S='/tmp/work'
for (n,d,ev) in [(19,6,False),(18,5,False),(19,6,True),(20,8,True),(23,6,True),(22,10,True),(25,6,True),(17,4,True),(18,4,True)]:
    P=schrijver.build(n,d,even=ev)
    info=P.solve(S)
    print(n,d,ev,info['phase.value'],info['objValDual'],flush=True)
