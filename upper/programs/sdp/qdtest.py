import split, sdpcore
P=split.build(12,4,0,name='qdt')
info=P.solve('agents/sdp/runs/qdt',threads=1)
c=P.make_certificate(info,K=100); b=P.check_certificate(c)
print(info['phase.value'],info['bound'],float(b))
