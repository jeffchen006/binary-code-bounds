import schrijver
S='/tmp/work'
known={4:{9:20,10:40,11:72,12:144,13:256,14:512,15:1024,16:2048},6:{12:24,13:32,14:64,15:128,16:256},8:{16:32,17:36,18:64,19:128,20:256}}
for d in (4,6,8):
  for n,a in known[d].items():
    P=schrijver.build(n,d); info=P.solve(S)
    print(n,d,'known',a,'sdp',info['phase.value'],-float(info['objValDual']),flush=True)
