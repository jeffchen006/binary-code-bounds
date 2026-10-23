import schrijver, sys, json
S='/tmp/work'
cur={4:{17:3276,18:6551,19:13104,20:26168,21:43688,22:87333,23:172361,24:344308,25:599184,26:1198368,27:2396736,28:4792950},
6:{17:340,18:673,19:1237,20:2279,21:4096,22:6941,23:13674,24:24106,25:47538,26:84260,27:157285,28:291269},
8:{25:5421,26:9275,27:17099,28:32151},10:{21:47,22:84,23:150,24:268,25:466,26:836,27:1585,28:2817},12:{25:55,26:96,27:169,28:288}}
out=open('results_schrijver.txt','a')
for d in (4,6,8,10,12):
  for n,u in cur[d].items():
    P=schrijver.build(n,d); info=P.solve(S)
    line='%d %d cur %d schrijver %s %.6f'%(n,d,u,info['phase.value'],info['bound'])
    print(line,flush=True); out.write(line+'\n'); out.flush()
