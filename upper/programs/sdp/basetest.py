import os, sys
sys.dont_write_bytecode = True
import bounds
import split, sdpcore
n, d = int(sys.argv[1]), int(sys.argv[2])
ALL = ['best18', 'cw22', 'mel19', 'mel21', 'best20', 'kt23']
for on in sys.argv[3:]:
    ON = [] if on == 'none' else on.split(',')
    P = split.build(n, d, 0, opts={k: (k in ON) for k in ALL}, name='base_%d_%d_%s' % (n, d, on.replace(',', '_')))
    (sdpcore.equilibrate(P) if os.environ.get("EQ") == "1" else None); info = P.solve('agents/sdp/runs/base', threads=1)
    cert = P.make_certificate(info, K=80); rig = P.check_certificate(cert)
    print(n, d, on, info['phase.value'], info['bound'], 'rig %.6f' % float(rig), flush=True)
