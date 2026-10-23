"""Run one spec: n,d,n1[,opt=0/1...]; solve with SDPA-DD, rigorous cert bound, log."""
import split, sys, time, os, json, pickle
from fractions import Fraction
BASE='agents/sdp'
def run(spec, tag, threads=4, K=60):
    p=spec.split(','); n,d,n1=map(int,p[:3]); opts={}
    for kv in p[3:]:
        k,v=kv.split('='); opts[k]=bool(int(v))
    name='j_%d_%d_%d'%(n,d,n1)+''.join('_%s%d'%(k,int(v)) for k,v in sorted(opts.items()))
    t0=time.time()
    P=split.build(n,d,n1,opts=opts,name=name)
    tb=time.time()-t0
    info=P.solve(os.path.join(BASE,'runs',tag),threads=threads)
    ts=time.time()-t0-tb
    rig=None
    try:
        cert=P.make_certificate(info,K=K)
        b=P.check_certificate(cert)
        rig=b
        with open(os.path.join(BASE,'runs',tag,name+'.cert.pkl'),'wb') as f:
            pickle.dump({'spec':spec,'cert':cert,'bound':b},f)
    except Exception as e:
        print('cert failed',repr(e),flush=True)
    line='%s | %s | n=%d d=%d n1=%d opts=%s vars=%d blocks=%d maxblk=%d lin=%d | %s num=%.6f rig=%s floor=%s | build %.0fs solve %.0fs total %.0fs'%(
        time.strftime('%m-%d %H:%M'),tag,n,d,n1,opts,len(P.vars),len(P.blocks),max(b['size'] for b in P.blocks),len(P.lin),
        info['phase.value'],info['bound'] or -1,('%.6f'%float(rig)) if rig is not None else 'NA',(rig.numerator//rig.denominator) if rig is not None else 'NA',tb,ts,time.time()-t0)
    print(line,flush=True)
    open(os.path.join(BASE,'results_all.txt'),'a').write(line+'\n')
if __name__=='__main__':
    tag=sys.argv[1]; th=int(sys.argv[2])
    for spec in sys.argv[3:]:
        run(spec,tag,th)
