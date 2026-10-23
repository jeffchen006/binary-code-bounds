"""Exclusion chain for A(n,d), d even: base SDP floor F, then try to exclude M=F,F-1,... with parity cases.
Every solve gets a rigorous rational certificate (bound < M needed)."""
import sys, time, os, pickle
import split, parity
from fractions import Fraction
BASE='agents/sdp'
def solve_cert(P, tag, K=60):
    info=P.solve(os.path.join(BASE,'runs',tag),threads=int(os.environ.get('TH','4')))
    try:
        cert=P.make_certificate(info,K=K); b=P.check_certificate(cert)
    except Exception as e:
        print('cert fail',repr(e),flush=True); cert=None; b=None
    return info,cert,b
def log(s):
    print(s,flush=True); open(os.path.join(BASE,'chain_results.txt'),'a').write(time.strftime('%m-%d %H:%M ')+s+'\n')
def run(n,d,n1=0,opts=None,tag='chain',maxsteps=6):
    os.makedirs(os.path.join(BASE,'runs',tag),exist_ok=True)
    P=split.build(n,d,n1,opts=opts,name='base_%d_%d_%d'%(n,d,n1))
    info,cert,b=solve_cert(P,tag)
    F=b.numerator//b.denominator
    log('%d %d n1=%d opts=%s base num=%.6f rig=%.6f floor=%d'%(n,d,n1,opts,info['bound'],float(b),F))
    certs={'base':(P.name,b)}
    pickle.dump(cert,open(os.path.join(BASE,'runs',tag,P.name+'.cert.pkl'),'wb'))
    M=F
    for step in range(maxsteps):
        allex=True; worst=None
        for case,w in parity.cases(M,n):
            Q=parity.build(n,d,n1,M,case,w,opts=opts)
            info,cert,b=solve_cert(Q,tag)
            ex = b is not None and b < M
            if cert is not None: pickle.dump(cert,open(os.path.join(BASE,'runs',tag,Q.name+'.cert.pkl'),'wb'))
            log('%d %d n1=%d M=%d %s w=%s num=%.6f rig=%s %s'%(n,d,n1,M,case,w,info['bound'] or -1,('%.6f'%float(b)) if b is not None else 'NA','EXCL' if ex else 'feasible'))
            if not ex:
                allex=False; break
        if not allex: break
        M-=1
    log('%d %d n1=%d RESULT: A <= %d (opts=%s)'%(n,d,n1,M,opts))
    return M
if __name__=='__main__':
    import bounds
    for kv in os.environ.get('OVR','').split():
        a,b,c=map(int,kv.split(','));  bounds.OVERRIDE[(a,b)]=c
    tag=sys.argv[1]
    for spec in sys.argv[2:]:
        p=spec.split(','); n,d,n1=map(int,p[:3]); opts={}
        for kv in p[3:]:
            k,v=kv.split('='); opts[k]=bool(int(v))
        run(n,d,n1,opts or None,tag)
