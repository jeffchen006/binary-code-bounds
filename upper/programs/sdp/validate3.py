import validate2, split, itertools, random
from validate2 import code_x, check
def span(gens):
    S={0}
    for g in gens: S|={s^g for s in S}
    return sorted(S)
def bits(s): return int(s[::-1],2)
# extended Hamming [8,4,4]
H8=span([bits('11110000'),bits('00111100'),bits('00001111'),bits('10101010')])
# RM(1,4) [16,5,8]
RM=span([0xFFFF,0x00FF,0x0F0F,0x3333,0x5555])
# extended Hamming [16,11,4] shortened to length 12 (fix 4 coords zero)
def ham16():
    gens=[]
    cols=[c for c in range(16)]
    # parity check: rows = bits of index + all-ones
    words=[]
    import itertools
    basis=[]
    for w in range(1<<16):
        pass
    return None
cases=[(H8,8,4),(RM,16,8)]
# [12,7,4]: ext Hamming 16 shortened: take codewords of RM(2,4) (=[16,11,4]) zero on last 4 coords
RM2=span([0xFFFF,0x00FF,0x0F0F,0x3333,0x5555]+[a&b for a,b in itertools.combinations([0x00FF,0x0F0F,0x3333,0x5555],2)])
C12=[c for c in RM2 if c>>12==0]
cases.append((C12,12,4))
C13=[c&0x1FFF for c in RM2 if c>>13==0]
cases.append((C13,13,4))
for C,n,d in cases:
    for n1 in sorted(set([0,1,n//3,n//2])):
        vals,M=code_x(C,n,n1)
        P=split.build(n,d,n1)
        print(n,d,n1,M,'ok %.2e'%check(P,vals,M),flush=True)
