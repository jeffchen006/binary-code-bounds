import random, sys, itertools
from validate2 import code_x, check
import covering
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 3)
def span(gens):
    S={0}
    for g in gens: S|={s^g for s in S}
    return sorted(S)
cases=[]
# Hamming [7,4,3] perfect, shortened Hamming codes, random greedy codes (lexicodes are good codes)
H7=span([0b1110000,0b1001100,0b0101010,0b1101001])
cases.append((H7,7,3))
cases.append(([c>>1 for c in H7 if c&1==0],6,3))
# lexicodes
for N,d in [(8,3),(9,3),(10,3),(11,3),(9,5),(10,5),(11,5),(12,5)]:
    C=[]
    for w in range(2**N):
        if all(bin(w^c).count('1')>=d for c in C): C.append(w)
    cases.append((C,N,d))
for trial in range(6):
    N=random.randint(7,11); d=random.choice([3,5])
    words=list(range(2**N)); random.shuffle(words); C=[]
    for w in words:
        if all(bin(w^c).count('1')>=d for c in C): C.append(w)
    cases.append((C,N,d))
for C,N,d in cases:
    vals,M=code_x(C,N,0)
    P=covering.build(N,d,T=M)
    print(N,d,M,'ok %.2e'%check(P,vals,M),flush=True)
