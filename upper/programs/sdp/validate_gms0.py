import random, sys
from validate_gms import gms_vals
from validate2 import check
import gms0
random.seed(11)
for trial in range(10):
    n=random.randint(6,9); d=random.choice([2,4,4,6])
    words=[w for w in range(2**n) if bin(w).count('1')%2==0]; random.shuffle(words)
    C=[]; lim=random.randint(3,30)
    for w in words:
        if all(bin(w^c).count('1')>=d for c in C): C.append(w)
        if len(C)>=lim: break
    vals,M=gms_vals(C,n)
    P=gms0.build_full(n,d,M,light=True)
    print(n,d,M,'min rel eig %.2e'%check(P,vals,M),flush=True)
