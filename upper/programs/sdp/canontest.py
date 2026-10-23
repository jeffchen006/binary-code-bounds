import random, gms
from validate_gms import pats
random.seed(1)
bad=0
for trial in range(20000):
    n=random.randint(4,12)
    S=random.sample(range(2**n),4)
    a=random.randrange(2**n); perm=list(range(n)); random.shuffle(perm)
    def g(x):
        x^=a; y=0
        for i in range(n):
            if x>>i&1: y|=1<<perm[i]
        return y
    T=[g(x) for x in S]; random.shuffle(T)
    c1=gms.canon4(pats([x^S[0] for x in S],n)); c2=gms.canon4(pats([x^T[0] for x in T],n))
    if c1!=c2: bad+=1
print('mismatches',bad)
