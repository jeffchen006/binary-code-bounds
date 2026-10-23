import random, sys
from validate2 import code_x, check
import parity
random.seed(5)
cnt=0
for trial in range(400):
    n=random.randint(6,10); d=random.choice([4,6])
    words=[w for w in range(2**n) if bin(w).count('1')%2==0]; random.shuffle(words)
    C=[]; lim=random.randint(2,60)
    for w in words:
        if all(bin(w^c).count('1')>=d for c in C): C.append(w)
        if len(C)>=lim: break
    M=len(C)
    s=0
    for c in C: s^=c
    w=bin(s).count('1')
    case=('odd' if M%2 else ('2mod4' if M%4==2 else '0mod4'))
    if case=='0mod4' and w==0: continue
    vals,_=code_x(C,n,0)
    P=parity.build(n,d,0,M,case,w)
    check(P,vals,M); cnt+=1
print('all ok',cnt)
