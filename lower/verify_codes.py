#!/usr/bin/env python3
"""Check a binary code file: one word per line (0/1 string). Usage: python3 verify_codes.py FILE D
Prints n, M, number of distinct words and the minimum distance; exits 1 if some pair is closer than D.
Standard library only, O(M^2): use it for small codes. For the large codes (n=27, 28, d=6) use verify_fast.c:
  gcc -O3 -fopenmp -o verify_fast verify_fast.c && ./verify_fast FILE -d D"""
import sys
words = [l.strip() for l in open(sys.argv[1]) if l.strip()]
d = int(sys.argv[2]); n = len(words[0])
assert all(len(w) == n and set(w) <= set('01') for w in words)
V = sorted(set(int(w, 2) for w in words))
M = len(V); mind = n + 1
# bucket by the first k bits is not needed for correctness; plain O(M^2) with early exit per pair
for i in range(M):
    a = V[i]
    for j in range(i + 1, M):
        x = (a ^ V[j]).bit_count() if hasattr(int, 'bit_count') else bin(a ^ V[j]).count('1')
        if x < mind:
            mind = x
print(f"n={n} M={len(words)} distinct={M} min_distance={mind}")
sys.exit(0 if (M == len(words) and mind >= d) else 1)
