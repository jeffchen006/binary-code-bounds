# A(24,10) >= 144 (evo, round 4, 2026-10-22 ~02:00); table lower 136, upper 268
Code: results/verified/A_24_10_144.txt (n=24, 144 words, dmin 10; distances 10,12,14,16,18,20).
= parity extension of oddform_23_144_9.txt (n=23, 144 words, dmin 9).
Structure: union of 18 cosets of a 3-dim linear kernel L (kernel.L, min weight >= 9) in F2^23; L is a subcode of the binary Golay
code [23,12,7] (parent_golay_23_12_7.G = novaE lib LIN_n23_d7_M4096_guava-bkl.base). NONLINEAR (18 cosets is not a power of 2).
Found by evo's generator program k6_parent2 (agents/evo/pop/progs/k6_parent2.py -> k5_parent.py), seed 16:
  step 1: agents/evo/cls3p -G parent -L kernel -d 9 -t 120 -s 16   (cosets of L INSIDE the Golay code: 4 cosets = 32 words)
  step 2: agents/evo/cos/cls3 -m 23 -d 9 -L kernel -i <step-1 code> -t 600 -s 16   (full space, seeded: 18 cosets = 144)
  (cls3 = d10's coset tabu, copy in agents/evo/cos; cls3p = evo's parent-restricted copy, src/cls3p.cpp)
Reproduce: cd agents/evo; python3 src/gloop.py run k6_parent2 1 --cpus <cpu> --seed0 16   (heuristic; the same seed should repeat)
Checks: tools/verify (n=24 M=144 dmin=10), tools/check_code.py RECORD, independent Python pairwise check (dmin 10, 144 distinct).
Corollaries (not separate results): A(23,9) >= 144 (same code); shortening gives A(23,10) >= 72 (< 84, nothing new).
