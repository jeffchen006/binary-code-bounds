# A(17,6) >= 260 (lns tie farm, 2026-10-23 07:53) — supersedes A(17,6) >= 259

File: results/verified/A_17_6_260.txt (md5 605eae3ccf06eb9338caed1d8179599f). n=17, M=260, dmin 6.
Distances 6:11002 8:7213 10:12888 12:2439 14:3 16:125. Nonlinear, no structure imposed.

## Chain (all steps = the exchange walk of ../A_17_6_259/README.md: equal-size moves + exact re-fill, k 4..60, K=0 full ball)
1. evo: (17,257,6) code harvested from its Z5 tie farm (tf/a17/z5_s7510.code) -> near-miss store tie deaa1529a77b
   (file seed_257_deaa1529a77b_evo.txt).
2. lns tie farm, `nearmiss.py walk next` on deaa1529a77b, lnsx seed 1898201803, puncture p=4: 60-s probe 257 -> 258, then a 20-min turn
   257 -> 259 (1.33 M moves, 388 k plateau). The farm stored the 259 code as new tie c06e91bdb1d1 (seed_259_c06e91bdb1d1.txt).
3. lns tie farm on c06e91bdb1d1, lnsx seed 414673923, p=12: 20-min turn (1.76 M moves, 360 k plateau) -> 260.
   Raw output: found_260_original.txt (= the verified file). walks_chain.jsonl / index_chain.jsonl = the store rows.
## Reproduce
lnsx -i seed_259_c06e91bdb1d1.txt -d 6 -K 0 -kmin 4 -kmax 60 -t 3000 -nodes 100000 -maxF 4000 -rec 259 -s 414673923
(deterministic xorshift; the farm ran the same seed twice: a 60-s probe, then the 20-min turn).
Re-run DONE: 260 reached again; repro_260.txt is IDENTICAL to found_260_original.txt (cmp); repro.log.
lnsx.cpp / nearmiss.py = the program versions used.
## Lesson
Growing a POOL of different tie codes (every walk's grown code is stored as a new seed) and walking each one is what gave 260:
the record seed and the first 259 codes walked > 10 M moves without 260; this 259 came from a 257 of a different evo family.
