# A(28,6) >= 139264 (team 135168, table 131072) — novaC, GSCU engine, 2026-10-21
Code: results/verified/A_28_6_139264.txt (n=28, 139264 words, dmin=6; tools/verify_fast + check_code.py RECORD).
Structure: union of 34 cosets of the 4096-word subgroup C' = C cap ker(psi) of alg's Z4^12 code C (z12_1, |C| = 2^14, Gray (24,16384,6)),
with e = 4 extra binary coordinates; psi = (1 2 1 1 2 1 3 0 0 0 2 1) — exactly alg's instance rec28.inst (agents/alg/records/A_28_6_135168).
Q = G/C' has 2^16 cosets. alg found 33 cosets (lengthening, coset-0 neighbourhood). GSCU searched the FULL Cayley graph Cay(Q, S)
(all 65536 cosets, no forced coset 0) with clique tabu: 34 cosets in 120 s.
Rebuild: agents/novaC/gscu/gscu expand -i rec28.inst -sel rec28_f.sel -o code.txt ; tools/verify_fast code.txt
Reproduce search: OMP_NUM_THREADS=16 gscu solve -i rec28.inst -t 0 -ls 120 -noexact -vmax 70000 -o x.sel   (seed 1)
