# A(21,6) >= 3072 (old lower bound 2560)

Code: results/verified/A_21_6_3072.txt (n=21, M=3072, dmin=6; tools/verify).
Odd form: os20_18_all_M3072.txt (n=20, M=3072, dmin=5); add an overall parity bit -> n=21, d=6.

## Construction
- Coordinates 0..15 = 8 Gray pairs (Z4 coordinates, Gray map 0->00, 1->01, 2->11, 3->10); coordinates 16..19 = 4 binary "label" coordinates.
- Octacode O (Z4, length 8, 256 words, Gray image = Nordstrom-Robinson). Each Z4-translation by an octacode word acts on F2^16 as a signed permutation (isometry).
- H = an order-64 subgroup of O (H = Z4^3), acting on coords 0..15 only (labels fixed). Generators (perm images + translation bits) in os20_18_all.gens
  (format: x -> perm(x) xor t, bit i goes to position perm[i]).
- The code is the union of 48 H-orbits (each of size 64) in F2^20. Orbit representatives: os20_18_all.reps (integer, bit i = coordinate i, plus orbit size);
  the chosen orbits are the 1-lines in os20_18_all.sol.
- Every label class (fixed coords 16..19) has 192 words with dmin 5.

## How it was found
orbgraph (src/orbgraph.cpp) builds the conflict graph of H-orbits (V=16384), ktabu (src/ktabu.cpp, k-fixed tabu) finds 48 independent orbits in < 2 min.
Reproduce:
  ./orbgraph -n 20 -d 5 -all -gens records/A_21_6_3072/os20_18_all.gens -out /tmp/x
  ./ktabu /tmp/x.graph -m 20 -wmin 64 -time 120 -seed 18 -maxni 50000 -out /tmp/x.sol   # re-checked: 48 orbits = 3072 at t=63 s
  ./orbgraph -n 20 -d 5 -all -gens records/A_21_6_3072/os20_18_all.gens -expand /tmp/x.reps /tmp/x.sol code20.txt

## Context (exact results, Cliquer / nbclique)
- Full O (order 256) on F2^20 with labels fixed: max invariant code = 2560 exactly (this is the old bound).
- Full O with any diagonal label homomorphism phi: O -> Z2^4 (all 3876 classes): max <= 2816.
- Subgroups of O up to conjugacy in Iso: 24 classes; only the class of H gives 3072; order-32 classes give <= 2816 (heuristic).

## Reusing the tools (agents/d6, all C++ built with g++ -O3)
- Group file format (one generator per line): `p0 p1 ... p(n-1) t b0b1...b(n-1)`; action x -> perm(x) xor t,
  bit i of x moves to position p_i, then translation t. Z4 coordinate j = bit pair (2j, 2j+1), Gray 0->00,1->01,2->11,3->10;
  Z4 translation by u on a pair: u=1 -> swap pair, t on bit 2j+1; u=3 -> swap, t on bit 2j; u=2 -> t on both bits.
- Octacode translations on 16 coords: groups/octT_16.txt (4 generators); padded with fixed coords: groups/octT_N.txt (N=17..23).
- Subgroup classes of the octacode translation group T up to conjugacy in Iso(NR) (24 classes, orders 16..128):
  groups/octcls_N.txt (N = 16..23, coords >= 16 fixed); same with independent generators: groups/octind_N.txt.
  GAP scripts: runs/oct_classes.g, runs/oct_ind.g (uses groups/nr17_gap.g = Iso(NR) subgroup of order 688128).
- orbgraph -n N -d D [-all] -gens g.txt -out pre : orbits, METIS conflict graph pre.graph (vertex weight = orbit size), pre.reps.
  (-all = all words; default = even words only, for even D). -expand pre.reps sol.txt code.txt ; -mark code.txt init.txt.
- ktabu pre.graph -m M -wmin W -time T -seed S [-init init.txt] : k-fixed tabu MIS on orbits of size >= W (+ greedy extension).
- nbclique pre.graph 0 nodelimit lb : exact max clique through orbit 0 (valid when the orbit graph is vertex-transitive).
- Drivers: src/scan_groups.py (ktabu + Cliquer exact for V<=4000), src/phi_scan.py / phi_scan2.py (diagonal label homs).
