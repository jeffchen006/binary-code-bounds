# alg records: Z4 coset unions with subgroups and extra coordinates

## Records (all verified with tools/verify_fast, results/verified/)
| bound | dir | e | j | psi | cosets |
|---|---|---|---|---|---|
| A(26,6) >= 36864 | A_26_6_36864 | 2 | 2 | 0 1 1 0 3 2 0 2 2 0 3 3 | 9 x 4096 |
| A(26,6) >= 40960 | A_26_6_40960 | 2 | 2 | 1 2 1 1 2 1 3 0 0 0 2 1 | 10 x 4096 |
| A(27,6) >= 69632 | A_27_6_69632 | 3 | 2 | 0 3 3 1 3 3 2 1 3 2 1 2 | 17 x 4096 |
| A(27,6) >= 81920 | A_27_6_81920 | 3 | 2 | 1 2 1 1 2 1 3 0 0 0 2 1 | 20 x 4096 |
| A(28,6) >= 135168 | A_28_6_135168 | 4 | 2 | 1 2 1 1 2 1 3 0 0 0 2 1 | 33 x 4096 |

All use the base code `runs/z12_1.h` (copied in some record dirs).

## Construction
- G = Z4^12 (Gray map 0->00, 1->01, 2->11, 3->10; bit 0 = first bit of the pair).
- C = ker(H) <= Z4^12, |C| = 2^14, min Lee distance 6 (Gray image is a (24, 16384, 6) code).
- Code length n = 24 + e: e extra binary coordinates ("tails") placed FIRST, then the 12 Gray pairs.
- psi: G -> Z2^j, psi(x) = sum over coordinates i with x_i odd of psi_i (psi_i in Z2^j, given as an integer 0..2^j-1).
  C' = C cap ker(psi) has index 2^j in C.
- The code is the union of the listed cosets of C' x {0} in Z2^e x Z4^12 (Gray images).

## File formats
- `*.h` (zsyn output): line 1 `a b x y d` (G = Z2^a x Z4^b, Q = Z2^x x Z4^y), then a+b integers: the columns h_i of
  H as packed Q elements: Z2 bits [0,x), Z4 low bits [x, x+y), Z4 high bits [x+y, x+2y).
- `rec*.inst` (zctabu instance): `a b d parity` / `X Y` / packed images phi(e_i) of all coordinates in
  Q' = Z2^X x Z4^Y (Z2 part = Q-Z2 bits, then psi bits, then tail bits; then Z4 low, Z4 high) / `0`.
- `rec*.reps`: one packed Q' element per line = the chosen cosets (0 = C' itself).
- Rebuild: `agents/alg/bin/zctabu -i rec.inst -expand rec.reps -o code.txt` then `tools/verify_fast code.txt`.
  Or from (psi, sel): `python3 agents/alg/src/mkrec.py code.h e j "psi" "sel" outprefix`.

## Tools (agents/alg/bin, sources in agents/alg/src)
- `zsyn -a A -b B -x X -y Y -d D -t sec -s seed -o out.h`: SA for an additive code with min Gray distance >= D
  and |C| = 2^(A+2B) / |Q|. "final cost=0" = success.
- `zsubscan -h code.h -e E -j J -d D -mode rand -n N -tgt T -tl sec -bb [-k4 1] [-j0 J0 -psi0 "..."] [-init0 f] [-fix]`:
  for N random psi, builds the graph of cosets of C' compatible with coset 0 and finds a max clique
  (BBMC colour B&B, time limit tl, incumbent if timeout). Prints `m=` (cosets incl. 0); prints psi and sel when m > T.
  -k4 1: extra Z4-valued psi component. -psi0: fix the low J0 bits of psi. -init0: warm start (coset ids in the
  packing with J0 psi bits). -fix: keep the warm-start cosets and search only among cosets compatible with all of them
  (this "lengthening" gave the n=27 and n=28 records from the n=26 one).
- `chain.py`: watches logs for (psi, sel) and extends each to more extra coordinates with -fix.
- `shortbest.py code.txt`: best 1- and 2-coordinate shortenings.
