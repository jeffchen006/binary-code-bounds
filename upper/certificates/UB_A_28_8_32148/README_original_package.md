# A(28,8) <= 32148  (table: 32151, Schrijver 2005)   equivalently A(27,7) <= 32148
Agent sdp, 2026-10-21 (round 4). Method: Schrijver 3-point SDP + B1 coset-count rows (novaB, agents/novaB/METHOD_B1.md),
exact rational dual certificates. No premise: no table value and no classical (Best/MEL/cw) row is used.

## Proof
C: binary code, length 28, min distance >= 8. Puncture + parity bit: assume C even-weight (same n, M, d).
1. base_n28_d8_mel.json: plain Schrijver SDP gives |C| <= 32151.7958 < 32152.
2. excl_n28_d8_M{32151,32150,32149}_b1mel.json: for |C| = M exactly, Schrijver SDP + B1 rows (all block subspaces with
   rho != 0, plus 99 general-V rows r = 2..8 from the separation loop, list in pool_generalV.json):
   M = 32151 -> <= 32148.7689;  M = 32150 -> <= 32148.9031;  M = 32149 -> <= 32148.8778.  No such codes.
So A(28,8) <= 32148 (a code of size 32149..32151 has a subcode of size exactly 32149, 32150 or 32151... each excluded
directly; larger sizes are excluded by step 1).
Limit: M = 32148 is feasible for this relaxation (value exactly M, no violated V found).

## Rows
B1 row for a linear V of dim r with m_k nonzero words of weight k, rho = M mod 2^r:
  sum_k m_k/C(n,k) sum_i A_i K_k(i) >= rho(2^r - rho)/M   (integer coset counts of V^perp).
The checker regenerates each row from its block list (oa_blocks) or generator columns (oa_cols); mel_row unused here.

## Numerics note
The cutting-plane search first used CSDP: it stopped with "reduced accuracy" at 32151.03 (a false 'feasible').
SDPA-DD (double-double, opt/sdpa-gmp/sdpa_dd2) gives 32148.7689 and the exact certificate confirms it.
Ablation: adding Best/cw/MEL/KT rows does not change the value (32148.76888 in every combination).

## Exactness (no reliance on solver status)
Each certificate holds integer matrices V (Y = V V^T / 4^K, so every Y is PSD exactly, by construction) and integer
multipliers z >= 0 (z / 2^K). check_cert.py computes the Lagrangian bound in exact fractions; every residual
coefficient is bounded with the box 0 <= x <= 1 (sum of positive parts). The solver status (pFEAS etc.) is not used.
Margins (target - exact bound): base 0.204; M=32151: 2.231; M=32150: 1.097; M=32149: 0.122.
positive_controls.log: the same relaxation does NOT exclude sizes of real codes (e.g. (28,8,16384), (28,6,139264)).

## Check
    python3 -I check_cert.py base_n28_d8_mel.json excl_n28_d8_M32151_b1mel.json excl_n28_d8_M32150_b1mel.json excl_n28_d8_M32149_b1mel.json
check_cert.py = the checker of UB_A_17_6_339 (novaB's audited check_oa_cert.py + MEL branch), with one change: the
weight enumerator of a block subspace is computed as prod_b (1 + z^{a_b}) instead of a loop over 2^r subsets (r up
to 28); same rows. Run time 14 min (ALL OK, check_all.log).
Generators: agents/sdp/cplin.py (USE_SDPA=1, LINOPTS=none), agents/sdp/make_cert_b1mel.py (--nomel --sdpa).
