# A(26,6) <= 84243  (table: 84260)   equivalently A(25,5) <= 84243
Agent sdp, 2026-10-22 (round 4). Method: Schrijver 3-point SDP + B1 coset-count rows (novaB, agents/novaB/METHOD_B1.md),
exact rational dual certificates. No premise: no table value and no classical (Best/MEL/cw) row is used.

## Proof
C: binary code, length 26, min distance >= 6. Puncture + parity bit: assume C even-weight (same n, M, d).
1. base_n26_d6_mel.json: plain Schrijver SDP gives |C| <= 84260.5714 < 84261.
2. excl_n26_d6_M{84260..84244}_b1mel.json (17 files): for |C| = M exactly, Schrijver SDP + B1 rows (all block
   subspaces with rho != 0, plus up to 615 general-V rows from the separation loop; each used row is stored in the
   certificate). Exact bounds: 84260: 84242.6862, 84259: 84242.9871, 84258: 84243.4070, 84257: 84243.9200, 84256: 84244.5222, 84255: 84243.9804, 84254: 84243.5570, 84253: 84243.2520, 84252: 84243.0654, 84251: 84242.9973, 84250: 84243.0475, 84249: 84242.9220, 84248: 84243.1111, 84247: 84243.1244, 84246: 84243.4487, 84245: 84243.5891, 84244: 84243.8303.
So A(26,6) <= 84243. Limit: M = 84243 is feasible for this relaxation (value exactly M).
Margins (M - exact bound): base 0.429; 84260: 17.314, 84259: 16.013, 84258: 14.593, 84257: 13.080, 84256: 11.478, 84255: 11.020, 84254: 10.443, 84253: 9.748, 84252: 8.935, 84251: 8.003, 84250: 6.952, 84249: 6.078, 84248: 4.889, 84247: 3.876, 84246: 2.551, 84245: 1.411, 84244: 0.170. Smallest: 0.170 (M=84244).

## Exactness (no reliance on solver status)
Each certificate holds integer matrices V (Y = V V^T / 4^K, PSD exactly by construction) and integer multipliers z >= 0
(z / 2^K). check_cert.py computes the Lagrangian bound in exact fractions; each residual coefficient is bounded with the box
0 <= x <= 1 (sum of positive parts). The solver status is not used. Solver: SDPA-DD (CSDP gave a false 'feasible').
positive_controls.log: the same relaxation does NOT exclude sizes of real codes (e.g. (27,6,81920), (28,6,139264)).

## Rows
B1 row for a linear V of dim r with m_k nonzero words of weight k, rho = M mod 2^r:
  sum_k m_k/C(n,k) sum_i A_i K_k(i) >= rho(2^r - rho)/M   (integer coset counts of V^perp).
The checker regenerates each row from its block list (oa_blocks) or generator columns (oa_cols).

## Check
    python3 -I check_cert.py base_n26_d6_mel.json excl_n26_d6_M842*_b1mel.json
-> ALL OK (check_all.log). check_cert.py = same file as in UB_A_28_8_32148.
Generators: agents/sdp/jobs/descend.sh (cplin.py, USE_SDPA=1, LINOPTS=none), agents/sdp/jobs/certs.sh (make_cert_b1mel.py).
Not a corollary: no team bound at n=25 or n=27 implies it (A(26,6) <= 2 A(25,6) = 95076 with table A(25,6) <= 47538).
