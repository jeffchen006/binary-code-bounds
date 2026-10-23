# A(28,6) <= 291265  (table: 291269)   equivalently A(27,5) <= 291265
Agent sdp, 2026-10-22 (round 4). Method: Schrijver 3-point SDP + B1 coset-count rows (novaB, agents/novaB/METHOD_B1.md),
exact rational dual certificates. No premise: no table value and no classical (Best/MEL/cw) row is used.

## Proof
C: binary code, length 28, min distance >= 6. Puncture + parity bit: assume C even-weight (same n, M, d).
1. base_n28_d6_mel.json: plain Schrijver SDP gives |C| <= 291271.1111 < 291272.
2. excl_n28_d6_M{291271..291266}_b1mel.json: for |C| = M exactly, Schrijver SDP + B1 rows (all block subspaces with
   rho != 0, plus general-V rows from the separation loop; each used row is stored in the certificate):
   exact bounds 291263.5047 (M=291271), 291263.8646 (291270), 291264.0464 (291269), 291264.4348 (291268),
   291264.8587 (291267), 291265.4825 (291266).
So A(28,6) <= 291265. Limit: M = 291265 is feasible for this relaxation (value exactly M).
Margins (target - exact bound): base 0.889; M=291271: 7.50; 291270: 6.14; 291269: 4.95; 291268: 3.57; 291267: 2.14; 291266: 0.517.

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
    python3 -I check_cert.py base_n28_d6_mel.json excl_n28_d6_M2912*_b1mel.json
-> ALL OK (check_all.log). check_cert.py = same file as in UB_A_28_8_32148.
Generators: agents/sdp/jobs/descend.sh (cplin.py, USE_SDPA=1, LINOPTS=none), agents/sdp/jobs/certs.sh (make_cert_b1mel.py).
Not a corollary: A(28,6) <= 2 A(27,6) gives only 314562 (with A(27,6) <= 157281).
