# A(27,6) <= 157281  (table: 157285)   equivalently A(26,5) <= 157281
Agent sdp, 2026-10-22 (round 4). Method: Schrijver 3-point SDP + B1 coset-count rows (novaB, agents/novaB/METHOD_B1.md),
exact rational dual certificates. No premise: no table value and no classical (Best/MEL/cw) row is used.

## Proof
C: binary code, length 27, min distance >= 6. Puncture + parity bit: assume C even-weight (same n, M, d).
1. base_n27_d6_mel.json: plain Schrijver SDP gives |C| <= 157286.4000 < 157287.
2. excl_n27_d6_M{157286..157282}_b1mel.json: for |C| = M exactly, Schrijver SDP + B1 rows (all block subspaces with
   rho != 0, plus general-V rows r = 2..10 from the separation loop; each used row is stored in the certificate):
   exact bounds 157280.7252 (M=157286), 157280.8356 (157285), 157281.1696 (157284), 157281.4693 (157283), 157281.9436 (157282).
So A(27,6) <= 157281. Limit: M = 157281 is feasible for this relaxation (value exactly M).
Margins (target - exact bound): base 0.600; M=157286: 5.27; 157285: 4.16; 157284: 2.83; 157283: 1.53; 157282: 0.056.

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
    python3 -I check_cert.py base_n27_d6_mel.json excl_n27_d6_M15728[2-6]_b1mel.json
-> ALL OK (check_all.log). check_cert.py = same file as in UB_A_28_8_32148.
Generators: agents/sdp/jobs/descend.sh (cplin.py, USE_SDPA=1, LINOPTS=none), agents/sdp/jobs/certs.sh (make_cert_b1mel.py).
Not a corollary: A(28,6) <= 2 A(27,6) etc. go the other way; A(27,6) <= 2 A(26,6) gives only ~168486.
