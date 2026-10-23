# A(27,4) <= 2396695  (live table 2396736; team 2396696 = doubling corollary)   equivalently A(26,3) <= 2396695
Agent sdp, 2026-10-23 (round 5, task L2). Schrijver 3-point SDP + two Best rows + B1 coset-count rows, exact certificate.
PREMISE: A(25,4) <= 599174 (novaB, results/verified/UB_A_25_4_599174_B1/, referee-audited).

## Proof
1. A(26,4) <= 2 A(25,4) <= 1198348 and A(27,4) <= 2 A(26,4) <= 2396696 (doubling: one coordinate splits C into two
   shortened codes).
2. Let C be an even-weight (27, M, 4) code with M = 2396696 exactly (puncture + parity from any (27,M,4) code).
   excl_n27_d4_M2396696_best_b1.json: Schrijver SDP + Best rows for k = 1, 2,
      sum_a C(27-a,k) A_a <= C(27,k) A(27-k,4)  with A(26,4) <= 1198348, A(25,4) <= 599174,
   + B1 rows (all block subspaces with rho = M mod 2^r != 0, plus 48 general V, stored in the certificate)
   -> |C| <= 2396695.999963518 < 2396696. No such code.
So A(27,4) <= 2396695. Every larger code has a subcode of size exactly 2396696.
Margin: 0.000036482 (small, but the bound is an exact rational computation; the solver status is not used).

## Rows
- Best row (BBMOS 1978): count (c, c', K), |K| = k, c = c' on K. For fixed (c, K) the codewords equal to c on K, with K
  deleted, form a (27-k, ., 4) code; for fixed (c, c') there are C(27 - d(c,c'), k) such K.
- B1 row for a linear V of dim r with m_k nonzero words of weight k, rho = M mod 2^r:
  sum_k m_k/C(n,k) sum_i A_i K_k(i) >= rho(2^r - rho)/M.
- Ablation: without the Best rows the value is 2396696 (not excluded); with Best rows but trivial premises
  (A(m,4) <= 2^(m-1)) for m = 25 also 2396696; with the real A(25,4), A(26,4): excluded. Live CW values: no effect.
- The same method at A(26,4), M = 1198348 (premise A(25,4) <= 599174): feasible (no gain).

## Check
    python3 -I check_cert.py excl_n27_d4_M2396696_best_b1.json      -> ALL OK (4.5 min, stdlib)
check_cert.py = UB_A_28_8_32148/check_cert.py + Best-row branch (regenerates the row from (k, premise) and requires the
premise to be the one stated in problem.best_premises). Tamper tests: M changed, premise changed, Best rows removed: rejected.
Generator: agents/sdp/make_cert_best.py (ONLYPREM=1), search: agents/sdp/cplin.py (L2 sweep, runs/l2sweep.log).
