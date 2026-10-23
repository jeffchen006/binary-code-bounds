# A(22,6) <= 6939   (Brouwer table 2026-10-21: 6941)

Agent novaB, 2026-10-21. Method B1 (agents/novaB/METHOD_B1.md; same as results/verified/UB_A_17_19_4_B1/README.md):
Schrijver 3-point SDP (only Schrijver's own constraints) + coset-count integrality rows for codes of exact size M
(for a linear V of dim r: sum_{v in V\0} S_v^2 >= rho(2^r - rho), rho = M mod 2^r; averaged over S_n -> linear row
sum_k m_k(V)/C(n,k) sum_i A_i K_k(i) >= rho(2^r-rho)/M). General V found by cutting-plane separation (agents/novaB/cp.py).
No premise (no cap row) is used here. Self-contained: does not use any earlier upper bound.

## Proof
1. WLOG all codewords have even weight (puncture one coordinate, add parity; |C| kept), so all distances are even >= 6.
2. base_n22_d6.json: plain Schrijver SDP gives |C| <= 6943.695  =>  A(22,6) <= 6943.
3. For each M = 6943 down to 6940: file excl_n22_d6_M<M>_oa.json proves that no code of size exactly M exists
   (SDP + B1 rows: bound < M). A code with more than M words has a subcode with exactly M words.
   =>  A(22,6) <= 6939.

## Check
    python3 -I check_oa_cert.py base_*.json excl_*.json        # exact rational; output in CHECK_novaB.txt
(check_oa_cert.py = agents/sdp/check_sdp_cert.py + a branch that regenerates every B1 row from its data and rejects any
other new row.) Row validity tests: agents/novaB/selftest_oa.py, selftest2.py (0 violations on ~48800 (code, V) pairs).
