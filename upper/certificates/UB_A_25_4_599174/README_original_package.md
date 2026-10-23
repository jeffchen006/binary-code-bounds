# A(25,4) <= 599174   (Brouwer table 2026-10-21: 599184)

Agent novaB, 2026-10-21. Method B1 (agents/novaB/METHOD_B1.md; rows explained in results/verified/UB_A_17_19_4_B1/README.md):
Schrijver 3-point SDP (Schrijver's own constraints only) + coset-count integrality rows (general linear V, no premise).
## Proof (WLOG even-weight code; a code with > M words has a subcode of exactly M words)
1. base_n25_d4.json: plain Schrijver SDP |C| <= 599186.288  =>  A(25,4) <= 599186.
2. For M = 599186 down to 599175: excl_n25_d4_M<M>_oa.json gives a rigorous bound < M (no code of size exactly M).
   Rigorous bounds: 599185.27, 599179.43, 599175.12, 599176.40, 599173.63, 599174.41, 599173.86, 599176.64, 599174.07,
   599174.37, 599174.70, 599173.72 (M = 599186 .. 599175).    =>  A(25,4) <= 599174.  (Equivalently A(24,3) <= 599174.)
Note: some CSDP solves ended with "reduced accuracy"; this only weakens the rounded certificate, and every rigorous bound
above is still < M.
## Check
    python3 -I check_oa_cert.py *.json        (exact rational; output CHECK_novaB.txt)
