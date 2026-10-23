# A(24,6) <= 24102   (Brouwer table 24106; team bound 24105, results/verified/UB_A_24_6_24105)

Agent novaB, 2026-10-21. Method B1 (agents/novaB/METHOD_B1.md; see results/verified/UB_A_17_19_4_B1/README.md for the rows).
## Proof (WLOG even-weight code; a code with > M words has a subcode of exactly M words)
1. base_n24_d6.json: plain Schrijver SDP: |C| <= 24107.8858  =>  A(24,6) <= 24107.
2. M = 24107 (odd) and M = 24106 (2 mod 4, all 13 even w): agent sdp's certificates (copied from UB_A_24_6_24105, Kim-Toan/AVZ
   parity rows)  =>  A(24,6) <= 24105.
3. M = 24105, 24104, 24103: Schrijver SDP + B1 coset-count rows (general V, no premise): rigorous bounds
   24102.635, 24103.109, 24102.898 < M  =>  A(24,6) <= 24102.   (Equivalently A(23,5) <= 24102.)
## Check
    python3 -I check_oa_cert.py *.json     (handles base, odd, 2mod4 and B1 rows; output CHECK_novaB.txt)
