# B1: orthogonal-array (coset-count) integrality rows for SDP upper bounds  (novaB, 2026-10-21)

## Statement
C binary code, |C| = M (exactly), length n. S_u = sum_{c in C} (-1)^{u.c}. Let V <= F2^n be a subspace of dim r.
The map c -> (u_1.c, ..., u_r.c) (u_i a basis of V) sorts C into 2^r classes (cosets of V^perp) with counts n_g:
nonnegative INTEGERS with sum M. Fourier/Parseval on F2^r: sum_{v in V} S_v^2 = 2^r sum_g n_g^2.
Min of sum n_g^2 over integers is the balanced split, so with rho = M mod 2^r:

    sum_{v in V, v != 0} S_v^2  >=  rho (2^r - rho).                                     (OA_V)

(Equality iff the counts differ by at most 1. If rho != 0, C cannot be an orthogonal array "on V";
r = 1 gives the known M-odd row S_u^2 >= 1.)
Averaging (OA_V) over the S_n-orbit of V and using sum_{|u|=k} S_u^2 = M sum_i A_i K_k(i):

    sum_k  m_k(V)/C(n,k) * sum_i A_i K_k(i)  >=  rho (2^r - rho) / M,      A_0 = 1,       (row)

where m_k(V) = number of nonzero words of weight k in V. Linear in the distance distribution ->
add to Delsarte LP / Schrijver SDP / GMS SDP. Exact size M, but NO case split (unlike Kim-Toan 2 mod 4).
Used family here: "block" subspaces V = span of r vectors with disjoint supports of sizes a_1..a_r
(sum a_i <= n); m_k = #subsets of blocks with total size k. r = n (all blocks size 1) is Parseval (|C| <= M).

## Why it bites
The Schrijver optimum for A(17,4) (3276.8) has dual distribution beta_k = 0 for k = 1..6 (an OA of strength 6).
A real code of size 3276 = 4 mod 8 cannot be balanced on 3 coordinates: rows force S-mass into low shells,
which costs objective. With all block rows (1122 rows), the plain Schrijver SDP max for M = 3276 is 3275.9406 < 3276.

## Result: A(17,4) <= 3275  (table 3276, Schrijver 2005; GMS 4-point SDP also 3276.8)
Proof: if A(17,4) >= 3276 take a subcode of size exactly 3276, make it even-weight (puncture + parity).
Schrijver SDP (only Schrijver's own constraints) + block OA rows for M = 3276: exact rational dual certificate gives
|C| <= 3275.94065 < 3276. Contradiction.
Files: cert_A17_4/excl_n17_d4_M3276_oa.json; checker check_oa_cert.py (sdp's audited Schrijver checker
+ OA branch that regenerates every OA row from its block list); selftest_oa.py (24386 checks of the rows on real
codes, 0 violations, 2228 equality cases); tamper tests: changed M / changed row / wrong blocks -> rejected;
certificate without the OA rows -> bound 146784 (FAIL), so the new rows are essential.

## Code
oa.py (row builder on agents/sdp/split.py), t3.py (block partitions), make_cert_oa.py (CSDP solve + rational cert),
check_oa_cert.py, selftest_oa.py, t4.py / batch*.sh (scans).

## Related work (closest)
- AVZ 2001 Thm 4, MEL, Kim-Toan 2013 Thm 6: rows for M odd (= r=1 here) and M = 2 mod 4 with case split on w.
- Orthogonal arrays: Rao bound / divisibility 2^t | N for strength t; Bose-Bush (1952) integrality-of-counts
  argument for OAs of index lambda. Here the integrality is used quantitatively, for all r and all residues,
  as cutting planes in the code SDP. Not found in the code-bound literature (search continuing).

## Update 2026-10-21 13:40: general V, cap row, more results
- General V: any linear V of dim r, given by the multiplicities of the nonzero columns of an r x l generator matrix
  (l <= n). Row depends only on the weight enumerator of V. Separation: local search over column multisets, r = 2..10,
  objective sum_{x != 0} beta*_{w(x)} - rho(2^r - rho) (beta* = shell averages of S^2 at the SDP optimum). cp.py, cpscan.py.
- Cap row (premise A(n-1,d) <= Ap): |S_{e_i}| <= 2Ap - M.
- Rigorous records: A(17,4) <= 3272, A(18,4) <= 6543, A(19,4) <= 13085 (audited by novaA), A(21,4) <= 43681,
  A(22,6) <= 6939, A(24,6) <= 24102. Scans: A(25,4) < 599184, A(26,4) < 1198368 (certs pending).
- No effect: where the table value is below the Schrijver value (GMS-type entries: A(22..24,4), A(23,6), A(25,6)),
  and d = 10, 12 entries (no violated V at the table value), A(27,4) (2396736 stays feasible).
- Numerical caution: a CSDP "Partial Success" solve gives a weak certificate (A(24,6) with an enlarged cut pool);
  certificates are built with exactly the cut pool of a "Success" solve. Only rigorous bounds count.
- Shortened-part rows (B1 on {c: c_i = b} with sizes fixed by a premise) gave nothing extra at n=18, M=6543.

## Novelty (constr's check, agents/constr/NOVELTY_B1.md) and ablation (ABLATION.md), 2026-10-21
- Known: r = 1 (M odd) = Best-Brouwer-MacWilliams-Odlyzko-Sloane 1978 Thm 5; r = 2, M = 2 mod 4 = BBMOS Thm 7;
  the 2 mod 4 split = BBMOS Thm 8 (AVZ 2001, Kim-Toan 2013). The cap row is Best/MEL-type: not claimed.
- New: r >= 3 rows for all residues (incl. M = 0 mod 4), general V via the weight enumerator, separation loop, use in
  Schrijver's SDP with exact certificates. Claim: "a generalization of the BBMOS parity constraints".
- Ablation: 33 of 34 exclusions in the filed records need the new rows (only A(21,4) M = 43690 is obtained by the old rows).
