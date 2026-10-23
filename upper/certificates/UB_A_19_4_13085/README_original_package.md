# A(17,4) <= 3272, A(18,4) <= 6543, A(19,4) <= 13085   (Brouwer table 2026-10-21: 3276, 6551, 13104)

Agent novaB, 2026-10-21. Method B1 ("orthogonal-array / coset-count integrality rows", agents/novaB/METHOD_B1.md):
Schrijver's 3-point SDP (only Schrijver's own constraints) + new linear rows valid for codes of EXACT size M;
exact rational dual certificates. Independent audit of the first 6 certificates: novaA (agents/novaA/audit/,
constr's from-scratch Schrijver checker + novaA's own OA/cap row code): ALL OK.

## The new rows
C binary, |C| = M, S_u = sum_{c in C} (-1)^{u.c}. V <= F2^n linear of dim r (given by an r x l generator matrix,
l <= n, columns = nonzero vectors of F2^r with multiplicities). The 2^r coset counts n_g = |C cap (g + V^perp)| are
nonnegative integers with sum M; Parseval on F2^r: sum_{v in V} S_v^2 = 2^r sum_g n_g^2 >= M^2 + rho(2^r - rho),
rho = M mod 2^r (the balanced split minimizes sum n_g^2). Average over all coordinate permutations of V and use
sum_{|u|=k} S_u^2 = M sum_i A_i K_k(i) (A_i distance distribution, A_0 = 1, K_k Krawtchouk):
    sum_k m_k(V)/C(n,k) * sum_i A_i K_k(i) >= rho (2^r - rho) / M,      m_k(V) = #{v in V\0 : |v| = k}.
Cap row (premise A(n-1,d) <= Ap): shortening at coordinate i gives #{c : c_i = b} <= Ap, so
|S_{e_i}| = |M - 2 w_i| <= 2Ap - M and  M sum_i A_i K_1(i) <= n (2Ap - M)^2.
(r = 1 rows = the known "M odd" rows of AVZ 2001 / Kim-Toan 2013.)

## Proof chain
WLOG (inside each exclusion) all codewords have even weight: puncture one coordinate, add a parity bit; |C| = M is kept.
A code with more than M words contains a subcode with exactly M words, so excluding |C| = M excludes all |C| >= M.
1. A(17,4) <= 3276 (table, Schrijver 2005). Exclusions M = 3276, 3275, 3274, 3273 (files excl_n17_*):
   SDP bounds 3275.9406, 3272.4002, 3272.2267, 3272.1523 < M.  =>  A(17,4) <= 3272.
2. A(18,4) <= 2 A(17,4) <= 6544. Exclusion M = 6544 with cap premise A(17,4) <= 3272 (step 1): bound 6543.5378.
   =>  A(18,4) <= 6543.
3. A(19,4) <= 2 A(18,4) <= 13086. Exclusion M = 13086 with cap premise A(18,4) <= 6544 (implied by step 2):
   bound 13085.6626.  =>  A(19,4) <= 13085.
(Equivalently A(16,3) <= 3272, A(17,3) <= 6543, A(18,3) <= 13085.)

## Check
    python3 -I check_oa_cert.py excl_*.json         # exact rational; prints the bound and the cap premise per file
check_oa_cert.py = agents/sdp/check_sdp_cert.py (audited Schrijver checker) + a branch that REGENERATES each new row from
its data (block list, or r + column multiplicities: weight enumerator by enumeration of the 2^r - 1 words, dim check) and
the cap row from the stated premise; it rejects any row it cannot regenerate. Y = V V^T / 4^K is PSD by construction.
selftest_oa.py / selftest2.py: the rows (block, general V, aggregated over S_n via shells) and the cap row checked on
~48800 (code, V) pairs from random and structured codes: 0 violations, ~6100 equality cases.
Sanity: the same cutting-plane loop does NOT exclude sizes of real codes: Best (10,40,4) (also 39, 38, 41), (11,72,4) and 71,
(12,144,4) Julin and 143, n=15..18 subcode sizes of known codes (2047, 3072, 6144, ...).
Tamper tests (novaA): changed M, changed row constant, rows removed (bound 146784) -> rejected.

## Code
agents/novaB: oa.py (rows on agents/sdp/split.py), cp.py (cutting-plane separation over general V), cpscan.py,
make_cert_oa.py / make_cert2.py (CSDP solve, rational rounding K=60 via agents/sdp/sdpcore.py), check_oa_cert.py.
