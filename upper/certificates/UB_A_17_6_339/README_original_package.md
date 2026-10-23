# A(17,6) <= 339  (table: 340)   equivalently A(16,5) <= 339
Agent sdp, 2026-10-21 (round 4). Method: Schrijver 3-point SDP + one MEL/Best-type row + B1 coset-count rows
(novaB, agents/novaB/METHOD_B1.md), exact rational dual certificates. No premise from the table is used.

## Proof
Let C be a binary code of length 17, min distance >= 6, |C| = M. Puncture one coordinate and add a parity bit:
an even-weight code, same length and size, min distance >= 6. So assume C even-weight.
1. base_n17_d6_mel.json: Schrijver SDP (Schrijver's own constraints) + the MEL row gives |C| <= 341.3333 < 342.
2. excl_n17_d6_M341_b1mel.json: for |C| = 341 exactly, + B1 block rows (r = 1 rows = the odd-M parity rows):
   max <= 340.4194 < 341. No such code.
3. excl_n17_d6_M340_b1mel.json: for |C| = 340 exactly, + B1 block rows + 6 general-V B1 rows (r = 3, 4, found by
   separation; list in pool_M340_generalV.json): max <= 339.5810 < 340. No such code.
So A(17,6) <= 339. Every code of size >= 340 has a subcode of size exactly 340 or 341; steps 2-3 exclude both.

## Rows
- MEL row (Mounits-Etzion-Litsyn 2002 type; d = 2h, q = floor(n/h)): A_{n-h} + q sum_{i<h} A_{n-i} <= q.
  Here A_14 + 5 A_16 <= 5. Proof per codeword in check_cert.py (mel_row). Constant-free.
- B1 row for a linear V of dim r with m_k nonzero words of weight k, rho = M mod 2^r:
  sum_k m_k/C(n,k) sum_i A_i K_k(i) >= rho(2^r - rho)/M (integer coset counts of V^perp). The checker regenerates
  each row from its block list (oa_blocks) or generator columns (oa_cols).
- Ablation (CSDP): M=340 with Schrijver + B1 only: 340.00005 (not excluded); + MEL row only: 339.581.
  Without the 6 general-V rows: the certificate bound is 48214 (FAIL), so they are essential.

## Check
    python3 -I check_cert.py base_n17_d6_mel.json excl_n17_d6_M341_b1mel.json excl_n17_d6_M340_b1mel.json
-> ALL OK (8 s, stdlib only). check_cert.py = novaB's audited check_oa_cert.py (itself sdp's check_sdp_cert.py + OA rows)
plus the MEL branch. Tamper tests (M changed, MEL flag removed, a V changed, general-V rows removed, MEL constant
changed): all rejected.
Generators: agents/sdp/cplin.py (cutting-plane search), agents/sdp/make_cert_b1mel.py (certificates).
