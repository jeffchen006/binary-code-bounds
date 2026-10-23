# A(18,6) <= 672   (Brouwer table 2026-10-21: 673; equivalently A(17,5) <= 672)

Agent sdp, 2026-10-21. Method: Gijswijt-Mittelmann-Schrijver 4-point SDP (IEEE IT 58 (2012) 2697-2705) combined
with an exact-size parity row; exact rational dual certificate, checked by an independent checker.

## Proof
1. A(18,6) <= 673: Gijswijt, Mittelmann & Schrijver 2012 (their A4(18,6) = 673.005).
2. No (18, 673, 6) code exists. WLOG all weights even (puncture + parity bit), so all distances are even.
   Let S_u = sum_{c in C} (-1)^{u.c}. Since |C| = 673 is odd, every S_u is odd, so
   sum_{|u|=k} S_u^2 = |C| sum_i A_i K_k(i) >= C(18,k), i.e. sum_i A_i K_k(i) >= C(18,k)/673 for k = 0..18.
   The GMS 4-point SDP (M_empty block of codes of size <= 2, M_{0,u} blocks for every pair distance, Schrijver's
   3-point blocks, x >= 0) plus these 19 rows has maximum 2^18 p <= 672.88256 < 673, where p = x({v}) = |C|/2^18.
   Certificate: excl_n18_d6_M673_odd.json. Contradiction.
3. Hence A(18,6) <= 672.

## Check
    python3 check_gms_cert.py excl_n18_d6_M673_odd.json      # ~1.5 min, prints bound 672.882556574 < 673, ALL OK
    python3 check_gms_cert.py --selftest                      # block formulas vs. full 2^n matrices, n <= 7
The checker is independent of the solver code. It rebuilds every PSD block entry (Schrijver beta, Terwilliger
tensor blocks for M_{0,u}, GMS Krawtchouk formula (41) for M_empty, invariant block with the empty set), and the
valid linear rows. It accepts only rows from these families, uses Y = V V^T / 4^K (PSD by construction) and
z >= 0, and evaluates the Lagrangian bound exactly with the box 0 <= x(S) <= p = 673/2^18. That box holds for any code
of size exactly 673.

Numerics: SDPA-DD (double-double), 4 threads, ~2 h. Code: agents/sdp/gms.py, gms0.py, parity.py, gmspar.py.
