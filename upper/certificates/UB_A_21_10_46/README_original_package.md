# A(21,10) <= 46   (Brouwer table 2026-10-21: 47; equivalently A(20,9) <= 46)

Agent sdp, 2026-10-21. Method: Gijswijt-Mittelmann-Schrijver 4-point SDP (IEEE IT 58 (2012) 2697-2705) combined
with an exact-size parity row; exact rational dual certificate, checked by an independent checker.

## Proof
1. A(21,10) <= 47: Gijswijt, Mittelmann & Schrijver 2012 (their A4(21,10) = 47.007). Our own certificate:
   no (21,48,10) code, GMS SDP max 2^21 p <= 47.0069 < 48 (base_n21_d10_T48.json).
2. No (21, 47, 10) code exists. WLOG all weights even (puncture + parity bit), so all distances are even.
   Let S_u = sum_{c in C} (-1)^{u.c}. Since |C| = 47 is odd, every S_u is odd, so
   sum_{|u|=k} S_u^2 = |C| sum_i A_i K_k(i) >= C(21,k), i.e. sum_i A_i K_k(i) >= C(21,k)/47 for k = 0..21.
   The GMS 4-point SDP (M_empty block of codes of size <= 2, M_{0,u} blocks for every pair distance, Schrijver's
   3-point blocks, x >= 0) plus these 22 rows has maximum 2^21 p <= 46.91509 < 47, where p = x({v}) = |C|/2^21.
   Certificate: excl_n21_d10_M47_odd.json. Contradiction.
3. Hence A(21,10) <= 46.

## Check
    python3 check_gms_cert.py *.json      # ~4 min: base 47.006895450 < 48, excl 46.915093170 < 47, ALL OK
    python3 check_gms_cert.py --selftest                      # block formulas vs. full 2^n matrices, n <= 7
The checker is independent of the solver code. It rebuilds every PSD block entry (Schrijver beta, Terwilliger
tensor blocks for M_{0,u}, GMS Krawtchouk formula (41) for M_empty, invariant block with the empty set), and the
valid linear rows. It accepts only rows from these families, uses Y = V V^T / 4^K (PSD by construction) and
z >= 0, and evaluates the Lagrangian bound exactly with the box 0 <= x(S) <= p = 47/2^21. That box holds for any code
of size exactly 47.

Numerics: SDPA-QD (quad-double; SDPA-DD fails on this problem, as GMS also reported), 166 iterations, ~3.5 h.
Code: agents/sdp/gms.py, gms0.py, parity.py, gmspar.py, export_gms_cert.py.
