# Data for "A cutting-plane method for exact-size code bounds, with fifteen new upper bounds for binary codes"
Zhiyang Chen, Ligeng Zhu, Fan Long (2026).

Files
- upper_bounds.pdf, upper_bounds.tex: the paper.
- certificates/UB_A_<n>_<d>_<U>/: for each new bound A(n,d) <= U, the decisive certificate (no code of size exactly U+1),
  the exact checker (Python standard library only) and the README of the original package.
  Large JSON files are compressed: run `unxz *.json.xz` first.
- programs/: the generator programs (SDP builders, cutting-plane loop, certificate export; agents sdp, novaB, exact).

Check one bound (example):
    cd certificates/UB_A_17_6_339 && python3 check_cert.py excl_n17_d6_M340_b1mel.json        # -> bound 339.58... < 340, ALL OK
Checkers: check_oa_cert.py (Schrijver SDP + coset-count rows, optional cap row), check_cert.py (same + MEL / Best rows),
check_gms_cert.py (GMS 4-point SDP + odd-size rows). Each checker rebuilds every constraint from (n, d, M) and the stored
subspaces, uses PSD multipliers Y = V V^T / 4^K with integer V, and computes the bound in exact rational arithmetic.
The premise-based bounds use: A(18,4): A(17,4) <= 3272; A(19,4): A(18,4) <= 6544; A(27,4): A(26,4) <= 1198348 and A(25,4) <= 599174.

Use of AI: AI tools were used in the computations and in the writing. The authors do not guarantee the correctness
of the results.
