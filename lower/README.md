# Data for "Orbit unions and exchange walks: seven new lower bounds for binary codes"
Zhiyang Chen, Ligeng Zhu, Fan Long (2026).

Files
- lower_bounds.pdf, lower_bounds.tex: the paper.
- codes/A_<n>_<d>_<M>.txt: the codes, one word per line (0/1 string, leftmost character = coordinate 0).
  A_17_6_260, A_21_6_3072, A_24_10_144, A_25_10_256, A_27_6_81920, A_28_6_139264, A_28_12_192 (new bounds);
  A_26_6_40960, A_25_6_20480 (shortenings of A_27_6_81920).
- verify_codes.py (standard library, for small codes) and verify_fast.c (C, OpenMP, for large codes):
    python3 verify_codes.py codes/A_24_10_144.txt 10
    gcc -O3 -fopenmp -o verify_fast verify_fast.c && ./verify_fast codes/A_28_6_139264.txt -d 6
- programs/: source code of the search programs (orbit graphs, coset unions, tabu/CHILS drivers, exchange walk).
- construction/: the data behind each construction (groups, kernels, Z4 base code z12_1.h with the map psi and the chosen
  cosets, the start codes and the program of the exchange walk). The README files there give the file formats.

The proof of each bound is the code file itself: n, M and the minimum distance can be checked directly.

Use of AI: AI tools were used in the search and in the writing. The correctness of the bounds does not depend on them:
each code is checked directly, and the authors guarantee the stated bounds.
