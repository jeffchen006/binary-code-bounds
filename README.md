# New bounds on A(n,d) for binary codes, n <= 28

A(n,d) is the largest size of a binary code of length n and minimum distance d.
Reference table: A. E. Brouwer, https://aeb.win.tue.nl/codes/binary-1.html

## Papers (Zenodo, open access)
- Z. Chen, L. Zhu, F. Long, *A cutting-plane method for exact-size code bounds, with fifteen new upper bounds for binary codes*, 2026.
  DOI [10.5281/zenodo.23599648](https://doi.org/10.5281/zenodo.23599648) - `paper/upper_bounds.pdf`
- Z. Chen, L. Zhu, F. Long, *Orbit unions and exchange walks: seven new lower bounds for binary codes*, 2026.
  DOI [10.5281/zenodo.23599650](https://doi.org/10.5281/zenodo.23599650) - `paper/lower_bounds.pdf`

Please cite the papers (see `CITATION.cff`).

## New lower bounds
| Bound | Previous | New |
|---|---|---|
| A(17,6) >= | 258 | 260 |
| A(21,6) >= | 2560 | 3072 |
| A(24,10) >= | 136 | 144 |
| A(25,10) >= | 196 | 256 |
| A(27,6) >= | 65536 | 81920 |
| A(28,6) >= | 131072 | 139264 |
| A(28,12) >= | 178 | 192 |

Corollaries by shortening: A(26,6) >= 40960, A(25,6) >= 20480.

## New upper bounds
| Bound | Previous | New |
|---|---|---|
| A(17,4) <= | 3276 | 3272 |
| A(18,4) <= | 6551 | 6543 |
| A(19,4) <= | 13104 | 13085 |
| A(21,4) <= | 43688 | 43681 |
| A(25,4) <= | 599184 | 599174 |
| A(27,4) <= | 2396736 | 2396695 |
| A(17,6) <= | 340 | 339 |
| A(18,6) <= | 673 | 672 |
| A(22,6) <= | 6941 | 6939 |
| A(24,6) <= | 24106 | 24102 |
| A(26,6) <= | 84260 | 84243 |
| A(27,6) <= | 157285 | 157281 |
| A(28,6) <= | 291269 | 291265 |
| A(28,8) <= | 32151 | 32148 |
| A(21,10) <= | 47 | 46 |

Corollary by doubling: A(26,4) <= 1198348.

## Contents
- `paper/`: both papers (PDF and LaTeX).
- `lower/codes/`: the codes, one word per line; `lower/verify_codes.py`, `lower/verify_fast.c`: checkers;
  `lower/construction/`: groups, kernels, Z4 base code and coset lists; `lower/programs/`: search programs.
- `upper/certificates/UB_A_<n>_<d>_<U>/`: for each new upper bound, the certificate that no code of size U+1 exists
  (compressed JSON: run `unxz` first) and its exact checker (Python standard library only);
  `upper/programs/`: SDP builders, cutting-plane loop and certificate export.

Check a lower bound: `python3 lower/verify_codes.py lower/codes/A_24_10_144.txt 10`
Check an upper bound: `cd upper/certificates/UB_A_17_6_339 && python3 check_cert.py excl_n17_d6_M340_b1mel.json`

## Use of AI
AI tools were used in both works. Upper bounds: the authors do not guarantee the correctness of the results.
Lower bounds: the correctness does not depend on the AI tools, since each code is checked directly; the authors guarantee
the stated lower bounds.

## License
CC BY 4.0 (same as the Zenodo records).
