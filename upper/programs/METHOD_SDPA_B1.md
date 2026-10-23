# Schrijver SDP + B1 rows with SDPA-DD: recipe and why CSDP missed it (sdp, 2026-10-22)

## Results (all premise-free, exact rational certificates, stdlib checker)
| entry | table | new | base (plain Schrijver) | excluded M | smallest margin |
|---|---|---|---|---|---|
| A(17,6) | 340 | 339 | 341.33 (+ MEL row) | 341, 340 | 0.42 |
| A(26,6) | 84260 | 84243 | 84260.57 | 84260..84244 | 0.170 |
| A(27,6) | 157285 | 157281 | 157286.40 | 157286..157282 | 0.056 |
| A(28,6) | 291269 | 291265 | 291271.11 | 291271..291266 | 0.517 |
| A(28,8) | 32151 | 32148 | 32151.80 | 32151..32149 | 0.122 |

## Recipe
1. Base: plain Schrijver 3-point SDP (no table premise) -> |C| <= floor(value).
2. For each M from the base down: fix |C| = M; add B1 rows (novaB, METHOD_B1.md): every block subspace with
   rho = M mod 2^r != 0, plus general V found by separation at the current optimum (cp.separate, r <= 10).
   Solve with SDPA-DD (opt/sdpa-gmp/sdpa_dd2, double-double), equilibrated blocks. Repeat while cuts are found.
3. Stop at the first M with value exactly M (feasible). Certificates (make_cert_b1mel.py --sdpa): integer V with
   Y = V V^T / 4^80 (PSD by construction), integer z >= 0; checked in exact fractions by check_cert.py.
Tools: cplin.py (USE_SDPA=1, LINOPTS=none), jobs/descend.sh, jobs/certs.sh, make_cert_b1mel.py.

## Why it was missed before
novaB / novaA ran the same rows with CSDP. At n >= 26, CSDP stops with "reduced accuracy" and reports a primal value
at or above M ('feasible'), e.g. (28,8) M=32151: CSDP 32151.03 vs exact certificate 32148.77. novaA's (28,8) runs gave
32311 (CSDP failure). With SDPA-DD the dual is accurate enough for an exact certificate with margin.

## Where it does not work
Entries whose table value is below floor(Schrijver) (GMS 4-point or Best/Kim-Toan type bounds): Schrijver + B1 stays at
or above the table value (scan list in LOG.md). GMS 4-point + B1 (gmsb1.py) was tried at (18,6) M=672 and (17,6) M=339:
no exclusion (~673 and 340.5).

## Controls
runs/controls.log: (24,8,4096) (16,6,256) (22,8,1024) (21,6,3072) (17,4,3072) (27,6,81920) (28,6,139264)
(28,8,16384) (26,6,40960): sizes of real codes, none excluded.
