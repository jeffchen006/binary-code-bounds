/* Parallel verifier (OpenMP). Same input as verify: one 0/1 word per line, n <= 64.
 * Usage: verify_fast <file> [-d D] [-t THREADS] [-p MAXPRINT]
 *   prints n, M, dmin, distance distribution, weight distribution.
 *   With -d D: also lists up to MAXPRINT (default 20) bad pairs (distance < D) as line numbers,
 *   and the number of words involved in bad pairs (for debugging near-codes).
 * Exit 0: ok (and dmin >= D if -d given); 1: duplicates, length mismatch, or dmin < D. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
int main(int argc, char **argv) {
  if (argc < 2) { fprintf(stderr, "usage: verify_fast file [-d D] [-t T] [-p P]\n"); return 2; }
  int D = 0, T = 0, P = 20;
  for (int i = 2; i + 1 < argc; i += 2) {
    if (!strcmp(argv[i], "-d")) D = atoi(argv[i+1]);
    else if (!strcmp(argv[i], "-t")) T = atoi(argv[i+1]);
    else if (!strcmp(argv[i], "-p")) P = atoi(argv[i+1]);
  }
  if (T > 0) omp_set_num_threads(T);
  FILE *f = fopen(argv[1], "r"); if (!f) { perror("open"); return 2; }
  size_t cap = 1 << 16, M = 0; uint64_t *w = malloc(cap * 8); long *ln = malloc(cap * sizeof(long));
  char line[4096]; int n = -1; long lineno = 0;
  while (fgets(line, sizeof line, f)) {
    lineno++;
    if (line[0] == '#') continue;
    uint64_t x = 0; int len = 0;
    for (char *p = line; *p; p++) {
      if (*p == '0' || *p == '1') { if (len >= 64) { fprintf(stderr, "n>64\n"); return 2; } x = (x << 1) | (uint64_t)(*p - '0'); len++; }
      else if (*p != ' ' && *p != '\t' && *p != '\n' && *p != '\r' && *p != ',') { fprintf(stderr, "bad char line %ld\n", lineno); return 2; }
    }
    if (len == 0) continue;
    if (n < 0) n = len; else if (len != n) { fprintf(stderr, "length mismatch at line %ld\n", lineno); return 1; }
    if (M == cap) { cap *= 2; w = realloc(w, cap * 8); ln = realloc(ln, cap * sizeof(long)); }
    ln[M] = lineno; w[M++] = x;
  }
  int nt = omp_get_max_threads();
  long long *dist = calloc((size_t)nt * 65, sizeof(long long));
  char *bad = D ? calloc(M, 1) : NULL;
  long long nbad = 0; int printed = 0;
  #pragma omp parallel for schedule(dynamic, 64) reduction(+:nbad)
  for (size_t i = 0; i < M; i++) {
    long long *dd = dist + (size_t)omp_get_thread_num() * 65;
    uint64_t a = w[i];
    for (size_t j = i + 1; j < M; j++) {
      int d = __builtin_popcountll(a ^ w[j]); dd[d]++;
      if (d < D) { nbad++; bad[i] = bad[j] = 1;
        if (printed < P) {
          #pragma omp critical
          { if (printed < P) { printf("bad pair: line %ld, line %ld, distance %d\n", ln[i], ln[j], d); printed++; } }
        }
      }
    }
  }
  long long tot[65] = {0}; for (int t = 0; t < nt; t++) for (int d = 0; d <= 64; d++) tot[d] += dist[t * 65 + d];
  int dmin = 65; for (int d = 0; d <= 64; d++) if (tot[d]) { dmin = d; break; }
  long long wt[65] = {0}; for (size_t i = 0; i < M; i++) wt[__builtin_popcountll(w[i])]++;
  printf("n=%d M=%zu dmin=%d\n", n, M, M > 1 ? dmin : -1);
  printf("distance distribution:"); for (int d = 0; d <= n; d++) if (tot[d]) printf(" %d:%lld", d, tot[d]); printf("\n");
  printf("weight distribution:"); for (int d = 0; d <= n; d++) if (wt[d]) printf(" %d:%lld", d, wt[d]); printf("\n");
  if (tot[0]) { printf("DUPLICATE words: %lld pairs\n", tot[0]); return 1; }
  if (D) {
    size_t nb = 0; for (size_t i = 0; i < M; i++) nb += bad[i];
    printf("pairs with distance < %d: %lld; words involved: %zu\n", D, nbad, nb);
    return nbad ? 1 : 0;
  }
  return 0;
}
