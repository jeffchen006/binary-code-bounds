// lnsx: word-level large-neighbourhood search for binary codes, exact re-fill by bitset max clique.
// Works in the odd form: an (n, M, d) code (d even) is punctured at coordinate p -> (n-1, M, D=d-1).
// cnt[x] = number of codewords at distance < D from x (uint16). Move: remove the k codewords nearest a centre,
// free set F = words with cnt == 0 inside the removed balls (contains the removed words), exact clique search for k+1
// in F (distance >= D); else plateau move: drop one removed word from F and search a k-clique (equal-size exchange).
// Output: improved codes are written in EVEN form (parity bit re-inserted at coordinate p) to <out>_<M>.txt.
#include <bits/stdc++.h>
using namespace std;
typedef uint32_t u32; typedef uint64_t u64;
static u64 rs = 88172645463325252ULL;
static inline u64 rnd() { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double urand() { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }

// ---------------- bitset max clique (coloring bound, BBMC-like) ----------------
struct Clq {
  int V = 0, W = 0; vector<u64> A;  // adjacency rows, W words each
  long long nodes = 0, limit = 0; bool aborted = false; int target = 0;
  vector<int> cur, best;
  u64* row(int v) { return &A[(size_t)v * W]; }
  void init(int v) { V = v; W = (v + 63) / 64; A.assign((size_t)V * W, 0); }
  void edge(int a, int b) { row(a)[b >> 6] |= 1ull << (b & 63); row(b)[a >> 6] |= 1ull << (a & 63); }
  void expand(vector<u64>& P, int pc) {
    if (aborted) return;
    if (++nodes > limit) { aborted = true; return; }
    // greedy colouring: order[] with colours col[]
    vector<int> order, col; order.reserve(pc); col.reserve(pc);
    vector<u64> U(P), Q(W);
    int c = 0, left = pc;
    while (left > 0) {
      c++; Q = U;
      for (int w = 0; w < W; w++) while (Q[w]) {
        int b = __builtin_ctzll(Q[w]); int v = w * 64 + b;
        Q[w] &= Q[w] - 1; U[w] &= ~(1ull << b);
        order.push_back(v); col.push_back(c); left--;
        const u64* r = row(v); for (int w2 = w; w2 < W; w2++) Q[w2] &= ~r[w2];
      }
    }
    vector<u64> NP(W);
    for (int i = (int)order.size() - 1; i >= 0; i--) {
      if ((int)cur.size() + col[i] <= (int)best.size()) return;
      int v = order[i]; const u64* r = row(v); int npc = 0;
      for (int w = 0; w < W; w++) { NP[w] = P[w] & r[w]; npc += __builtin_popcountll(NP[w]); }
      cur.push_back(v);
      if (npc == 0) { if (cur.size() > best.size()) best = cur; }
      else expand(NP, npc);
      cur.pop_back();
      if ((int)best.size() >= target) { aborted = true; return; }
      if (aborted) return;
      P[v >> 6] &= ~(1ull << (v & 63));
    }
  }
  // find a clique of size >= tgt (stop early); lb = known size not worth reporting (best starts empty, pruning uses tgt-1)
  int solve(int tgt, long long lim) {
    target = tgt; limit = lim; nodes = 0; aborted = false; cur.clear(); best.clear();
    // pruning threshold: pretend best has tgt-1 elements until a real one is found
    vector<int> fake(max(0, tgt - 1), -1); best = fake;
    vector<u64> P(W, 0); for (int v = 0; v < V; v++) P[v >> 6] |= 1ull << (v & 63);
    expand(P, V);
    if ((int)best.size() >= tgt && (best.empty() || best[0] >= 0)) return best.size();
    best.clear(); return 0;
  }
};

static int n0, n, D, pcoord;
static u32 NW;
static vector<u32> masks;
static vector<uint16_t> cnt;
static vector<u32> code;

static string wr_even(u32 x) {  // re-insert parity bit at pcoord
  string s; int par = __builtin_popcount(x) & 1; int j = 0;
  for (int c = 0; c < n0; c++) { if (c == pcoord) s += char('0' + par); else { s += char('0' + ((x >> j) & 1)); j++; } }
  return s;
}
static void addw(u32 c, int dlt) { for (u32 m : masks) cnt[c ^ m] += dlt; }
static int dmin_check() {
  int dm = 99; for (size_t i = 0; i < code.size(); i++) for (size_t j = i + 1; j < code.size(); j++) dm = min(dm, __builtin_popcount(code[i] ^ code[j]));
  return dm;
}
static void save(const string& out) {
  string f = out + "_" + to_string(code.size()) + ".txt";
  FILE* fp = fopen(f.c_str(), "w"); for (u32 c : code) fprintf(fp, "%s\n", wr_even(c).c_str()); fclose(fp);
  fprintf(stderr, "SAVED %s\n", f.c_str());
}

// low-count list L: words with 1 <= cnt <= K (codewords included, cnt 1). Lazy deletion via inL bitmap.
static int K = 0;
static vector<u32> L; static vector<u64> inL; static long long Lstale = 0;
static inline bool getb(const vector<u64>& b, u32 x) { return (b[x >> 6] >> (x & 63)) & 1; }
static inline void setb(vector<u64>& b, u32 x, bool v) { if (v) b[x >> 6] |= 1ull << (x & 63); else b[x >> 6] &= ~(1ull << (x & 63)); }
static void buildL() {
  L.clear(); inL.assign(NW / 64 + 1, 0); Lstale = 0;
  for (u32 x = 0; x < NW; x++) if (cnt[x] >= 1 && cnt[x] <= K) { L.push_back(x); setb(inL, x, 1); }
}
static void addwL(u32 c, int dlt) {  // update cnt and L membership
  for (u32 m : masks) {
    u32 x = c ^ m; cnt[x] += dlt; bool want = cnt[x] >= 1 && cnt[x] <= K; bool has = getb(inL, x);
    if (want && !has) { L.push_back(x); setb(inL, x, 1); } else if (!want && has) { setb(inL, x, 0); Lstale++; }
  }
}
static void compactL() {
  vector<u32> L2; L2.reserve(L.size()); vector<u64> seen(NW / 64 + 1, 0);
  for (u32 x : L) if (getb(inL, x) && !getb(seen, x)) { L2.push_back(x); setb(seen, x, 1); }
  L.swap(L2); Lstale = 0;
}

// ---- near-miss mode: conflict graph (odd form) and minimum vertex cover ----
static vector<pair<int,int>> conf_pairs(const vector<u32>& W, int Dd) {
  vector<pair<int,int>> P;
  if (W.size() <= 8000) {
    for (size_t i = 0; i < W.size(); i++) for (size_t j = i + 1; j < W.size(); j++) if (__builtin_popcount(W[i] ^ W[j]) < Dd) P.push_back({(int)i, (int)j});
  } else {  // uses masks/cnt (built for W)
    unordered_map<u32, vector<int>> pos; for (size_t i = 0; i < W.size(); i++) pos[W[i]].push_back(i);
    for (size_t i = 0; i < W.size(); i++) if (cnt[W[i]] > 1)
      for (u32 m : masks) { auto it = pos.find(W[i] ^ m); if (it == pos.end()) continue; for (int j : it->second) if (j > (int)i) P.push_back({(int)i, j}); }
  }
  return P;
}
struct VC {
  int V; vector<vector<int>> g; vector<char> in, alive; vector<int> best; int bestsz; long long nodes, lim; bool exact = true;
  int matchLB() { vector<char> used(V, 0); int m = 0; for (int v = 0; v < V; v++) if (alive[v] && !used[v]) for (int u : g[v]) if (alive[u] && !used[u]) { used[u] = used[v] = 1; m++; break; } return m; }
  void rec(int cur) {
    if (++nodes > lim) { exact = false; return; }
    int bv = -1, bd = 0;
    for (int v = 0; v < V; v++) if (alive[v]) { int dg = 0; for (int u : g[v]) if (alive[u]) dg++; if (dg > bd) { bd = dg; bv = v; } }
    if (bv < 0) { if (cur < bestsz) { bestsz = cur; best.clear(); for (int v = 0; v < V; v++) if (in[v]) best.push_back(v); } return; }
    if (cur + matchLB() >= bestsz) return;
    // branch 1: take N(bv)
    vector<int> nb; for (int u : g[bv]) if (alive[u]) nb.push_back(u);
    for (int u : nb) { alive[u] = 0; in[u] = 1; } alive[bv] = 0;
    rec(cur + nb.size());
    for (int u : nb) { alive[u] = 1; in[u] = 0; } alive[bv] = 1;
    // branch 2: take bv
    alive[bv] = 0; in[bv] = 1; rec(cur + 1); alive[bv] = 1; in[bv] = 0;
  }
  vector<int> greedy() { vector<char> al(V, 1); vector<int> c; while (true) { int bv = -1, bd = 0; for (int v = 0; v < V; v++) if (al[v]) { int dg = 0; for (int u : g[v]) if (al[u]) dg++; if (dg > bd) { bd = dg; bv = v; } } if (bv < 0) break; al[bv] = 0; c.push_back(bv); } return c; }
  vector<int> solve(long long L) { lim = L; nodes = 0; best = greedy(); bestsz = best.size(); in.assign(V, 0); alive.assign(V, 1); rec(0); return best; }
};
static vector<int> cover_of(int M, const vector<pair<int,int>>& P, bool exact, bool& isexact) {
  vector<int> id(M, -1), back; for (auto& e : P) for (int v : {e.first, e.second}) if (id[v] < 0) { id[v] = back.size(); back.push_back(v); }
  VC vc; vc.V = back.size(); vc.g.assign(vc.V, {}); for (auto& e : P) { vc.g[id[e.first]].push_back(id[e.second]); vc.g[id[e.second]].push_back(id[e.first]); }
  vector<int> c = exact ? vc.solve(2000000) : vc.greedy(); isexact = exact && vc.exact;
  vector<int> r; for (int v : c) r.push_back(back[v]); return r;
}

int main(int argc, char** argv) {
  string inf, out = "lns"; int d = 0; double tlim = 600; u64 seed = 1; int kmin = 4, kmax = 30; long long nlim = 2000000;
  int maxF = 12000; double pplat = 1.0; int recM = 0; int cmode = 2; int rmode = 0; long long kickT = 0; int kpmax = 0; int scratchN = 0; int Kmin = 0; int nm = 0; vector<u32> centers; double rep = 300;
  pcoord = -1;
  for (int i = 1; i < argc; i++) {
    string a = argv[i]; auto nx = [&]() { return string(argv[++i]); };
    if (a == "-i") inf = nx(); else if (a == "-o") out = nx(); else if (a == "-d") d = stoi(nx());
    else if (a == "-t") tlim = stod(nx()); else if (a == "-s") seed = stoull(nx()); else if (a == "-kmin") kmin = stoi(nx());
    else if (a == "-kmax") kmax = stoi(nx()); else if (a == "-nodes") nlim = stoll(nx()); else if (a == "-maxF") maxF = stoi(nx());
    else if (a == "-p") pcoord = stoi(nx()); else if (a == "-plat") pplat = stod(nx()); else if (a == "-rec") recM = stoi(nx());
    else if (a == "-cmode") cmode = stoi(nx()); else if (a == "-rmode") rmode = stoi(nx()); else if (a == "-K") K = stoi(nx());
    else if (a == "-kick") kickT = stoll(nx()); else if (a == "-kpmax") kpmax = stoi(nx());
    else if (a == "-scratch") scratchN = stoi(nx()); else if (a == "-rep") rep = stod(nx()); else if (a == "-Kmin") Kmin = stoi(nx()); else if (a == "-nm") nm = stoi(nx());
    else { fprintf(stderr, "unknown %s\n", a.c_str()); return 1; }
  }
  rs ^= seed * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 50; i++) rnd();
  vector<string> S; if (scratchN == 0) { ifstream in(inf); string s; while (in >> s) S.push_back(s); }
  n0 = scratchN ? scratchN : S[0].size(); n = n0 - 1; D = d - 1; NW = 1u << n;
  auto punct = [&](int p) { vector<u32> W; for (auto& s : S) { u32 v = 0; int j = 0; for (int c = 0; c < n0; c++) { if (c == p) continue; if (s[c] == '1') v |= 1u << j; j++; } W.push_back(v); } return W; };
  for (u32 e = 0; e < NW; e++) if (__builtin_popcount(e) < D) masks.push_back(e);
  cnt.assign(NW, 0);
  if (nm) {
    // near miss: choose the puncture coordinate with the smallest (greedy) vertex cover of the odd-form clash graph,
    // remove an exact minimum cover; the cover words become LNS centres. Parity re-insertion on output fixes d-1 pairs.
    int M = S.size(); vector<int> ps; if (pcoord >= 0 || M > 8000) ps.push_back(pcoord >= 0 ? pcoord : n0 - 1); else for (int p = 0; p < n0; p++) ps.push_back(p);
    int bp = -1; size_t bc = SIZE_MAX, bpairs = 0;
    for (int p : ps) {
      vector<u32> W = punct(p);
      if (M > 8000) { for (u32 c : W) addw(c, 1); }
      auto P = conf_pairs(W, D); bool ex; auto cv = cover_of(M, P, false, ex);
      if (M > 8000) { for (u32 c : W) addw(c, -1); }
      if (cv.size() < bc || (cv.size() == bc && P.size() < bpairs)) { bc = cv.size(); bp = p; bpairs = P.size(); }
    }
    pcoord = bp; vector<u32> W = punct(pcoord);
    if (M > 8000) { for (u32 c : W) addw(c, 1); }
    auto P = conf_pairs(W, D); bool ex; auto cv = cover_of(M, P, true, ex);
    if (M > 8000) { for (u32 c : W) addw(c, -1); }
    set<int> bw; for (auto& e : P) { bw.insert(e.first); bw.insert(e.second); }
    vector<char> rm(M, 0); for (int v : cv) { rm[v] = 1; centers.push_back(W[v]); }
    for (int i = 0; i < M; i++) if (!rm[i]) code.push_back(W[i]);
    fprintf(stderr, "NM p=%d bad_pairs=%zu bad_words=%zu cover=%zu cover_exact=%d clean=%zu\n", pcoord, P.size(), bw.size(), cv.size(), (int)ex, code.size());
  } else {
    if (pcoord < 0) pcoord = rnd() % n0;
    code = punct(pcoord);
  }
  if (dmin_check() < D) { fprintf(stderr, "input dmin < %d after puncture\n", D); return 1; }
  for (u32 c : code) addw(c, 1);
  int M0 = code.size();
  fprintf(stderr, "n=%d (odd form %d) D=%d M=%d p=%d |ball|=%zu seed=%llu K=%d k=%d..%d cmode=%d rmode=%d\n", n0, n, D, M0, pcoord, masks.size(),
          (unsigned long long)seed, K, kmin, kmax, cmode, rmode);
  if (scratchN) { u32 a = (u32)rnd() | 1, b = (u32)rnd(); for (u32 i = 0; i < NW; i++) { u32 x = (a * i + b) & (NW - 1); if (cnt[x] == 0) { code.push_back(x); addw(x, 1); } } }
  for (u32 x = 0; x < NW; x++) if (cnt[x] == 0) { code.push_back(x); addw(x, 1); }
  if ((int)code.size() > M0) { fprintf(stderr, "greedy fill -> %zu\n", code.size()); if ((int)code.size() > recM) save(out); }
  { vector<long long> h(65536, 0); for (u32 x = 0; x < NW; x++) h[cnt[x]]++;
    fprintf(stderr, "cnt histogram:"); for (int c = 0; c < 40; c++) if (h[c]) fprintf(stderr, " %d:%lld", c, h[c]); fprintf(stderr, "\n"); }
  if (tlim <= 0) return 0;
  if (K > 0) { buildL(); fprintf(stderr, "|L(K=%d)|=%zu\n", K, L.size()); }
  int best = code.size();
  auto t0 = chrono::steady_clock::now();
  auto el = [&]() { return chrono::duration<double>(chrono::steady_clock::now() - t0).count(); };
  long long it = 0, nplat = 0, nab = 0, nbig = 0, nimp = 0, nsame = 0; double sumF = 0, sumK = 0, sumC = 0;
  double lastrep = 0;
  vector<u32> Fw; Clq Q;
  // L3: best code (for revert) and kick state
  vector<u32> bestcode = code; long long stall = 0, nkick = 0, nrev = 0; int curtop = code.size();
  auto upd = [&](u32 w, int dl) { if (K > 0) addwL(w, dl); else addw(w, dl); };
  auto setcode = [&](const vector<u32>& tgt) {  // move current code to tgt by symmetric difference
    unordered_set<u32> T(tgt.begin(), tgt.end()), C(code.begin(), code.end());
    for (u32 c : code) if (!T.count(c)) upd(c, -1);
    for (u32 c : tgt) if (!C.count(c)) upd(c, 1);
    code = tgt; if (K > 0) compactL();
  };
  auto kick = [&]() {
    // face: codewords agreeing with a random codeword c0 on a random coordinate set S; isometry on the complement
    int M = code.size(); vector<int> P; u32 maskS = 0;
    for (int tries = 0; tries < 200; tries++) {
      u32 c0 = code[rnd() % M]; int s = 1 + rnd() % (n / 2); maskS = 0;
      while (__builtin_popcount(maskS) < s) maskS |= 1u << (rnd() % n);
      P.clear(); for (int i = 0; i < M; i++) if (((code[i] ^ c0) & maskS) == 0) P.push_back(i);
      if ((int)P.size() >= 2 && (int)P.size() <= (kpmax > 0 ? kpmax : M / 4)) break;
    }
    if (P.size() < 2) return;
    vector<int> fr; for (int c = 0; c < n; c++) if (!((maskS >> c) & 1)) fr.push_back(c);
    // mild isometry of the face: kw random transpositions of free coordinates + translation of weight <= kw
    vector<int> pf = fr; int kw = 1 + rnd() % 2;
    for (int q = 0; q < kw; q++) swap(pf[rnd() % pf.size()], pf[rnd() % pf.size()]);
    u32 t = 0; int tw = rnd() % (kw + 1); for (int q = 0; q < tw; q++) t |= 1u << fr[rnd() % fr.size()];
    vector<u32> Pw; for (int i : P) Pw.push_back(code[i]);
    sort(P.begin(), P.end(), greater<int>()); for (int i : P) { code[i] = code.back(); code.pop_back(); }
    for (u32 w : Pw) upd(w, -1);
    vector<u32> img;
    for (u32 w : Pw) { u32 y = w & maskS; for (size_t j = 0; j < fr.size(); j++) if ((w >> fr[j]) & 1) y |= 1u << pf[j]; img.push_back(y ^ t); }
    // force-insert the image; remove outside codewords that clash with it (repair)
    sort(img.begin(), img.end()); img.erase(unique(img.begin(), img.end()), img.end());
    int ins = img.size(), lost = 0;
    vector<char> kill(code.size(), 0);
    for (size_t i = 0; i < code.size(); i++) for (u32 y : img) if (__builtin_popcount(code[i] ^ y) < D) { kill[i] = 1; break; }
    vector<u32> kept; for (size_t i = 0; i < code.size(); i++) { if (kill[i]) { upd(code[i], -1); lost++; } else kept.push_back(code[i]); }
    vector<u32> gone; for (size_t i = 0; i < code.size(); i++) if (kill[i]) gone.push_back(code[i]);
    code = kept;
    for (u32 y : img) { code.push_back(y); upd(y, 1); }
    // greedy fill in the vacated balls, random start, original words of P excluded
    unordered_set<u32> tabu(Pw.begin(), Pw.end());
    int fill = 0;
    vector<u32> zs = Pw; zs.insert(zs.end(), gone.begin(), gone.end());
    for (int i = zs.size() - 1; i > 0; i--) swap(zs[i], zs[rnd() % (i + 1)]);
    size_t off = rnd() % masks.size();
    for (u32 w : zs) for (size_t q = 0; q < masks.size(); q++) { u32 x = w ^ masks[(q + off) % masks.size()]; if (cnt[x] == 0 && !tabu.count(x)) { code.push_back(x); upd(x, 1); fill++; } }
    if (K > 0) compactL();
    nkick++;
    if (nkick <= 20 || nkick % 100 == 0) fprintf(stderr, "kick %lld: face |P|=%zu s=%d -> inserted %d, lost %d, fill %d, M=%zu\n", nkick, Pw.size(), __builtin_popcount(maskS), ins, lost, fill, code.size());
  };
  while (el() < tlim) {
    it++;
    if (kickT > 0 && ++stall > kickT) {
      stall = 0;
      if ((int)code.size() < best) { setcode(bestcode); nrev++; }
      else kick();
      curtop = code.size();
    }
    int M = code.size();
    int k = kmin + rnd() % (kmax - kmin + 1); if (k > M) k = M;
    // centre: 0 random word, 1 codeword + few flips, 2 mix, 3 word with cnt <= k
    u32 z;
    int cm = cmode == 2 ? (int)(rnd() % 3) : cmode;
    if (!centers.empty() && (rnd() & 1)) cm = 4;
    if (cm == 1) { z = code[rnd() % M]; int f = 1 + rnd() % (D / 2 + 1); for (int q = 0; q < f; q++) z ^= 1u << (rnd() % n); }
    else if (cm == 2 || cm == 3) { int tries = 0; do { z = rnd() & (NW - 1); tries++; } while (cnt[z] > k && tries < 2000); }
    else if (cm == 4) { z = centers[rnd() % centers.size()]; int f = rnd() % 3; for (int q = 0; q < f; q++) z ^= 1u << (rnd() % n); }
    else z = rnd() & (NW - 1);
    // metric: 0 Hamming; 1 projected on a random coordinate subset (sub-cube face); 2 mix; 3 random subset
    int rm = rmode == 2 ? (int)(rnd() % 2) : rmode;
    u32 pm = NW - 1;
    if (rm == 1) { int s = 2 + rnd() % (n - 2); pm = 0; while (__builtin_popcount(pm) < s) pm |= 1u << (rnd() % n); }
    vector<pair<u32, int>> ds(M);
    for (int i = 0; i < M; i++) ds[i] = {(rm == 3 ? 0u : (u32)__builtin_popcount((code[i] ^ z) & pm)) << 24 | (u32)(rnd() & 0xFFFFFF), i};
    nth_element(ds.begin(), ds.begin() + k, ds.end());
    vector<int> rem(k); for (int q = 0; q < k; q++) rem[q] = ds[q].second;
    vector<u32> R(k); for (int q = 0; q < k; q++) R[q] = code[rem[q]];
    Fw.clear(); bool big = false;
    if (K > 0) {
      int Km = (Kmin > 0 && Kmin < K) ? Kmin + (int)(rnd() % (K - Kmin + 1)) : K;
      for (u32 x : L) {
        if (!getb(inL, x)) continue;
        int c = cnt[x]; if (c > k || c > Km) continue;
        int h = 0; for (u32 r : R) if (__builtin_popcount(x ^ r) < D) { if (++h > c) break; }
        if (h == c) { if (cnt[x] & 0x8000) continue; cnt[x] |= 0x8000; Fw.push_back(x); }
      }
      for (u32 x : Fw) cnt[x] &= 0x7FFF;
      if ((int)Fw.size() > maxF) big = true;
    } else {
      for (u32 r : R) addw(r, -1);
      for (u32 r : R) { for (u32 m : masks) { u32 x = r ^ m; if (cnt[x] == 0) { cnt[x] = 0xFFFF; Fw.push_back(x); } } if ((int)Fw.size() > maxF) { big = true; break; } }
      for (u32 x : Fw) cnt[x] = 0;
      for (u32 r : R) addw(r, 1);
    }
    if (big) { nbig++; continue; }
    int F = Fw.size(); sumF += F; sumK += k;
    vector<int> deg(F, 0);
    vector<char> adj((size_t)F * F, 0);
    for (int a = 0; a < F; a++) for (int b = a + 1; b < F; b++) if (__builtin_popcount(Fw[a] ^ Fw[b]) >= D) { adj[(size_t)a * F + b] = adj[(size_t)b * F + a] = 1; deg[a]++; deg[b]++; }
    if (getenv("LNSDBG") && it < 40) { long long e = 0; for (int a = 0; a < F; a++) e += deg[a]; fprintf(stderr, "dbg k=%d F=%d compat density %.3f avg conflict deg %.1f\n", k, F, F > 1 ? (double)e / F / (F - 1) : 0, F - 1 - (double)e / max(1, F)); }
    vector<int> ord(F); iota(ord.begin(), ord.end(), 0);
    vector<u32> key(F); for (int a = 0; a < F; a++) key[a] = rnd() & 1023;
    sort(ord.begin(), ord.end(), [&](int a, int b) { return deg[a] != deg[b] ? deg[a] > deg[b] : key[a] < key[b]; });
    auto build = [&](int skip, int tgt) {
      // (tgt-1)-core peeling: a vertex of a tgt-clique has >= tgt-1 neighbours inside the set
      vector<char> alive(F, 1); if (skip >= 0) alive[skip] = 0;
      vector<int> dg(F, 0);
      for (int a = 0; a < F; a++) if (alive[a]) for (int b = 0; b < F; b++) if (alive[b] && adj[(size_t)a * F + b]) dg[a]++;
      vector<int> st; for (int a = 0; a < F; a++) if (alive[a] && dg[a] < tgt - 1) { alive[a] = 0; st.push_back(a); }
      while (!st.empty()) { int a = st.back(); st.pop_back(); for (int b = 0; b < F; b++) if (alive[b] && adj[(size_t)a * F + b] && --dg[b] < tgt - 1) { alive[b] = 0; st.push_back(b); } }
      vector<int> vs; for (int a : ord) if (alive[a]) vs.push_back(a);
      if ((int)vs.size() < tgt) vs.clear();
      Q.init(vs.size());
      for (size_t i = 0; i < vs.size(); i++) for (size_t j = i + 1; j < vs.size(); j++) if (adj[(size_t)vs[i] * F + vs[j]]) Q.edge(i, j);
      return vs;
    };
    vector<u32> add;
    vector<int> vs;
    if (F > k) {  // F == k: only R itself
      vs = build(-1, k + 1); sumC += vs.size();
      int got = Q.solve(k + 1, nlim);
      if (got >= k + 1) { for (int v : Q.best) add.push_back(Fw[vs[v]]); nimp++; }
      else if (Q.aborted) nab++;
      if (add.empty() && urand() < pplat) {
        int drop = -1; u32 r0 = R[rnd() % k];
        for (int a = 0; a < F; a++) if (Fw[a] == r0) drop = a;
        vs = build(drop, k);
        int g2 = Q.solve(k, nlim / 4 + 1);
        if (g2 >= k) { for (int v : Q.best) add.push_back(Fw[vs[v]]); nplat++; }
        else if (Q.aborted) nab++;
      }
    } else nsame++;
    if (!add.empty()) {
      // symmetric difference update
      unordered_set<u32> As(add.begin(), add.end()), Rs(R.begin(), R.end());
      sort(rem.begin(), rem.end(), greater<int>());
      for (int i : rem) { code[i] = code.back(); code.pop_back(); }
      for (u32 r : R) if (!As.count(r)) { if (K > 0) addwL(r, -1); else addw(r, -1); }
      for (u32 x : add) { code.push_back(x); if (!Rs.count(x)) { if (K > 0) addwL(x, 1); else addw(x, 1); } }
      if (K > 0 && Lstale > (long long)L.size() / 2) compactL();
      if ((int)code.size() > curtop) { curtop = code.size(); stall = 0; }
      if ((int)code.size() >= best && kickT > 0) bestcode = code;
      if ((int)code.size() > best) {
        best = code.size();
        fprintf(stderr, "IMPROVED M=%d it=%lld t=%.0f k=%d F=%d\n", best, it, el(), k, F);
        int dm = dmin_check();
        if (dm < D) { fprintf(stderr, "BUG dmin %d\n", dm); return 2; }
        if (best > recM) save(out);
      }
    }
    if (el() - lastrep > rep) {
      lastrep = el();
      fprintf(stderr, "t=%.0f it=%lld M=%zu best=%d imp=%lld plat=%lld abort=%lld big=%lld onlyR=%lld avgK=%.1f avgF=%.1f core=%.1f |L|=%zu kicks=%lld rev=%lld\n", el(), it, code.size(), best, nimp, nplat, nab, nbig, nsame,
              sumK / max(1LL, it - nbig), sumF / max(1LL, it - nbig), sumC / max(1LL, it - nbig), L.size(), nkick, nrev);
    }
  }
  fprintf(stderr, "END t=%.0f it=%lld M=%zu best=%d imp=%lld plat=%lld abort=%lld big=%lld onlyR=%lld avgK=%.1f avgF=%.1f\n", el(), it, code.size(), best, nimp, nplat, nab, nbig, nsame,
          sumK / max(1LL, it - nbig), sumF / max(1LL, it - nbig));
  save(out + "_final");
  return 0;
}
