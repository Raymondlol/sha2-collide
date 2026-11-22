// Step A: Re-solve t=9..16 (or t=10..16) of 2026/1120 Table 14 with t=0..8 fixed,
// shrinking popcount(E13 ^ E13') from 34 bits down to ~5-14 bits to unlock 2^{38+} (W14,W15) freedom!
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <array>
#include <chrono>
#include "cadical.hpp"

using u64 = uint64_t;
static inline u64 ror64(u64 x, int n) { return (x >> n) | (x << (64 - n)); }
static inline u64 S0(u64 x) { return ror64(x,28)^ror64(x,34)^ror64(x,39); }
static inline u64 S1(u64 x) { return ror64(x,14)^ror64(x,18)^ror64(x,41); }
static inline u64 s0(u64 x) { return ror64(x,1)^ror64(x,8)^(x>>7); }
static inline u64 s1(u64 x) { return ror64(x,19)^ror64(x,61)^(x>>6); }
static inline u64 Ch(u64 x, u64 y, u64 z) { return (x&y)^((~x)&z); }
static inline u64 Maj(u64 x, u64 y, u64 z) { return (x&y)^(x&z)^(y&z); }

static const u64 K512[36] = {
    0x428a2f98d728ae22ULL, 0x7137449123ef65cdULL, 0xb5c0fbcfec4d3b2fULL, 0xe9b5dba58189dbbcULL,
    0x3956c25bf348b538ULL, 0x59f111f1b605d019ULL, 0x923f82a4af194f9bULL, 0xab1c5ed5da6d8118ULL,
    0xd807aa98a3030242ULL, 0x12835b0145706fbeULL, 0x243185be4ee4b28cULL, 0x550c7dc3d5ffb4e2ULL,
    0x72be5d74f27b896fULL, 0x80deb1fe3b1696b1ULL, 0x9bdc06a725c71235ULL, 0xc19bf174cf692694ULL,
    0xe49b69c19ef14ad2ULL, 0xefbe4786384f25e3ULL, 0x0fc19dc68b8cd5b5ULL, 0x240ca1cc77ac9c65ULL,
    0x2de92c6f592b0275ULL, 0x4a7484aa6ea6e483ULL, 0x5cb0a9dcbd41fbd4ULL, 0x76f988da831153b5ULL,
    0x983e5152ee66dfabULL, 0xa831c66d2db43210ULL, 0xb00327c898fb213fULL, 0xbf597fc7beef0ee4ULL,
    0xc6e00bf33da88fc2ULL, 0xd5a79147930aa725ULL, 0x06ca6351e003826fULL, 0x142929670a0e6e70ULL,
    0x27b70a8546d22ffcULL, 0x2e1b21385c26c926ULL, 0x4d2c6dfc5ac42aedULL, 0x53380d139d95b3dfULL
};

static const u64 REF_A1[17] = {
    0,0,0,0,0,
    0x4608d7572ef3642aULL, 0xb662caca8de8796eULL, 0xa96aba2137015543ULL, 0xc0720f9042caf688ULL,
    0xfcc893e8654115e2ULL, 0xbc157d076390492cULL, 0xe3698ee4fe1319afULL, 0x863aebeaa5768210ULL,
    0x44f0c7b6bd1e6710ULL, 0x123803f0e014e59aULL, 0x7b5a8e1c3bb8af4bULL, 0xa6e27582bc9f52b6ULL
};
static const u64 REF_A2[17] = {
    0,0,0,0,0,
    0x4e08d7572ef3642aULL, 0xb26bc8ca0fec796eULL, 0xa96aba2137015543ULL, 0xc0720f9042caf688ULL,
    0xfcc893e8654115e2ULL, 0xbc157d656390db38ULL, 0xeb698ee4fe1319afULL, 0x863ad1ea27828210ULL,
    0x46f8cfb6bd1e6710ULL, 0x123803f0e014e59aULL, 0x735a8e1c3bb8af4bULL, 0xa6e27582bc9f52b6ULL
};
static const u64 REF_E1[17] = {
    0,0,0,0,0,
    0xbfa6e9c2f5238cdbULL, 0x0dfb919806ae7544ULL, 0x457b9fe677c1632fULL, 0x4054cf54917b012bULL,
    0x85d69ffc65af65c3ULL, 0xec8212b4d3f4e17dULL, 0xa15c89731a29737dULL, 0x3319e0efefabffffULL,
    0x6fffff2fd821023fULL, 0xf3f8dd0307d580adULL, 0x80001002018a91f5ULL, 0x4d912b7cab6ff654ULL
};
static const u64 REF_E2[17] = {
    0,0,0,0,0,
    0xc7a6e9c2f5238cdbULL, 0x02048f9806a27544ULL, 0x453ea9e6375395bbULL, 0x3c54cf54917b012bULL,
    0x8dce9ffc65af65c3ULL, 0xe88b111255f97389ULL, 0xed5584323228f355ULL, 0x3319c6efefa80003ULL,
    0x70000132579d023fULL, 0xf3f8dd6125d484a9ULL, 0x80001002018a91f5ULL, 0x4d91097cab6bf654ULL
};
static const u64 REF_W1[17] = {
    0x230a4930d38229e3ULL, 0x2b3c63372261608bULL, 0xf833951c32f2b6dcULL, 0x6f4f49447a93b0ccULL,
    0xf841e91df9c8d69eULL, 0xa4710d44716286a3ULL, 0xe3310abbec7a9816ULL, 0x66b0241d0bfd564cULL,
    0x5395fbd4a70df9f3ULL, 0xb0f8774e493b43bcULL, 0x13287fa0a341d877ULL, 0xf0d19e5ae32f6abfULL,
    0xf26d943cf8041873ULL, 0x06f2b17bb37eb2b7ULL, 0x2cccf19ec39c66ddULL, 0xf284ae681c043be9ULL,
    0xdcf55067def22084ULL
};
static const u64 REF_W2[17] = {
    0x230a4930d38229e3ULL, 0x2b3c63372261608bULL, 0xf833951c32f2b6dcULL, 0x6f4f49447a93b0ccULL,
    0xf841e91df9c8d69eULL, 0xac710d44716286a3ULL, 0xe7390abbec7a9816ULL, 0x6eb0241d0bfd564cULL,
    0x73c5db55a70df9f3ULL, 0xf0d8764e493b43bcULL, 0x13287fa0a341d877ULL, 0xf0d19e5ae32f6abfULL,
    0xf26d943cf8041873ULL, 0x02fab17bb37eb2b7ULL, 0x24ccf19ec39c66ddULL, 0xf284ae681c043be9ULL,
    0xdcf55067def22084ULL
};

struct Builder {
    CaDiCaL::Solver &S;
    int nv = 0;
    int ncl = 0;
    int zero_lit, one_lit;

    Builder(CaDiCaL::Solver &s) : S(s) {
        one_lit = new_var();
        add({one_lit});
        zero_lit = -one_lit;
    }
    int new_var() { return ++nv; }
    std::array<int, 64> new_word() {
        std::array<int, 64> w;
        for (int i = 0; i < 64; i++) w[i] = new_var();
        return w;
    }
    std::array<int, 64> const_word(u64 val) {
        std::array<int, 64> w;
        for (int i = 0; i < 64; i++) w[i] = ((val >> i) & 1) ? one_lit : zero_lit;
        return w;
    }
    void add(std::initializer_list<int> lits) {
        for (int l : lits) S.add(l);
        S.add(0);
        ncl++;
    }
    void fix_word(const std::array<int,64> &w, u64 val) {
        for (int i = 0; i < 64; i++) add({((val >> i) & 1) ? w[i] : -w[i]});
    }
    void eq(int a, int b) { add({-a, b}); add({a, -b}); }
    void eq_word(const std::array<int,64> &a, const std::array<int,64> &b) {
        for (int i = 0; i < 64; i++) eq(a[i], b[i]);
    }
    int xor2(int a, int b) {
        int z = new_var();
        add({-a, -b, -z}); add({a, b, -z}); add({-a, b, z}); add({a, -b, z});
        return z;
    }
    int xor3(int a, int b, int c) {
        int z = new_var();
        add({-a, -b, -c,  z}); add({-a,  b,  c,  z}); add({ a, -b,  c,  z}); add({ a,  b, -c,  z});
        add({-a, -b,  c, -z}); add({-a,  b, -c, -z}); add({ a, -b, -c, -z}); add({ a,  b,  c, -z});
        return z;
    }
    int maj3(int a, int b, int c) {
        int z = new_var();
        add({-a, -b, z}); add({-a, -c, z}); add({-b, -c, z});
        add({ a,  b,-z}); add({ a,  c,-z}); add({ b,  c,-z});
        return z;
    }
    int ch3(int x, int y, int z) {
        int r = new_var();
        add({-x, -y,  r}); add({-x,  y, -r});
        add({ x, -z,  r}); add({ x,  z, -r});
        add({-y, -z,  r}); add({ y,  z, -r});
        return r;
    }
    std::array<int,64> w_S0(const std::array<int,64> &x) {
        std::array<int,64> r;
        for (int i = 0; i < 64; i++) r[i] = xor3(x[(i+28)&63], x[(i+34)&63], x[(i+39)&63]);
        return r;
    }
    std::array<int,64> w_S1(const std::array<int,64> &x) {
        std::array<int,64> r;
        for (int i = 0; i < 64; i++) r[i] = xor3(x[(i+14)&63], x[(i+18)&63], x[(i+41)&63]);
        return r;
    }
    std::array<int,64> w_s0(const std::array<int,64> &x) {
        std::array<int,64> r;
        for (int i = 0; i < 64; i++) {
            int a = x[(i+1)&63], b = x[(i+8)&63];
            if (i + 7 < 64) r[i] = xor3(a, b, x[i+7]);
            else r[i] = xor2(a, b);
        }
        return r;
    }
    std::array<int,64> w_s1(const std::array<int,64> &x) {
        std::array<int,64> r;
        for (int i = 0; i < 64; i++) {
            int a = x[(i+19)&63], b = x[(i+61)&63];
            if (i + 6 < 64) r[i] = xor3(a, b, x[i+6]);
            else r[i] = xor2(a, b);
        }
        return r;
    }
    std::array<int,64> w_Ch(const std::array<int,64> &x, const std::array<int,64> &y, const std::array<int,64> &z) {
        std::array<int,64> r;
        for (int i = 0; i < 64; i++) r[i] = ch3(x[i], y[i], z[i]);
        return r;
    }
    std::array<int,64> w_Maj(const std::array<int,64> &x, const std::array<int,64> &y, const std::array<int,64> &z) {
        std::array<int,64> r;
        for (int i = 0; i < 64; i++) r[i] = maj3(x[i], y[i], z[i]);
        return r;
    }
    std::array<int,64> add2(const std::array<int,64> &a, const std::array<int,64> &b) {
        std::array<int,64> s;
        int c = zero_lit;
        for (int i = 0; i < 64; i++) {
            if (i == 0) {
                s[i] = xor2(a[i], b[i]);
                int nc = new_var();
                add({-a[i], -b[i], nc}); add({a[i], -nc}); add({b[i], -nc});
                c = nc;
            } else {
                s[i] = xor3(a[i], b[i], c);
                if (i < 63) c = maj3(a[i], b[i], c);
            }
        }
        return s;
    }
    void csa3(const std::array<int,64> &a, const std::array<int,64> &b, const std::array<int,64> &c,
              std::array<int,64> &u, std::array<int,64> &v) {
        v[0] = zero_lit;
        for (int i = 0; i < 64; i++) {
            u[i] = xor3(a[i], b[i], c[i]);
            if (i < 63) v[i+1] = maj3(a[i], b[i], c[i]);
        }
    }
    std::array<int,64> add_many(std::vector<std::array<int,64>> ws) {
        while (ws.size() > 2) {
            std::vector<std::array<int,64>> nxt;
            size_t i = 0;
            for (; i + 2 < ws.size(); i += 3) {
                std::array<int,64> u, v;
                csa3(ws[i], ws[i+1], ws[i+2], u, v);
                nxt.push_back(u); nxt.push_back(v);
            }
            for (; i < ws.size(); i++) nxt.push_back(ws[i]);
            ws = std::move(nxt);
        }
        if (ws.size() == 1) return ws[0];
        return add2(ws[0], ws[1]);
    }
    std::vector<int> build_seq_counter(const std::vector<int> &xs, int maxK) {
        int n = (int)xs.size();
        std::vector<std::vector<int>> s(n, std::vector<int>(maxK, 0));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < maxK; j++)
                s[i][j] = new_var();
        add({-xs[0], s[0][0]});
        for (int j = 1; j < maxK; j++) add({-s[0][j]});
        for (int i = 1; i < n; i++) {
            add({-xs[i], s[i][0]});
            add({-s[i-1][0], s[i][0]});
            for (int j = 1; j < maxK; j++) {
                add({-s[i-1][j], s[i][j]});
                add({-xs[i], -s[i-1][j-1], s[i][j]});
            }
        }
        return s[n-1];
    }
    u64 val_word(const std::array<int,64> &w) {
        u64 r = 0;
        for (int i = 0; i < 64; i++) if (S.val(w[i]) > 0) r |= (1ULL << i);
        return r;
    }
};

static void verify_and_report(const u64 W1_sol[16], const u64 W2_sol[16], const char *tag) {
    u64 A1[17], A2[17], E1[17], E2[17];
    for (int t = 5; t <= 8; t++) {
        A1[t] = REF_A1[t]; A2[t] = REF_A2[t];
        E1[t] = REF_E1[t]; E2[t] = REF_E2[t];
    }
    for (int t = 9; t <= 16; t++) {
        u64 w1 = (t < 16) ? W1_sol[t] : (s1(W1_sol[14]) + W1_sol[9] + s0(W1_sol[1]) + W1_sol[0]);
        u64 w2 = (t < 16) ? W2_sol[t] : (s1(W2_sol[14]) + W2_sol[9] + s0(W2_sol[1]) + W2_sol[0]);
        u64 t1_1 = E1[t-4] + S1(E1[t-1]) + Ch(E1[t-1], E1[t-2], E1[t-3]) + K512[t] + w1;
        u64 t2_1 = S0(A1[t-1]) + Maj(A1[t-1], A1[t-2], A1[t-3]);
        A1[t] = t1_1 + t2_1; E1[t] = A1[t-4] + t1_1;
        u64 t1_2 = E2[t-4] + S1(E2[t-1]) + Ch(E2[t-1], E2[t-2], E2[t-3]) + K512[t] + w2;
        u64 t2_2 = S0(A2[t-1]) + Maj(A2[t-1], A2[t-2], A2[t-3]);
        A2[t] = t1_2 + t2_2; E2[t] = A2[t-4] + t1_2;
    }
    bool ok = true;
    auto chk = [&](const char *name, u64 got, u64 exp) {
        if (got != exp) { printf("  [FAIL %s] got=%016llx exp=%016llx\n", name, (unsigned long long)got, (unsigned long long)exp); ok = false; }
    };
    chk("dW9",   W2_sol[9]-W1_sol[9],   REF_W2[9]-REF_W1[9]);
    chk("ds0W9", s0(W2_sol[9])-s0(W1_sol[9]), s0(REF_W2[9])-s0(REF_W1[9]));
    chk("dW10",  W2_sol[10]-W1_sol[10], 0);
    chk("dW11",  W2_sol[11]-W1_sol[11], 0);
    chk("dW12",  W2_sol[12]-W1_sol[12], 0);
    chk("dW13",  W2_sol[13]-W1_sol[13], REF_W2[13]-REF_W1[13]);
    chk("ds0W13",s0(W2_sol[13])-s0(W1_sol[13]), s0(REF_W2[13])-s0(REF_W1[13]));
    chk("dW14",  W2_sol[14]-W1_sol[14], REF_W2[14]-REF_W1[14]);
    chk("ds0W14",s0(W2_sol[14])-s0(W1_sol[14]), s0(REF_W2[14])-s0(REF_W1[14]));
    chk("ds1W14",s1(W2_sol[14])-s1(W1_sol[14]), s1(REF_W2[14])-s1(REF_W1[14]));
    chk("dW15",  W2_sol[15]-W1_sol[15], 0);
    chk("dA13",  A2[13]-A1[13], REF_A2[13]-REF_A1[13]);
    chk("dE13",  E2[13]-E1[13], REF_E2[13]-REF_E1[13]);
    chk("dA14",  A2[14]-A1[14], 0);
    chk("dE14",  E2[14]-E1[14], REF_E2[14]-REF_E1[14]);
    chk("dA15",  A2[15]-A1[15], REF_A2[15]-REF_A1[15]);
    chk("dE15",  E2[15]-E1[15], 0);
    chk("dA16",  A2[16]-A1[16], 0);
    chk("dE16",  E2[16]-E1[16], REF_E2[16]-REF_E1[16]);

    u64 xE13 = E1[13] ^ E2[13];
    u64 xE14 = E1[14] ^ E2[14];
    u64 xA13 = A1[13] ^ A2[13];
    int pop13 = __builtin_popcountll(xE13);
    int pop14 = __builtin_popcountll(xE14);
    int popU  = __builtin_popcountll(xE13 | xE14);

    u64 dAE15 = E1[15] - A1[15];
    u64 dA15 = A2[15] - A1[15];
    u64 orig_dT2_16 = (S0(A2[15]) + Maj(A2[15], A2[14], A2[13])) - (S0(A1[15]) + Maj(A1[15], A1[14], A1[13]));
    u64 orig_dCh16  = Ch(E2[15], E2[14], E2[13]) - Ch(E1[15], E1[14], E1[13]);
    int a15_free = 0;
    for (int b = 0; b < 64; b++) {
        u64 a = A1[15] ^ (1ULL << b), ap = a + dA15, e = a + dAE15;
        if ((ap ^ a) != (1ULL << 59)) continue;
        if (((e ^ E1[15]) & ((1ULL<<18)|(1ULL<<41)|(1ULL<<45))) != 0) continue;
        u64 dt2 = (S0(ap) + Maj(ap, A2[14], A2[13])) - (S0(a) + Maj(a, A1[14], A1[13]));
        u64 dch = Ch(e, E2[14], E2[13]) - Ch(e, E1[14], E1[13]);
        if (dt2 == orig_dT2_16 && dch == orig_dCh16) a15_free++;
    }
    u64 dE14 = E2[14] - E1[14], dW14 = W2_sol[14] - W1_sol[14];
    u64 ods1 = s1(W2_sol[14])-s1(W1_sol[14]), ods0 = s0(W2_sol[14])-s0(W1_sol[14]);
    u64 odt15 = (S1(E2[14])-S1(E1[14])) + (Ch(E2[14],E2[13],E2[12])-Ch(E1[14],E1[13],E1[12]));
    u64 odm15 = Maj(A2[14],A2[13],A2[12]) - Maj(A1[14],A1[13],A1[12]);
    int a14_free = 0;
    for (int b = 0; b < 64; b++) {
        u64 delta = (A1[14] ^ (1ULL << b)) - A1[14];
        u64 a14 = A1[14] + delta, e14 = E1[14] + delta, e14p = e14 + dE14;
        u64 w14 = W1_sol[14] + delta, w14p = w14 + dW14;
        if ((s1(w14p)-s1(w14)) != ods1) continue;
        if ((s0(w14p)-s0(w14)) != ods0) continue;
        if (((S1(e14p)-S1(e14)) + (Ch(e14p,E2[13],E2[12])-Ch(e14,E1[13],E1[12]))) != odt15) continue;
        if ((Maj(a14,A2[13],A2[12]) - Maj(a14,A1[13],A1[12])) != odm15) continue;
        if ((e14 ^ e14p) != xE14) continue;
        if (((e14 ^ E1[14]) & ((1ULL<<18)|(1ULL<<41)|(1ULL<<45))) != 0) continue;
        a14_free++;
    }
    u64 e16 = REF_E1[16], e16p = REF_E2[16];
    u64 dt1_17 = (E2[13]-E1[13]) + (S1(e16p)-S1(e16)) + (Ch(e16p, E1[15], E2[14]) - Ch(e16, E1[15], E1[14]));

    printf("[%s] VERIFY=%s | |xE13|=%d |xE14|=%d |xE13|xE14|=%d |xA13|=%d | A15_free=%d A14_free=%d (est 2^%d) | s17_compat=%s\n",
           tag, ok ? "PASS" : "FAIL", pop13, pop14, popU, __builtin_popcountll(xA13),
           a15_free, a14_free, a15_free + a14_free, (dt1_17 == 0) ? "EXACT_0" : "DIFF");
    if (ok && popU <= 35) {
        printf("  W1[9..15] =");
        for (int t = 9; t <= 15; t++) printf(" %016llx", (unsigned long long)W1_sol[t]);
        printf("\n  E1[13]=%016llx E2[13]=%016llx xE13=%016llx\n",
               (unsigned long long)E1[13], (unsigned long long)E2[13], (unsigned long long)xE13);
    }
    fflush(stdout);
}

int main(int argc, char **argv) {
    int fix9 = 0, mode_verify = 0, start_B = 34, min_B = 5, seed = 1, enum_k = 0, free_e14_carry = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--verify")) mode_verify = 1;
        else if (!strcmp(argv[i], "--fix9")) fix9 = 1;
        else if (!strncmp(argv[i], "--startB=", 9)) start_B = atoi(argv[i]+9);
        else if (!strncmp(argv[i], "--minB=", 7)) min_B = atoi(argv[i]+7);
        else if (!strncmp(argv[i], "--seed=", 7)) seed = atoi(argv[i]+7);
        else if (!strncmp(argv[i], "--enum=", 7)) enum_k = atoi(argv[i]+7);
        else if (!strcmp(argv[i], "--free14")) free_e14_carry = 1;
    }

    CaDiCaL::Solver S;
    S.set("seed", seed);
    S.set("quiet", 1);
    Builder B(S);

    std::array<int,64> A1[17], A2[17], E1[17], E2[17], W1[16], W2[16];
    for (int t = 5; t <= 8; t++) {
        A1[t] = B.const_word(REF_A1[t]); A2[t] = B.const_word(REF_A2[t]);
        E1[t] = B.const_word(REF_E1[t]); E2[t] = B.const_word(REF_E2[t]);
    }
    if (fix9) {
        A1[9] = B.const_word(REF_A1[9]); A2[9] = B.const_word(REF_A2[9]);
        E1[9] = B.const_word(REF_E1[9]); E2[9] = B.const_word(REF_E2[9]);
        W1[9] = B.const_word(REF_W1[9]); W2[9] = B.const_word(REF_W2[9]);
    }

    int t0 = fix9 ? 10 : 9;
    for (int t = t0; t <= 15; t++) {
        W1[t] = B.new_word();
        if (t == 10 || t == 11 || t == 12 || t == 15) W2[t] = W1[t];
        else W2[t] = B.add2(W1[t], B.const_word(REF_W2[t] - REF_W1[t]));

        auto T1_1 = B.add_many({E1[t-4], B.w_S1(E1[t-1]), B.w_Ch(E1[t-1], E1[t-2], E1[t-3]), B.const_word(K512[t]), W1[t]});
        auto T2_1 = B.add2(B.w_S0(A1[t-1]), B.w_Maj(A1[t-1], A1[t-2], A1[t-3]));
        A1[t] = B.add2(T1_1, T2_1);
        E1[t] = B.add2(A1[t-4], T1_1);

        auto T1_2 = B.add_many({E2[t-4], B.w_S1(E2[t-1]), B.w_Ch(E2[t-1], E2[t-2], E2[t-3]), B.const_word(K512[t]), W2[t]});
        auto T2_2 = B.add2(B.w_S0(A2[t-1]), B.w_Maj(A2[t-1], A2[t-2], A2[t-3]));
        A2[t] = B.add2(T1_2, T2_2);
        E2[t] = B.add2(A2[t-4], T1_2);
    }

    if (!fix9) {
        auto s0_9_1 = B.w_s0(W1[9]), s0_9_2 = B.w_s0(W2[9]);
        B.eq_word(s0_9_2, B.add2(s0_9_1, B.const_word(s0(REF_W2[9]) - s0(REF_W1[9]))));
    }
    {
        auto s0_13_1 = B.w_s0(W1[13]), s0_13_2 = B.w_s0(W2[13]);
        B.eq_word(s0_13_2, B.add2(s0_13_1, B.const_word(s0(REF_W2[13]) - s0(REF_W1[13]))));
    }
    {
        auto s0_14_1 = B.w_s0(W1[14]), s0_14_2 = B.w_s0(W2[14]);
        B.eq_word(s0_14_2, B.add2(s0_14_1, B.const_word(s0(REF_W2[14]) - s0(REF_W1[14]))));
        auto s1_14_1 = B.w_s1(W1[14]), s1_14_2 = B.w_s1(W2[14]);
        B.eq_word(s1_14_2, B.add2(s1_14_1, B.const_word(s1(REF_W2[14]) - s1(REF_W1[14]))));
    }

    B.eq_word(A2[13], B.add2(A1[13], B.const_word(REF_A2[13] - REF_A1[13])));
    B.add({-A1[13][43]}); B.add({-A1[13][51]}); B.add({-A1[13][57]});
    B.eq_word(E2[13], B.add2(E1[13], B.const_word(REF_E2[13] - REF_E1[13])));
    B.eq_word(A2[14], A1[14]);
    B.eq_word(E2[14], B.add2(E1[14], B.const_word(REF_E2[14] - REF_E1[14])));
    if (!free_e14_carry) {
        u64 x14 = REF_E1[14] ^ REF_E2[14];
        for (int b = 0; b < 64; b++) {
            if ((x14 >> b) & 1) {
                B.add({((REF_E1[14] >> b) & 1) ? E1[14][b] : -E1[14][b]});
                B.add({((REF_E2[14] >> b) & 1) ? E2[14][b] : -E2[14][b]});
            } else {
                B.eq(E1[14][b], E2[14][b]);
            }
        }
    }
    B.eq_word(A2[15], B.add2(A1[15], B.const_word(REF_A2[15] - REF_A1[15])));
    B.add({A1[15][59]});
    B.eq_word(E2[15], E1[15]);
    B.add({E1[14][18]}); B.add({-E1[15][18]});
    B.eq(E1[15][41], E1[14][41]);
    B.eq(E1[15][45], E1[14][45]);

    auto ch1_16 = B.w_Ch(E1[15], E1[14], E1[13]);
    auto ch2_16 = B.w_Ch(E2[15], E2[14], E2[13]);
    auto m1_16  = B.w_Maj(A1[15], A1[14], A1[13]);
    auto m2_16  = B.w_Maj(A2[15], A2[14], A2[13]);
    auto lhs_a16 = B.add_many({E1[12], ch1_16, B.w_S0(A1[15]), m1_16});
    auto rhs_a16 = B.add_many({E2[12], ch2_16, B.w_S0(A2[15]), m2_16});
    B.eq_word(lhs_a16, rhs_a16);
    auto lhs_e16 = B.add_many({A1[12], E1[12], ch1_16, B.const_word(REF_E2[16] - REF_E1[16])});
    auto rhs_e16 = B.add_many({A2[12], E2[12], ch2_16});
    B.eq_word(lhs_e16, rhs_e16);

    if (mode_verify) {
        for (int t = t0; t <= 15; t++) B.fix_word(W1[t], REF_W1[t]);
        printf("CNF built: vars=%d clauses=%d (verify mode, fix9=%d)\n", B.nv, B.ncl, fix9);
        auto t_s = std::chrono::steady_clock::now();
        int res = S.solve();
        double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_s).count();
        printf("Verify solve: res=%d in %.4fs\n", res, dt);
        if (res == 10) {
            u64 w1[16], w2[16];
            for (int t = 0; t < 16; t++) { w1[t] = (t < t0) ? REF_W1[t] : B.val_word(W1[t]); w2[t] = (t < t0) ? REF_W2[t] : B.val_word(W2[t]); }
            verify_and_report(w1, w2, "VERIFY_T14");
        }
        return 0;
    }

    std::vector<int> x13_bits;
    for (int b = 18; b < 64; b++) x13_bits.push_back(B.xor2(E1[13][b], E2[13][b]));
    auto cnt = B.build_seq_counter(x13_bits, start_B + 1);
    printf("CNF built (fix9=%d, free14=%d, seed=%d): vars=%d clauses=%d\n", fix9, free_e14_carry, seed, B.nv, B.ncl);
    fflush(stdout);

    if (enum_k > 0) {
        if (start_B < (int)cnt.size()) B.add({-cnt[start_B]});
        auto t_all = std::chrono::steady_clock::now();
        for (int k = 0; k < enum_k; k++) {
            auto t_s = std::chrono::steady_clock::now();
            int res = S.solve();
            double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_s).count();
            if (res != 10) { printf("Enum stopped at k=%d (res=%d, %.2fs)\n", k, res, dt); break; }
            u64 w1[16], w2[16];
            for (int t = 0; t < 16; t++) { w1[t] = (t < t0) ? REF_W1[t] : B.val_word(W1[t]); w2[t] = (t < t0) ? REF_W2[t] : B.val_word(W2[t]); }
            if (k < 5 || (k+1) % 10 == 0 || k+1 == enum_k) {
                char tag[64]; snprintf(tag, sizeof(tag), "ENUM_%d %.3fs", k+1, dt);
                verify_and_report(w1, w2, tag);
            }
            for (int b = 0; b < 64; b++) S.add(((w1[12] >> b) & 1) ? -W1[12][b] : W1[12][b]);
            S.add(0);
        }
        double dt_all = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_all).count();
        printf("Total enum time: %.2fs (%.2f cores/sec)\n", dt_all, enum_k / dt_all);
        return 0;
    }

    int curB = start_B;
    while (curB >= min_B) {
        printf("--- Trying |xE13| <= %d ---\n", curB); fflush(stdout);
        S.assume(-cnt[curB]);
        auto t_s = std::chrono::steady_clock::now();
        int res = S.solve();
        double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_s).count();
        if (res == 10) {
            u64 w1[16], w2[16];
            for (int t = 0; t < 16; t++) { w1[t] = (t < t0) ? REF_W1[t] : B.val_word(W1[t]); w2[t] = (t < t0) ? REF_W2[t] : B.val_word(W2[t]); }
            u64 e13_1 = B.val_word(E1[13]), e13_2 = B.val_word(E2[13]);
            int got = __builtin_popcountll(e13_1 ^ e13_2);
            char tag[64]; snprintf(tag, sizeof(tag), "B<=%d got=%d %.2fs", curB, got, dt);
            verify_and_report(w1, w2, tag);
            B.add({-cnt[got]});
            curB = got - 1;
        } else {
            printf("  |xE13| <= %d -> res=%d in %.2fs\n", curB, res, dt);
            break;
        }
    }
    return 0;
}
