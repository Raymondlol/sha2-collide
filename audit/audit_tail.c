// Independent CPU audit of the sha512_36_t14 tail model (CPU only, no GPU).
// Idea: take the characteristic's per-step modular state differences from the Table 14 pair,
// generate REAL (lower, upper) message pairs exactly as the GPU search does, and for each pair
// test every step t generically: "given char diff before step t and message-1's real state,
// does step t produce the char diff after step t?"  c[t] = that boolean.
// Then measure P(c23) (the GPU's "W21" counter = dW23 check), P(c25|c23) (GPU's "W23" counter),
// conditional rates P(c23 | c_t) for t=14..22, P(c23 | c19..c22), and whether c24..c35 hold
// once c23&c25 hold (structural: collision guaranteed after the last counter?).
// Build: gcc -O3 -march=native -fopenmp audit_tail.c -o audit_tail -lm
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

static inline uint64_t rotr(uint64_t x, int n) { return (x >> n) | (x << (64 - n)); }
static inline uint64_t S0(uint64_t x) { return rotr(x,28) ^ rotr(x,34) ^ rotr(x,39); }
static inline uint64_t S1(uint64_t x) { return rotr(x,14) ^ rotr(x,18) ^ rotr(x,41); }
static inline uint64_t s0(uint64_t x) { return rotr(x,1) ^ rotr(x,8) ^ (x >> 7); }
static inline uint64_t s1(uint64_t x) { return rotr(x,19) ^ rotr(x,61) ^ (x >> 6); }
static inline uint64_t Ch(uint64_t e, uint64_t f, uint64_t g) { return (e & f) ^ ((~e) & g); }
static inline uint64_t Maj(uint64_t a, uint64_t b, uint64_t c) { return (a & b) ^ (a & c) ^ (b & c); }

static const uint64_t K[36] = {
    0x428a2f98d728ae22ULL, 0x7137449123ef65cdULL, 0xb5c0fbcfec4d3b2fULL, 0xe9b5dba58189dbbcULL,
    0x3956c25bf348b538ULL, 0x59f111f1b605d019ULL, 0x923f82a4af194f9bULL, 0xab1c5ed5da6d8118ULL,
    0xd807aa98a3030242ULL, 0x12835b0145706fbeULL, 0x243185be4ee4b28cULL, 0x550c7dc3d5ffb4e2ULL,
    0x72be5d74f27b896fULL, 0x80deb1fe3b1696b1ULL, 0x9bdc06a725c71235ULL, 0xc19bf174cf692694ULL,
    0xe49b69c19ef14ad2ULL, 0xefbe4786384f25e3ULL, 0x0fc19dc68b8cd5b5ULL, 0x240ca1cc77ac9c65ULL,
    0x2de92c6f592b0275ULL, 0x4a7484aa6ea6e483ULL, 0x5cb0a9dcbd41fbd4ULL, 0x76f988da831153b5ULL,
    0x983e5152ee66dfabULL, 0xa831c66d2db43210ULL, 0xb00327c898fb213fULL, 0xbf597fc7beef0ee4ULL,
    0xc6e00bf33da88fc2ULL, 0xd5a79147930aa725ULL, 0x06ca6351e003826fULL, 0x142929670a0e6e70ULL,
    0x27b70a8546d22ffcULL, 0x2e1b21385c26c926ULL, 0x4d2c6dfc5ac42aedULL, 0x53380d139d95b3dfULL,
};
static const uint64_t IV[8] = {
    0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL, 0x3c6ef372fe94f82bULL, 0xa54ff53a5f1d36f1ULL,
    0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL, 0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL,
};

// ---- Table 14 SFS pair (ePrint 2026/1120), as used by the project's sentinel ----
static const uint64_t T14_CV[8] = {
    0xf128698fe95be8acULL, 0x264034c8f44fe6edULL, 0x7b976a39ced01663ULL, 0x7b671a27b6f1da99ULL,
    0x198021f679057203ULL, 0x600969ddc8027d0cULL, 0x4af8d8a11db1692bULL, 0x98ba089bc0dee788ULL };
static const uint64_t T14_M1[16] = {
    0x230a4930d38229e3ULL, 0x2b3c63372261608bULL, 0xf833951c32f2b6dcULL, 0x6f4f49447a93b0ccULL,
    0xf841e91df9c8d69eULL, 0xa4710d44716286a3ULL, 0xe3310abbec7a9816ULL, 0x66b0241d0bfd564cULL,
    0x5395fbd4a70df9f3ULL, 0xb0f8774e493b43bcULL, 0x13287fa0a341d877ULL, 0xf0d19e5ae32f6abfULL,
    0xf26d943cf8041873ULL, 0x06f2b17bb37eb2b7ULL, 0x2cccf19ec39c66ddULL, 0xf284ae681c043be9ULL };
static const uint64_t T14_M2[16] = {
    0x230a4930d38229e3ULL, 0x2b3c63372261608bULL, 0xf833951c32f2b6dcULL, 0x6f4f49447a93b0ccULL,
    0xf841e91df9c8d69eULL, 0xac710d44716286a3ULL, 0xe7390abbec7a9816ULL, 0x6eb0241d0bfd564cULL,
    0x73c5db55a70df9f3ULL, 0xf0d8764e493b43bcULL, 0x13287fa0a341d877ULL, 0xf0d19e5ae32f6abfULL,
    0xf26d943cf8041873ULL, 0x02fab17bb37eb2b7ULL, 0x24ccf19ec39c66ddULL, 0xf284ae681c043be9ULL };

// ---- Generator constants (copied from src/gpu/sha512_36_t14_gpu.cu) ----
static const uint64_t C_T1_4=0xf3aea7c98ce111b8ULL, C_T1_5=0xeed3968f2fef0f6fULL, C_T1_6=0x0582af92b9b1c77dULL;
static const uint64_t C_T1_7=0x1a7b8214f1370046ULL, C_T1_8=0x9a4c318c1745a72fULL, C_ORIG_dCh6=0xf000000000000000ULL;
static const uint64_t C_dW5=0x0800000000000000ULL, C_dW6=0x0408000000000000ULL, C_dW7=0x0800000000000000ULL;
static const uint64_t C_dW8=0x202fdf8100000000ULL, C_ds0_W6=0xfdf7f80000000000ULL, C_ds0_W7=0xfbf8000000000000ULL;
static const uint64_t C_ds0_W8=0xf0385fe0ff000000ULL;
static const uint64_t C_A1_1=0xd0d35333c5347d6cULL, C_A1_2=0x0878e2054cfcadc7ULL, C_A1_3=0x2b001dd1868a62e9ULL;
static const uint64_t C_A1_4=0xa6089dc87a3559fcULL, C_E1_4=0x8761759d36cebdecULL, C_E1_5=0xbfa6e9c2f5238cdbULL;
static const uint64_t C_E2_5=0xc7a6e9c2f5238cdbULL, C_E1_6=0x0dfb919806ae7544ULL, C_E1_7=0x457b9fe677c1632fULL;
static const uint64_t C_A1_11=0xe3698ee4fe1319afULL, C_A1_12=0x863aebeaa5768210ULL, C_A1_13=0x44f0c7b6bd1e6710ULL;
static const uint64_t C_A2_13=0x46f8cfb6bd1e6710ULL, C_A1_14=0x123803f0e014e59aULL, C_A1_15=0x7b5a8e1c3bb8af4bULL;
static const uint64_t C_E1_11=0xa15c89731a29737dULL, C_E1_12=0x3319e0efefabffffULL, C_E1_13=0x6fffff2fd821023fULL;
static const uint64_t C_E2_13=0x70000132579d023fULL, C_E1_14=0xf3f8dd0307d580adULL, C_E2_14=0xf3f8dd6125d484a9ULL;
static const uint64_t C_E1_15=0x80001002018a91f5ULL;
static const uint64_t C_W1_9=0xb0f8774e493b43bcULL, C_W1_10=0x13287fa0a341d877ULL, C_W1_11=0xf0d19e5ae32f6abfULL;
static const uint64_t C_W1_12=0xf26d943cf8041873ULL, C_W1_13=0x06f2b17bb37eb2b7ULL, C_W1_14=0x2cccf19ec39c66ddULL;
static const uint64_t C_W1_15_REF=0xf284ae681c043be9ULL;

static inline void stepf(uint64_t *s, uint64_t w, int t) {
    uint64_t t1 = s[7] + S1(s[4]) + Ch(s[4], s[5], s[6]) + K[t] + w;
    uint64_t t2 = S0(s[0]) + Maj(s[0], s[1], s[2]);
    s[7] = s[6]; s[6] = s[5]; s[5] = s[4]; s[4] = s[3] + t1;
    s[3] = s[2]; s[2] = s[1]; s[1] = s[0]; s[0] = t1 + t2;
}
static inline void expand(const uint64_t *M, uint64_t *W) {
    for (int i = 0; i < 16; i++) W[i] = M[i];
    for (int i = 16; i < 36; i++) W[i] = s1(W[i-2]) + W[i-7] + s0(W[i-15]) + W[i-16];
}
static inline uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static uint64_t DS[37][8];   // char modular state diff: DS[t] = before step t (DS[0] = 0 at CV)
static uint64_t DM[16];      // char message diffs
static uint64_t DW[36];      // char expanded diffs

// ---- generic per-step check from a known message-1 state ----
// s: message-1 state before step t (updated in place). Returns 1 if step t maps DS[t] -> DS[t+1].
static inline int check_step(uint64_t *s, uint64_t w1, uint64_t w2, int t) {
    uint64_t p[8];
    for (int r = 0; r < 8; r++) p[r] = s[r] + DS[t][r];
    stepf(s, w1, t);
    stepf(p, w2, t);
    return ((p[0] - s[0]) == DS[t+1][0]) && ((p[4] - s[4]) == DS[t+1][4]);
}

static inline uint64_t e4_from_idx(uint32_t idx) {
    uint64_t d = ((uint64_t)(idx & 0xFULL)) |
        ((uint64_t)((idx >> 4) & 0x1FFFULL) << 5) | ((uint64_t)((idx >> 17) & 0x1FULL) << 20) |
        ((uint64_t)((idx >> 22) & 0x7ULL) << 27) | ((uint64_t)((idx >> 25) & 0x7ULL) << 34) |
        ((uint64_t)((idx >> 28) & 1ULL) << 63);
    return C_E1_4 ^ d;
}

typedef struct { uint64_t seed, cv[8], C3, C2, C1, C0, W0, W1, W2, A2A1_and, A2A1_xor, A1Am1_and, A1Am1_xor, Am1Am2_and, Am1Am2_xor; } BCV;

static void compress36(const uint64_t iv[8], const uint64_t M[16], uint64_t out[8]) {
    uint64_t W[36], s[8];
    expand(M, W);
    for (int i = 0; i < 8; i++) s[i] = iv[i];
    for (int t = 0; t < 36; t++) stepf(s, W[t], t);
    for (int i = 0; i < 8; i++) out[i] = iv[i] + s[i];
}

static int is_golden(uint64_t seed, uint64_t cv[8]) {
    uint64_t M0[16] = {0};
    M0[0] = seed; M0[1] = 0x5348413531323336ULL; M0[2] = 0x434f4c4c49444521ULL;
    compress36(IV, M0, cv);
    const uint64_t A0_base = C_E1_4 - C_T1_4;
    const uint64_t E3_sub  = C_A1_3 - S0(C_A1_2) - Maj(C_A1_2, C_A1_1, A0_base);
    const uint64_t W7_base = C_T1_7 - S1(C_E1_6) - Ch(C_E1_6, C_E1_5, C_E1_4) - K[7];
    const uint64_t W6_base = C_T1_6 - S1(C_E1_5) - K[6];
    const uint64_t E2_sub  = C_A1_2 - S0(C_A1_1);
    uint64_t Am1 = cv[0], Am2 = cv[1];
    uint64_t E3_n = Am1 + E3_sub;
    if ((Ch(C_E2_5, C_E1_4, E3_n) - Ch(C_E1_5, C_E1_4, E3_n)) != C_ORIG_dCh6) return 0;
    uint64_t W7_n = W7_base - E3_n;
    if ((s0(W7_n + C_dW7) - s0(W7_n)) != C_ds0_W7) return 0;
    uint64_t E2_n = Am2 + E2_sub - Maj(C_A1_1, A0_base, Am1);
    uint64_t W6_n = W6_base - E2_n - Ch(C_E1_5, C_E1_4, E3_n);
    if ((s0(W6_n + C_dW6) - s0(W6_n)) != C_ds0_W6) return 0;
    return 1;
}

static BCV make_bcv(uint64_t seed, const uint64_t cv[8]) {
    BCV B; B.seed = seed; for (int i = 0; i < 8; i++) B.cv[i] = cv[i];
    uint64_t Am1=cv[0], Am2=cv[1], Am3=cv[2], Am4=cv[3], Em1=cv[4], Em2=cv[5], Em3=cv[6], Em4=cv[7];
    B.C3 = Am1 + C_A1_3 - S0(C_A1_2);
    B.C2 = Am2 + C_A1_2 - S0(C_A1_1);
    B.C1 = Am3 + C_A1_1;
    B.C0 = Am4 - S0(Am1) - Maj(Am1, Am2, Am3);
    B.W0 = 0ULL - S0(Am1) - Maj(Am1, Am2, Am3) - Em4 - S1(Em1) - Ch(Em1, Em2, Em3) - K[0];
    B.W1 = C_A1_1 - Em3 - K[1];
    B.W2 = C_A1_2 - S0(C_A1_1) - Em2 - K[2];
    B.A2A1_and = C_A1_2 & C_A1_1; B.A2A1_xor = C_A1_2 ^ C_A1_1;
    B.A1Am1_and = C_A1_1 & Am1;   B.A1Am1_xor = C_A1_1 ^ Am1;
    B.Am1Am2_and = Am1 & Am2;     B.Am1Am2_xor = Am1 ^ Am2;
    return B;
}

// Lower words W0..W8 for idx; returns GPU-filter pass flag; also reports W5 identity vs true condition.
static uint64_t TGT_ds0_W5;
static inline int make_lower(const BCV *B, uint32_t idx, uint64_t w[9], int *id_pass, int *true_pass) {
    uint64_t e4 = e4_from_idx(idx);
    uint64_t a0 = e4 - C_T1_4;
    uint64_t e3 = B->C3 - (B->A2A1_and | (B->A2A1_xor & a0));
    uint64_t e2 = B->C2 - (B->A1Am1_and | (B->A1Am1_xor & a0));
    uint64_t e1 = B->C1 - S0(a0) - (B->Am1Am2_and | (B->Am1Am2_xor & a0));
    uint64_t e0 = B->C0 + a0;
    w[0] = B->W0 + a0;
    w[1] = B->W1 - S0(a0) - (B->Am1Am2_and | (B->Am1Am2_xor & a0)) - S1(e0) - Ch(e0, B->cv[4], B->cv[5]);
    w[2] = B->W2 - (B->A1Am1_and | (B->A1Am1_xor & a0)) - S1(e1) - Ch(e1, e0, B->cv[4]);
    w[3] = C_A1_3 - S0(C_A1_2) - Maj(C_A1_2, C_A1_1, a0) - B->cv[4] - S1(e2) - Ch(e2, e1, e0) - K[3];
    w[4] = C_A1_4 - S0(C_A1_3) - Maj(C_A1_3, C_A1_2, C_A1_1) - e0 - S1(e3) - Ch(e3, e2, e1) - K[4];
    w[5] = C_T1_5 - e1 - S1(e4) - Ch(e4, e3, e2) - K[5];
    w[6] = C_T1_6 - e2 - S1(C_E1_5) - Ch(C_E1_5, e4, e3) - K[6];
    w[7] = C_T1_7 - e3 - S1(C_E1_6) - Ch(C_E1_6, C_E1_5, e4) - K[7];
    w[8] = C_T1_8 - e4 - S1(C_E1_7) - Ch(C_E1_7, C_E1_6, C_E1_5) - K[8];
    uint64_t w5 = w[5];
    int idp = ((w5 & 0x0800000000000004ULL) == 0ULL) &&
              ((((w5 >> 52) ^ (w5 >> 58)) & 1ULL) == 0ULL) && ((((w5 >> 53) ^ (w5 >> 60)) & 1ULL) == 1ULL);
    *id_pass = idp;
    *true_pass = ((s0(w5 + C_dW5) - s0(w5)) == TGT_ds0_W5);
    if (!idp) return 0;
    if ((s0(w[6] + C_dW6) - s0(w[6])) != C_ds0_W6) return 0;
    if ((s0(w[8] + C_dW8) - s0(w[8])) != C_ds0_W8) return 0;
    if ((Ch(C_E2_5, e4, e3) - Ch(C_E1_5, e4, e3)) != C_ORIG_dCh6) return 0;
    if ((s0(w[7] + C_dW7) - s0(w[7])) != C_ds0_W7) return 0;
    return 1;
}

// ---- Upper set (port of build_upper_2_25 DP) ----
static uint64_t h_a14[16], h_e14[16], h_w14[16], h_dae15[16], h_w15_sub[16];
static uint64_t *A15[16];
static uint32_t nA15[16];

static void build_upper(void) {
    const int a14_bits[4] = {28, 40, 62, 63};
    uint64_t e15_fix_mask = (1ULL<<18) | (1ULL<<41) | (1ULL<<45);
    for (int b = 0; b < 64; b++) {
        int c0 = (int)((Ch(0ULL, C_E2_14, C_E2_13) >> b) & 1ULL) - (int)((Ch(0ULL, C_E1_14, C_E1_13) >> b) & 1ULL);
        int c1 = (int)((Ch(~0ULL, C_E2_14, C_E2_13) >> b) & 1ULL) - (int)((Ch(~0ULL, C_E1_14, C_E1_13) >> b) & 1ULL);
        if (c0 != c1) e15_fix_mask |= (1ULL << b);
    }
    const int tgt_p = (int)(((C_A1_15 >> 1) ^ (C_A1_15 >> 6)) & 1ULL) |
                      ((int)(((C_A1_15 >> 0) ^ (C_A1_15 >> 53)) & 1ULL) << 1) |
                      ((int)(((C_A1_15 >> 48) ^ (C_A1_15 >> 54)) & 1ULL) << 2);
    for (int m = 0; m < 16; m++) {
        uint64_t a14 = C_A1_14;
        for (int i = 0; i < 4; i++) if ((m >> i) & 1) a14 ^= (1ULL << a14_bits[i]);
        uint64_t delta = a14 - C_A1_14;
        uint64_t e14 = C_E1_14 + delta, w14 = C_W1_14 + delta;
        uint64_t t2_15 = S0(a14) + Maj(a14, C_A1_13, C_A1_12);
        h_a14[m] = a14; h_e14[m] = e14; h_w14[m] = w14;
        h_dae15[m] = C_A1_11 - t2_15;
        h_w15_sub[m] = t2_15 + C_E1_11 + S1(e14) + Ch(e14, C_E1_13, C_E1_12) + K[15];
    }
    #pragma omp parallel for schedule(static)
    for (int m = 0; m < 16; m++) {
        uint64_t a14 = h_a14[m], dae15 = h_dae15[m];
        uint64_t am = (1ULL << 59);
        for (int b = 0; b < 64; b++) {
            if (b == 59) continue;
            int m0 = (int)((Maj(0ULL, a14, C_A2_13) >> b) & 1ULL) - (int)((Maj(0ULL, a14, C_A1_13) >> b) & 1ULL);
            int m1 = (int)((Maj(~0ULL, a14, C_A2_13) >> b) & 1ULL) - (int)((Maj(~0ULL, a14, C_A1_13) >> b) & 1ULL);
            if (m0 != m1) am |= (1ULL << b);
        }
        static __thread uint8_t reach[65][16];
        memset(reach, 0, sizeof(reach));
        reach[64][tgt_p] = 1; reach[64][8 + tgt_p] = 1;
        for (int b = 63; b >= 0; b--) {
            int da_b = (int)((dae15 >> b) & 1ULL);
            int a_r = ((am >> b) & 1ULL) ? (int)((C_A1_15 >> b) & 1ULL) : -1;
            int e_r = ((e15_fix_mask >> b) & 1ULL) ? (int)((C_E1_15 >> b) & 1ULL) : -1;
            int pm = ((b==1||b==6)?1:0) | ((b==0||b==53)?2:0) | ((b==48||b==54)?4:0);
            for (int st = 0; st < 16; st++) {
                int c = st >> 3, p = st & 7;
                for (int ab = 0; ab <= 1; ab++) {
                    if (a_r != -1 && ab != a_r) continue;
                    int s = ab + da_b + c;
                    if (e_r != -1 && (s & 1) != e_r) continue;
                    int nst = ((s >> 1) << 3) | (p ^ (ab ? pm : 0));
                    if (reach[b + 1][nst]) { reach[b][st] = 1; break; }
                }
            }
        }
        uint64_t *sols = (uint64_t*)malloc(sizeof(uint64_t) * (1U << 22));
        uint32_t ns = 0;
        struct Frame { int b, st; uint64_t val; } stk[256]; int sp = 0;
        stk[sp++] = (struct Frame){0, 0, 0ULL};
        while (sp > 0) {
            struct Frame fr = stk[--sp];
            if (fr.b == 64) { if (ns < (1U << 22)) sols[ns++] = fr.val; continue; }
            int b = fr.b, c = fr.st >> 3, p = fr.st & 7;
            int da_b = (int)((dae15 >> b) & 1ULL);
            int a_r = ((am >> b) & 1ULL) ? (int)((C_A1_15 >> b) & 1ULL) : -1;
            int e_r = ((e15_fix_mask >> b) & 1ULL) ? (int)((C_E1_15 >> b) & 1ULL) : -1;
            int pm = ((b==1||b==6)?1:0) | ((b==0||b==53)?2:0) | ((b==48||b==54)?4:0);
            for (int ab = 0; ab <= 1; ab++) {
                if (a_r != -1 && ab != a_r) continue;
                int s = ab + da_b + c;
                if (e_r != -1 && (s & 1) != e_r) continue;
                int nst = ((s >> 1) << 3) | (p ^ (ab ? pm : 0));
                if (reach[b + 1][nst]) stk[sp++] = (struct Frame){b + 1, nst, fr.val | ((uint64_t)ab << b)};
            }
        }
        A15[m] = sols; nA15[m] = ns;
    }
}

int main(int argc, char **argv) {
    double n_pairs_target = argc > 1 ? atof(argv[1]) : 1e9;
    int uppers_per_lower = argc > 2 ? atoi(argv[2]) : 4096;
    int n_golden = argc > 3 ? atoi(argv[3]) : 2048;

    // 1) characteristic from Table 14 pair
    {
        uint64_t W1[36], W2[36], a[8], b[8];
        expand(T14_M1, W1); expand(T14_M2, W2);
        for (int i = 0; i < 16; i++) DM[i] = T14_M2[i] - T14_M1[i];
        for (int i = 0; i < 36; i++) DW[i] = W2[i] - W1[i];
        for (int i = 0; i < 8; i++) { a[i] = T14_CV[i]; b[i] = T14_CV[i]; DS[0][i] = 0; }
        for (int t = 0; t < 36; t++) { stepf(a, W1[t], t); stepf(b, W2[t], t); for (int i = 0; i < 8; i++) DS[t+1][i] = b[i] - a[i]; }
        int coll = 1; for (int i = 0; i < 8; i++) if (DS[36][i] != 0) coll = 0;
        printf("[char] Table 14 pair: final state diff zero = %s\n", coll ? "YES" : "NO");
        printf("[char] nonzero dW_t: ");
        for (int t = 0; t < 36; t++) if (DW[t]) printf("%d ", t);
        printf("\n[char] dW23=%016llx dW25=%016llx\n", (unsigned long long)DW[23], (unsigned long long)DW[25]);
        TGT_ds0_W5 = s0(T14_M2[5]) - s0(T14_M1[5]);
        printf("[char] ds0(W5) target=%016llx | dM: ", (unsigned long long)TGT_ds0_W5);
        for (int i = 0; i < 16; i++) if (DM[i]) printf("W%d=%016llx ", i, (unsigned long long)DM[i]);
        printf("\n");
        int ok = DM[5]==C_dW5 && DM[6]==C_dW6 && DM[7]==C_dW7 && DM[8]==C_dW8;
        printf("[char] message diffs match GPU constants: %s\n", ok ? "YES" : "NO");
        if (!coll || !ok) return 1;
    }

    // 2) upper set
    double t0 = omp_get_wtime();
    build_upper();
    uint64_t totU = 0; for (int m = 0; m < 16; m++) totU += nA15[m];
    printf("[upper] built %llu a15 values in %.1fs (per-m:", (unsigned long long)totU, omp_get_wtime() - t0);
    for (int m = 0; m < 16; m++) printf(" %u", nA15[m]); printf(")\n");
    int sent_ok = 0;
    for (uint32_t i = 0; i < nA15[0]; i++) if (A15[0][i] == C_A1_15 && (A15[0][i] - h_w15_sub[0]) == C_W1_15_REF) { sent_ok = 1; break; }
    printf("[upper] Table 14 sentinel (m=0, a15=A1_15 -> W15_ref) present: %s\n", sent_ok ? "YES" : "NO");
    // per-m W14 structural conditions needed later: dW16 = s1-diff(W14)+dW9 == 0, dW29 needs ds0(W14) == ds0 target
    uint64_t tgt_ds1_w14 = s1(T14_M2[14]) - s1(T14_M1[14]);
    uint64_t tgt_ds0_w14 = s0(T14_M2[14]) - s0(T14_M1[14]);
    int w14_bad = 0;
    for (int m = 0; m < 16; m++) {
        uint64_t w = h_w14[m];
        int okA = (s1(w + DM[14]) - s1(w)) == tgt_ds1_w14;
        int okB = (s0(w + DM[14]) - s0(w)) == tgt_ds0_w14;
        if (!okA || !okB) { w14_bad++; printf("[upper] m=%d W14 ds1 %s ds0 %s\n", m, okA ? "ok" : "BAD", okB ? "ok" : "BAD"); }
    }
    printf("[upper] W14 sigma-difference conditions: %d/16 bad\n", w14_bad);

    // 3) golden M0s on CPU
    t0 = omp_get_wtime();
    BCV *G = (BCV*)malloc(sizeof(BCV) * n_golden);
    int ng = 0;
    #pragma omp parallel
    {
        uint64_t seed = 0x7a00000000000000ULL + (uint64_t)omp_get_thread_num() * 0x100000000ULL;
        uint64_t cv[8];
        while (1) {
            int done;
            #pragma omp atomic read
            done = ng;
            if (done >= n_golden) break;
            for (int k = 0; k < 4096; k++, seed++) {
                if (is_golden(seed, cv)) {
                    int slot;
                    #pragma omp atomic capture
                    slot = ng++;
                    if (slot < n_golden) G[slot] = make_bcv(seed, cv);
                }
            }
        }
    }
    ng = n_golden;
    printf("[golden] %d golden M0 in %.1fs\n", ng, omp_get_wtime() - t0);

    // 4) main sampling
    uint64_t n_lower_tried = 0, n_lower_ok = 0, n_id = 0, n_true = 0, n_id_not_true = 0, n_true_not_id = 0;
    uint64_t low_fail[16] = {0};
    uint64_t n_pairs = 0, pass[36] = {0}, and23[36] = {0}, nJ = 0, nJ23 = 0, n23 = 0, n2325 = 0, ntail = 0;
    uint64_t tailfail[36] = {0}, pm_pairs[16] = {0}, pm_23[16] = {0}, pm_2325[16] = {0}, pm_tail[16] = {0};
    uint64_t nJ2 = 0, nJ223 = 0, nJa = 0, nJa23 = 0, nJb = 0, nJb23 = 0;
    uint64_t target_pairs = (uint64_t)n_pairs_target;
    t0 = omp_get_wtime();
    #pragma omp parallel reduction(+:n_lower_tried,n_lower_ok,n_id,n_true,n_id_not_true,n_true_not_id,n_pairs,nJ,nJ23,n23,n2325,ntail,nJ2,nJ223,nJa,nJa23,nJb,nJb23)
    {
        uint64_t rng = 0xc0ffee1234567890ULL ^ ((uint64_t)(omp_get_thread_num() + 1) * 0x9e3779b97f4a7c15ULL);
        uint64_t l_pass[36] = {0}, l_and23[36] = {0}, l_tailfail[36] = {0}, l_low_fail[16] = {0};
        uint64_t l_pm_pairs[16] = {0}, l_pm_23[16] = {0}, l_pm_2325[16] = {0}, l_pm_tail[16] = {0};
        int nth = omp_get_num_threads();
        uint64_t my_target = target_pairs / nth;
        uint64_t mine = 0;
        while (mine < my_target) {
            const BCV *B = &G[splitmix64(&rng) % ng];
            uint64_t lw[9]; int idp, trp;
            uint32_t idx = (uint32_t)(splitmix64(&rng) & ((1U << 29) - 1));
            n_lower_tried++;
            int ok = make_lower(B, idx, lw, &idp, &trp);
            n_id += idp; n_true += trp; n_id_not_true += (idp && !trp); n_true_not_id += (!idp && trp);
            if (!ok) continue;
            n_lower_ok++;
            uint64_t M1[16], M2[16];
            for (int i = 0; i < 9; i++) M1[i] = lw[i];
            M1[9] = C_W1_9; M1[10] = C_W1_10; M1[11] = C_W1_11; M1[12] = C_W1_12; M1[13] = C_W1_13;
            // steps 0..13 once per lower (they only use W0..W13)
            uint64_t s[8];
            for (int i = 0; i < 8; i++) s[i] = B->cv[i];
            for (int t = 0; t < 14; t++) {
                uint64_t w2 = M1[t] + DM[t];
                if (!check_step(s, M1[t], w2, t)) l_low_fail[t]++;
            }
            for (int t = 0; t < 14; t++) l_pass[t]++;   // placeholder per-lower count (fail counted separately)
            uint64_t s13[8]; memcpy(s13, s, sizeof(s));
            for (int u = 0; u < uppers_per_lower; u++) {
                uint64_t r = splitmix64(&rng);
                int m = (int)(r & 15);
                uint64_t a15 = A15[m][(r >> 4) % nA15[m]];
                M1[14] = h_w14[m]; M1[15] = a15 - h_w15_sub[m];
                for (int i = 0; i < 16; i++) M2[i] = M1[i] + DM[i];
                uint64_t W1[36], W2[36];
                expand(M1, W1); expand(M2, W2);
                uint64_t st[8]; memcpy(st, s13, sizeof(st));
                uint64_t mask = 0;
                for (int t = 14; t < 36; t++) if (check_step(st, W1[t], W2[t], t)) mask |= (1ULL << t);
                mine++;
                l_pm_pairs[m]++;
                for (int t = 14; t < 36; t++) if ((mask >> t) & 1) l_pass[t]++;
                int c23 = (mask >> 23) & 1;
                if (c23) {
                    n23++; l_pm_23[m]++;
                    for (int t = 14; t < 36; t++) if ((mask >> t) & 1) l_and23[t]++;
                    if ((mask >> 25) & 1) {
                        n2325++; l_pm_2325[m]++;
                        uint64_t need = 0; for (int t = 24; t < 36; t++) need |= (1ULL << t);
                        if ((mask & need) == need) { ntail++; l_pm_tail[m]++; }
                        else for (int t = 24; t < 36; t++) if (!((mask >> t) & 1)) l_tailfail[t]++;
                    }
                }
                uint64_t J = (1ULL<<19)|(1ULL<<20)|(1ULL<<21)|(1ULL<<22);
                if ((mask & J) == J) { nJ++; if (c23) nJ23++; }
                uint64_t J2 = (1ULL<<16)|(1ULL<<17);
                if ((mask & J2) == J2) { nJ2++; if (c23) nJ223++; }
                uint64_t Ja = (1ULL<<19)|(1ULL<<20), Jb = (1ULL<<21)|(1ULL<<22);
                if ((mask & Ja) == Ja) { nJa++; if (c23) nJa23++; }
                if ((mask & Jb) == Jb) { nJb++; if (c23) nJb23++; }
            }
        }
        n_pairs += mine;
        #pragma omp critical
        {
            for (int t = 0; t < 36; t++) { pass[t] += l_pass[t]; and23[t] += l_and23[t]; tailfail[t] += l_tailfail[t]; }
            for (int t = 0; t < 16; t++) low_fail[t] += l_low_fail[t];
            for (int m = 0; m < 16; m++) { pm_pairs[m] += l_pm_pairs[m]; pm_23[m] += l_pm_23[m]; pm_2325[m] += l_pm_2325[m]; pm_tail[m] += l_pm_tail[m]; }
        }
    }
    double dt = omp_get_wtime() - t0;
    printf("\n[sample] %llu pairs from %llu valid lowers (%llu tried, pass %.3f%%) in %.1fs (%.1f M pairs/s)\n",
           (unsigned long long)n_pairs, (unsigned long long)n_lower_ok, (unsigned long long)n_lower_tried,
           100.0 * n_lower_ok / n_lower_tried, dt, n_pairs / dt / 1e6);
    printf("[lower] W5 4-bit identity pass=%llu, true ds0(W5) pass=%llu, identity&!true=%llu, true&!identity=%llu\n",
           (unsigned long long)n_id, (unsigned long long)n_true, (unsigned long long)n_id_not_true, (unsigned long long)n_true_not_id);
    printf("[lower] char steps 0..13 failures among valid lowers:");
    int anyf = 0; for (int t = 0; t < 14; t++) if (low_fail[t]) { printf(" t%d=%llu", t, (unsigned long long)low_fail[t]); anyf = 1; }
    printf("%s\n", anyf ? "" : " none");

    printf("\n[steps] unconditional per-step pass rate on real pairs (bits):\n");
    for (int t = 14; t < 36; t++)
        printf("  step %2d: %.4f%s\n", t, pass[t] ? -log2((double)pass[t] / n_pairs) : 99.0, t == 23 ? "   <-- GPU 'W21' counter (dW23 check)" : "");
    double b23 = -log2((double)n23 / n_pairs), se23 = 1.0 / (sqrt((double)n23) * log(2.0));
    printf("\n[row23] P(c23) = %llu / %llu = 2^-%.4f +/- %.4f   (model: 2^-11.5114)\n",
           (unsigned long long)n23, (unsigned long long)n_pairs, b23, se23);
    printf("[row23|t] conditional P(c23 | c_t) in bits (independence => ~%.3f):\n", b23);
    for (int t = 14; t < 23; t++)
        if (pass[t] && and23[t]) printf("  given step %2d: 2^-%.4f  (n_t=%llu, hits=%llu, se %.3f)\n", t,
            -log2((double)and23[t] / pass[t]), (unsigned long long)pass[t], (unsigned long long)and23[t], 1.0 / (sqrt((double)and23[t]) * log(2.0)));
    if (nJ) printf("  given steps 19..22 jointly: %llu / %llu = 2^-%.4f (se %.3f)\n", (unsigned long long)nJ23, (unsigned long long)nJ,
        nJ23 ? -log2((double)nJ23 / nJ) : 99.0, nJ23 ? 1.0 / (sqrt((double)nJ23) * log(2.0)) : 0.0);
    if (nJ2) printf("  given steps 16..17 jointly: %llu / %llu = 2^-%.4f (se %.3f)\n", (unsigned long long)nJ223, (unsigned long long)nJ2,
        nJ223 ? -log2((double)nJ223 / nJ2) : 99.0, nJ223 ? 1.0 / (sqrt((double)nJ223) * log(2.0)) : 0.0);
    if (nJa) printf("  given steps 19&20: %llu / %llu = 2^-%.4f (se %.3f)\n", (unsigned long long)nJa23, (unsigned long long)nJa,
        nJa23 ? -log2((double)nJa23 / nJa) : 99.0, nJa23 ? 1.0 / (sqrt((double)nJa23) * log(2.0)) : 0.0);
    if (nJb) printf("  given steps 21&22: %llu / %llu = 2^-%.4f (se %.3f)\n", (unsigned long long)nJb23, (unsigned long long)nJb,
        nJb23 ? -log2((double)nJb23 / nJb) : 99.0, nJb23 ? 1.0 / (sqrt((double)nJb23) * log(2.0)) : 0.0);
    printf("  joint steps 19..22 rate: 2^-%.4f (model 16.19 incl. E17 carries)\n", nJ ? -log2((double)nJ / n_pairs) : 99.0);
    printf("\n[row25] P(c25 | c23) = %llu / %llu = %.5f  (model 9/128 = %.5f)\n",
           (unsigned long long)n2325, (unsigned long long)n23, (double)n2325 / n23, 9.0 / 128.0);
    printf("[tail] P(steps 24..35 all pass | c23 & c25) = %llu / %llu\n", (unsigned long long)ntail, (unsigned long long)n2325);
    for (int t = 24; t < 36; t++) if (tailfail[t]) printf("  step %d failed %llu times\n", t, (unsigned long long)tailfail[t]);
    printf("\n[per-m] m: pairs  c23_bits  P(c25|c23)  tail_ok/(c23&c25)\n");
    for (int m = 0; m < 16; m++)
        printf("  %2d: %10llu  %7.3f  %7.4f  %llu/%llu\n", m, (unsigned long long)pm_pairs[m],
               pm_23[m] ? -log2((double)pm_23[m] / pm_pairs[m]) : 99.0, pm_23[m] ? (double)pm_2325[m] / pm_23[m] : 0.0,
               (unsigned long long)pm_tail[m], (unsigned long long)pm_2325[m]);
    return 0;
}
