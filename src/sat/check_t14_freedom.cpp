// Quick 5-second probe of ALL degrees of freedom in 2026/1120 Table 14!
// Mode 1 (--enumAll): fix t=5..fixUpto, pin signDiffs on t=fixUpto+1..13, block full (W[t0..13]) tuple!
// Mode 2 (--freeE10_12): fix t=5..fixUpto, pin ONLY dA[10..12] (or only dA11=+2^59), leave E10..E12 diffs FREE, shrink |xE13|!
// Mode 3 (--bottom): fix t=10..15 (or 11..15) concrete, unfreeze t=5..9 (2026/1080 p.20 trick!) with signDiffs!
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
// Also t=1..4 concrete values for Mode 3 (unfreezing t=5..13!):
// Wait: let's compute REF_A1[0..4], REF_E1[0..4] directly from Table 14 CV & M1!
static const u64 REF_CV[8] = {
    0xf128698fe95be8acULL, 0x264034c8f44fe6edULL, 0x7b976a39ced01663ULL, 0x7b671a27b6f1da99ULL,
    0x198021f679057203ULL, 0x600969ddc8027d0cULL, 0x4af8d8a11db1692bULL, 0x98ba089bc0dee788ULL
};

struct Builder {
    CaDiCaL::Solver &S;
    int nv = 0, ncl = 0, zero_lit, one_lit;
    Builder(CaDiCaL::Solver &s) : S(s) { one_lit = new_var(); add({one_lit}); zero_lit = -one_lit; }
    int new_var() { return ++nv; }
    std::array<int, 64> new_word() { std::array<int,64> w; for (int i=0;i<64;i++) w[i]=new_var(); return w; }
    std::array<int, 64> const_word(u64 val) { std::array<int,64> w; for (int i=0;i<64;i++) w[i]=((val>>i)&1)?one_lit:zero_lit; return w; }
    void add(std::initializer_list<int> lits) { for (int l:lits) S.add(l); S.add(0); ncl++; }
    void eq(int a, int b) { add({-a,b}); add({a,-b}); }
    void eq_word(const std::array<int,64> &a, const std::array<int,64> &b) { for (int i=0;i<64;i++) eq(a[i],b[i]); }
    void pin_signed_diff(const std::array<int,64> &x1, const std::array<int,64> &x2, u64 r1, u64 r2) {
        u64 xr = r1 ^ r2;
        for (int b=0;b<64;b++) {
            if ((xr>>b)&1) { add({((r1>>b)&1)?x1[b]:-x1[b]}); add({((r2>>b)&1)?x2[b]:-x2[b]}); }
            else eq(x1[b], x2[b]);
        }
    }
    int xor2(int a, int b) { int z=new_var(); add({-a,-b,-z}); add({a,b,-z}); add({-a,b,z}); add({a,-b,z}); return z; }
    int xor3(int a, int b, int c) { int z=new_var(); add({-a,-b,-c,z}); add({-a,b,c,z}); add({a,-b,c,z}); add({a,b,-c,z}); add({-a,-b,c,-z}); add({-a,b,-c,-z}); add({a,-b,-c,-z}); add({a,b,c,-z}); return z; }
    int maj3(int a, int b, int c) { int z=new_var(); add({-a,-b,z}); add({-a,-c,z}); add({-b,-c,z}); add({a,b,-z}); add({a,c,-z}); add({b,c,-z}); return z; }
    int ch3(int x, int y, int z) { int r=new_var(); add({-x,-y,r}); add({-x,y,-r}); add({x,-z,r}); add({x,z,-r}); add({-y,-z,r}); add({y,z,-r}); return r; }
    std::array<int,64> w_S0(const std::array<int,64> &x) { std::array<int,64> r; for (int i=0;i<64;i++) r[i]=xor3(x[(i+28)&63],x[(i+34)&63],x[(i+39)&63]); return r; }
    std::array<int,64> w_S1(const std::array<int,64> &x) { std::array<int,64> r; for (int i=0;i<64;i++) r[i]=xor3(x[(i+14)&63],x[(i+18)&63],x[(i+41)&63]); return r; }
    std::array<int,64> w_s0(const std::array<int,64> &x) { std::array<int,64> r; for (int i=0;i<64;i++) r[i]=(i+7<64)?xor3(x[(i+1)&63],x[(i+8)&63],x[i+7]):xor2(x[(i+1)&63],x[(i+8)&63]); return r; }
    std::array<int,64> w_s1(const std::array<int,64> &x) { std::array<int,64> r; for (int i=0;i<64;i++) r[i]=(i+6<64)?xor3(x[(i+19)&63],x[(i+61)&63],x[i+6]):xor2(x[(i+19)&63],x[(i+61)&63]); return r; }
    std::array<int,64> w_Ch(const std::array<int,64> &x, const std::array<int,64> &y, const std::array<int,64> &z) { std::array<int,64> r; for (int i=0;i<64;i++) r[i]=ch3(x[i],y[i],z[i]); return r; }
    std::array<int,64> w_Maj(const std::array<int,64> &x, const std::array<int,64> &y, const std::array<int,64> &z) { std::array<int,64> r; for (int i=0;i<64;i++) r[i]=maj3(x[i],y[i],z[i]); return r; }
    std::array<int,64> add2(const std::array<int,64> &a, const std::array<int,64> &b) {
        std::array<int,64> s; int c=zero_lit;
        for (int i=0;i<64;i++) {
            if (i==0) { s[i]=xor2(a[i],b[i]); int nc=new_var(); add({-a[i],-b[i],nc}); add({a[i],-nc}); add({b[i],-nc}); c=nc; }
            else { s[i]=xor3(a[i],b[i],c); if (i<63) c=maj3(a[i],b[i],c); }
        }
        return s;
    }
    void csa3(const std::array<int,64> &a, const std::array<int,64> &b, const std::array<int,64> &c, std::array<int,64> &u, std::array<int,64> &v) {
        v[0]=zero_lit; for (int i=0;i<64;i++) { u[i]=xor3(a[i],b[i],c[i]); if (i<63) v[i+1]=maj3(a[i],b[i],c[i]); }
    }
    std::array<int,64> add_many(std::vector<std::array<int,64>> ws) {
        while (ws.size()>2) {
            std::vector<std::array<int,64>> nxt; size_t i=0;
            for (;i+2<ws.size();i+=3) { std::array<int,64> u,v; csa3(ws[i],ws[i+1],ws[i+2],u,v); nxt.push_back(u); nxt.push_back(v); }
            for (;i<ws.size();i++) nxt.push_back(ws[i]);
            ws = std::move(nxt);
        }
        return (ws.size()==1) ? ws[0] : add2(ws[0],ws[1]);
    }
    u64 val_word(const std::array<int,64> &w) { u64 r=0; for (int i=0;i<64;i++) if (S.val(w[i])>0) r|=(1ULL<<i); return r; }
};

int main(int argc, char **argv) {
    // Check freedom of the ENTIRE core t=5..15 when ALL signed diffs of Table 14 are pinned!
    // Wait: fix t=1..4 FREE (except E4..A4 concrete? No, let A1..A4, E1..E4 be FREE too — ONLY fix signed diffs of t=5..15!)
    // Mode A (--freeFullCore): pin ONLY signed diffs of Table 14 on t=5..15 + s0/s1 diffs on W5..W9,W13,W14!
    //   Measure how fast CaDiCaL generates COMPLETELY INDEPENDENT concrete cores T (all sharing Table 14's exact differential!)
    //   AND ALSO check if fixing t=5..8 concrete allows W9..W11 to vary!
    int mode = 0; // 0: fix t=5..8, enum ANY diff in W9..W13; 1: full free core (only signed diffs fixed!); 2: fix t=1..8, pin ONLY dA[9..12], free E10..E12!
    int enumK = 20;
    for (int i=1;i<argc;i++) {
        if (!strcmp(argv[i], "--fullCore")) mode = 1;
        else if (!strcmp(argv[i], "--freeE")) mode = 2;
        else if (!strncmp(argv[i], "--enum=", 7)) enumK = atoi(argv[i]+7);
    }

    // Compute exact REF_A1/E1 for t=-4..4
    u64 fA[21], fE[21]; // index t+4
    fA[3]=REF_CV[0]; fA[2]=REF_CV[1]; fA[1]=REF_CV[2]; fA[0]=REF_CV[3];
    fE[3]=REF_CV[4]; fE[2]=REF_CV[5]; fE[1]=REF_CV[6]; fE[0]=REF_CV[7];
    for (int t=0;t<16;t++) {
        u64 t1 = fE[t] + S1(fE[t+3]) + Ch(fE[t+3],fE[t+2],fE[t+1]) + K512[t] + REF_W1[t];
        u64 t2 = S0(fA[t+3]) + Maj(fA[t+3],fA[t+2],fA[t+1]);
        fA[t+4] = t1 + t2; fE[t+4] = fA[t] + t1;
    }

    CaDiCaL::Solver S; S.set("quiet", 1); S.set("seed", 42);
    Builder B(S);
    std::array<int,64> A1[17], A2[17], E1[17], E2[17], W1[16], W2[16];

    if (mode == 0) {
        // Fix t=5..8 concrete, pin signDiffs on t=9..14, enumerate ANY distinct (W9,W10,W11,W12,W13)!
        for (int t=5;t<=8;t++) {
            A1[t]=B.const_word(REF_A1[t]); A2[t]=B.const_word(REF_A2[t]);
            E1[t]=B.const_word(REF_E1[t]); E2[t]=B.const_word(REF_E2[t]);
            W1[t]=B.const_word(REF_W1[t]); W2[t]=B.const_word(REF_W2[t]);
        }
        for (int t=9;t<=15;t++) {
            W1[t]=B.new_word();
            W2[t]=(t==10||t==11||t==12||t==15)?W1[t]:B.add2(W1[t], B.const_word(REF_W2[t]-REF_W1[t]));
            auto T1_1=B.add_many({E1[t-4],B.w_S1(E1[t-1]),B.w_Ch(E1[t-1],E1[t-2],E1[t-3]),B.const_word(K512[t]),W1[t]});
            auto T2_1=B.add2(B.w_S0(A1[t-1]),B.w_Maj(A1[t-1],A1[t-2],A1[t-3]));
            A1[t]=B.add2(T1_1,T2_1); E1[t]=B.add2(A1[t-4],T1_1);
            auto T1_2=B.add_many({E2[t-4],B.w_S1(E2[t-1]),B.w_Ch(E2[t-1],E2[t-2],E2[t-3]),B.const_word(K512[t]),W2[t]});
            auto T2_2=B.add2(B.w_S0(A2[t-1]),B.w_Maj(A2[t-1],A2[t-2],A2[t-3]));
            A2[t]=B.add2(T1_2,T2_2); E2[t]=B.add2(A2[t-4],T1_2);
            B.pin_signed_diff(A1[t],A2[t],REF_A1[t],REF_A2[t]);
            B.pin_signed_diff(E1[t],E2[t],REF_E1[t],REF_E2[t]);
        }
        B.eq_word(B.w_s0(W2[9]), B.add2(B.w_s0(W1[9]), B.const_word(s0(REF_W2[9])-s0(REF_W1[9]))));
        B.eq_word(B.w_s0(W2[13]), B.add2(B.w_s0(W1[13]), B.const_word(s0(REF_W2[13])-s0(REF_W1[13]))));
        B.eq_word(B.w_s0(W2[14]), B.add2(B.w_s0(W1[14]), B.const_word(s0(REF_W2[14])-s0(REF_W1[14]))));
        B.eq_word(B.w_s1(W2[14]), B.add2(B.w_s1(W1[14]), B.const_word(s1(REF_W2[14])-s1(REF_W1[14]))));
        auto ch1=B.w_Ch(E1[15],E1[14],E1[13]), ch2=B.w_Ch(E2[15],E2[14],E2[13]);
        auto m1=B.w_Maj(A1[15],A1[14],A1[13]), m2=B.w_Maj(A2[15],A2[14],A2[13]);
        B.eq_word(B.add_many({E1[12],ch1,B.w_S0(A1[15]),m1}), B.add_many({E2[12],ch2,B.w_S0(A2[15]),m2}));
        B.eq_word(B.add_many({A1[12],E1[12],ch1,B.const_word(REF_E2[16]-REF_E1[16])}), B.add_many({A2[12],E2[12],ch2}));

        printf("[MODE 0: fix t=5..8, pin signDiffs, block (W9..W13)] vars=%d cls=%d\n", B.nv, B.ncl);
        for (int k=0;k<enumK;k++) {
            auto ts=std::chrono::steady_clock::now();
            int res=S.solve();
            double dt=std::chrono::duration<double>(std::chrono::steady_clock::now()-ts).count();
            if (res!=10) { printf("  Stopped at k=%d res=%d (%.3fs)\n", k, res, dt); break; }
            u64 w9=B.val_word(W1[9]), w10=B.val_word(W1[10]), w11=B.val_word(W1[11]), w12=B.val_word(W1[12]), w13=B.val_word(W1[13]);
            printf("  #%2d (%.3fs): W9=%016llx W10=%016llx W11=%016llx W12=%016llx W13=%016llx\n",
                   k+1, dt, (unsigned long long)w9, (unsigned long long)w10, (unsigned long long)w11, (unsigned long long)w12, (unsigned long long)w13);
            u64 ws[5]={w9,w10,w11,w12,w13}; for (int t=9;t<=13;t++) { u64 v=ws[t-9]; for (int b=0;b<64;b++) S.add(((v>>b)&1)?-W1[t][b]:W1[t][b]); }
            S.add(0);
        }
    } else if (mode == 1) {
        // FULL CORE t=5..15 with ONLY Table 14's exact signed diffs pinned (A1..A4, E1..E4 FREE!)
        // Wait: how fast does CaDiCaL generate a BRAND NEW concrete core T (matching 2026/1120 Step 1)?
        // Wait: in 2026/1120, they got 1 core in 2^{39.4} (~seconds-minutes!). Let's check how fast CaDiCaL solves it AND how many distinct (W5..W13) exist!
        for (int t=1;t<=4;t++) { A1[t]=B.new_word(); A2[t]=A1[t]; E1[t]=B.new_word(); E2[t]=E1[t]; }
        for (int t=5;t<=15;t++) {
            W1[t]=B.new_word();
            W2[t]=(t==10||t==11||t==12||t==15)?W1[t]:B.add2(W1[t], B.const_word(REF_W2[t]-REF_W1[t]));
            auto T1_1=B.add_many({E1[t-4],B.w_S1(E1[t-1]),B.w_Ch(E1[t-1],E1[t-2],E1[t-3]),B.const_word(K512[t]),W1[t]});
            auto T2_1=B.add2(B.w_S0(A1[t-1]),B.w_Maj(A1[t-1],A1[t-2],A1[t-3]));
            A1[t]=B.add2(T1_1,T2_1); E1[t]=B.add2(A1[t-4],T1_1);
            auto T1_2=B.add_many({E2[t-4],B.w_S1(E2[t-1]),B.w_Ch(E2[t-1],E2[t-2],E2[t-3]),B.const_word(K512[t]),W2[t]});
            auto T2_2=B.add2(B.w_S0(A2[t-1]),B.w_Maj(A2[t-1],A2[t-2],A2[t-3]));
            A2[t]=B.add2(T1_2,T2_2); E2[t]=B.add2(A2[t-4],T1_2);
            B.pin_signed_diff(A1[t],A2[t],REF_A1[t],REF_A2[t]);
            B.pin_signed_diff(E1[t],E2[t],REF_E1[t],REF_E2[t]);
        }
        for (int t : {5,6,7,8,9,13,14})
            B.eq_word(B.w_s0(W2[t]), B.add2(B.w_s0(W1[t]), B.const_word(s0(REF_W2[t])-s0(REF_W1[t]))));
        B.eq_word(B.w_s1(W2[14]), B.add2(B.w_s1(W1[14]), B.const_word(s1(REF_W2[14])-s1(REF_W1[14]))));
        auto ch1=B.w_Ch(E1[15],E1[14],E1[13]), ch2=B.w_Ch(E2[15],E2[14],E2[13]);
        auto m1=B.w_Maj(A1[15],A1[14],A1[13]), m2=B.w_Maj(A2[15],A2[14],A2[13]);
        B.eq_word(B.add_many({E1[12],ch1,B.w_S0(A1[15]),m1}), B.add_many({E2[12],ch2,B.w_S0(A2[15]),m2}));
        B.eq_word(B.add_many({A1[12],E1[12],ch1,B.const_word(REF_E2[16]-REF_E1[16])}), B.add_many({A2[12],E2[12],ch2}));

        printf("[MODE 1: FULL FREE CORE t=5..15, Table 14 signDiffs] vars=%d cls=%d\n", B.nv, B.ncl);
        for (int k=0;k<enumK;k++) {
            auto ts=std::chrono::steady_clock::now();
            int res=S.solve();
            double dt=std::chrono::duration<double>(std::chrono::steady_clock::now()-ts).count();
            if (res!=10) { printf("  Stopped at k=%d res=%d (%.3fs)\n", k, res, dt); break; }
            printf("  #%2d (%.3fs): A5=%016llx W9=%016llx W12=%016llx W13=%016llx\n",
                   k+1, dt, (unsigned long long)B.val_word(A1[5]), (unsigned long long)B.val_word(W1[9]),
                   (unsigned long long)B.val_word(W1[12]), (unsigned long long)B.val_word(W1[13]));
            u64 v=B.val_word(W1[12]); for (int b=0;b<64;b++) S.add(((v>>b)&1)?-W1[12][b]:W1[12][b]);
            S.add(0);
        }
    }
    return 0;
}
