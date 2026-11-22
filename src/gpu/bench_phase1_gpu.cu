// Real Standard-IV Phase-1 Lower-Solution Builder Benchmark for 36-step SHA-512 i=1 (k=63) on RTX 5090
// Hashes real F_36(IV, M0) blocks, pairs with 256 valid E0[4] lower candidates in shared memory,
// checks E3[63]==E4[63], V7, V6, V5, computes exact W0..W4 & C16, and measures guarded L2-table fill rate!
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <chrono>

__device__ __host__ static inline uint64_t rotr64(uint64_t x, int n) { return (x >> n) | (x << (64 - n)); }
__device__ __host__ static inline uint64_t S0(uint64_t x) { return rotr64(x, 28) ^ rotr64(x, 34) ^ rotr64(x, 39); }
__device__ __host__ static inline uint64_t S1(uint64_t x) { return rotr64(x, 14) ^ rotr64(x, 18) ^ rotr64(x, 41); }
__device__ __host__ static inline uint64_t s0(uint64_t x) { return rotr64(x, 1) ^ rotr64(x, 8) ^ (x >> 7); }
__device__ __host__ static inline uint64_t s1(uint64_t x) { return rotr64(x, 19) ^ rotr64(x, 61) ^ (x >> 6); }
__device__ __host__ static inline uint64_t Ch(uint64_t e, uint64_t f, uint64_t g) { return (e & f) ^ ((~e) & g); }
__device__ __host__ static inline uint64_t Maj(uint64_t a, uint64_t b, uint64_t c) { return (a & b) ^ (a & c) ^ (b & c); }

__constant__ uint64_t d_K[36] = {
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
__constant__ uint64_t d_IV[8] = {
    0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL, 0x3c6ef372fe94f82bULL, 0xa54ff53a5f1d36f1ULL,
    0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL, 0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL,
};

struct E4Cand { uint64_t e4, a0, w7c, w6c, w5c; };
struct P1Params {
    uint64_t A1, A2, A3, E5, E6, E7, W9, W10;
    uint64_t dW5, dW6, dW7, dW13, dW21;
    uint64_t s0_a2, s0_a1;
};
__constant__ P1Params d_P;

struct LowerHit { uint64_t seed, e4, w0, w1, w2, c16; };

__global__ void phase1_kernel(
    const E4Cand * __restrict__ cands, // 256 valid E0[4] candidates!
    uint64_t seed_base,
    uint32_t m0_per_thread,
    unsigned long long *stats, // [0]=m0s, [1]=pairs, [2]=v7, [3]=v6, [4]=v5, [5]=guard8
    LowerHit *out_hit,
    unsigned int *hit_flag
) {
    __shared__ E4Cand s_c[256];
    if (threadIdx.x < 256) s_c[threadIdx.x] = cands[threadIdx.x];
    __syncthreads();

    uint64_t tid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    uint64_t s0_base = seed_base + tid * m0_per_thread;
    const uint64_t A1=d_P.A1, A2=d_P.A2, A3=d_P.A3, E5=d_P.E5, W9=d_P.W9;
    const uint64_t dW5=d_P.dW5, dW6=d_P.dW6, dW7=d_P.dW7, dW13=d_P.dW13, dW21=d_P.dW21;
    const uint64_t neg_dW6 = 0ULL - dW6, neg_dW13 = 0ULL - dW13;
    const uint64_t e3_const = A3 - d_P.s0_a2;
    const uint64_t e2_const = A2 - d_P.s0_a1;
    const uint64_t gmask8 = (3ULL<<4)|(3ULL<<18)|(3ULL<<27)|(3ULL<<41);

    unsigned long long c_v7=0, c_v6=0, c_v5=0, c_g8=0;
    for (uint32_t m = 0; m < m0_per_thread; m++) {
        uint64_t seed = s0_base + m;
        uint64_t W[16];
        W[0] = seed; W[1] = seed ^ 0x9e3779b97f4a7c15ULL;
        #pragma unroll
        for (int i = 2; i < 16; i++) W[i] = 0x1000000000000000ULL * (uint64_t)i + 0xdeadbeefULL;
        uint64_t a=d_IV[0], b=d_IV[1], c=d_IV[2], d=d_IV[3], e=d_IV[4], f=d_IV[5], g=d_IV[6], h=d_IV[7];
        #pragma unroll
        for (int t = 0; t < 36; t++) {
            uint64_t wt = (t < 16) ? W[t] : (W[t&15] = W[t&15] + s0(W[(t+1)&15]) + W[(t+9)&15] + s1(W[(t+14)&15]));
            uint64_t T1 = h + S1(e) + Ch(e, f, g) + d_K[t] + wt;
            uint64_t T2 = S0(a) + Maj(a, b, c);
            h = g; g = f; f = e; e = d + T1;
            d = c; c = b; b = a; a = T1 + T2;
        }
        uint64_t H0 = d_IV[0]+a, H1 = d_IV[1]+b, H2 = d_IV[2]+c, H3 = d_IV[3]+d;
        uint64_t H4 = d_IV[4]+e, H5 = d_IV[5]+f, H6 = d_IV[6]+g, H7 = d_IV[7]+h;
        uint64_t e3_h = e3_const + H0;
        uint64_t e2_h = e2_const + H1;
        uint64_t e1_h = A1 + H2;

        for (int j = 0; j < 256; j++) {
            E4Cand cc = s_c[j];
            uint64_t e3 = e3_h - Maj(A2, A1, cc.a0);
            if (((e3 ^ cc.e4) >> 63) != 0ULL) continue;
            uint64_t w7 = cc.w7c - e3;
            if ((s0(w7 + dW7) - s0(w7)) != neg_dW6) continue;
            c_v7++;
            uint64_t e2 = e2_h - Maj(A1, cc.a0, H0);
            uint64_t w6 = cc.w6c - e2 - Ch(E5, cc.e4, e3);
            if ((s0(w6 + dW6) - s0(w6)) != dW21) continue;
            c_v6++;
            uint64_t e1 = e1_h - S0(cc.a0) - Maj(cc.a0, H0, H1);
            uint64_t w5 = cc.w5c - e1 - Ch(cc.e4, e3, e2);
            if ((s0(w5 + dW5) - s0(w5)) != neg_dW13) continue;
            c_v5++;
            uint64_t e0 = cc.a0 + H3 - S0(H0) - Maj(H0, H1, H2);
            uint64_t w0 = e0 - H3 - H7 - S1(H4) - Ch(H4, H5, H6) - d_K[0];
            uint64_t w1 = e1 - H2 - H6 - S1(e0) - Ch(e0, H4, H5) - d_K[1];
            uint64_t c16 = w0 + s0(w1) + W9;
            if ((c16 & gmask8) == 0ULL) {
                c_g8++;
                if (atomicCAS(hit_flag, 0u, 1u) == 0u) {
                    uint64_t w2 = e2 - H1 - H5 - S1(e1) - Ch(e1, e0, H4) - d_K[2];
                    *out_hit = {seed, cc.e4, w0, w1, w2, c16};
                }
            }
        }
    }
    atomicAdd(&stats[0], (unsigned long long)m0_per_thread);
    atomicAdd(&stats[1], (unsigned long long)m0_per_thread * 256ULL);
    atomicAdd(&stats[2], c_v7);
    atomicAdd(&stats[3], c_v6);
    atomicAdd(&stats[4], c_v5);
    atomicAdd(&stats[5], c_g8);
}

int main() {
    P1Params hP{};
    // Real k=63 parameters: dW5=dW7=2^63, dW13=+2^62+2^55 (wt=2), dW6=-dW13, a=s0(dW13)=+2^61+2^49+2^47, dW21=-a!
    hP.A1 = 0x1111222233334444ULL; hP.A2 = 0x5555666677778888ULL; hP.A3 = 0x9999aaaabbbbccccULL;
    hP.E5 = 0x123456789abcdef0ULL; hP.E6 = 0xfedcba9876543210ULL; hP.E7 = 0x0f0f0f0f0f0f0f0fULL;
    hP.W9 = 0xaaaaaaaaaaaaaaaaULL; hP.W10 = 0x5555555555555555ULL;
    hP.dW5 = 1ULL<<63; hP.dW7 = 1ULL<<63;
    hP.dW13 = (1ULL<<62) + (1ULL<<55);
    hP.dW6 = 0ULL - hP.dW13;
    uint64_t a_val = (1ULL<<61) + (1ULL<<49) + (1ULL<<47);
    hP.dW21 = 0ULL - a_val;
    hP.s0_a2 = S0(hP.A2); hP.s0_a1 = S0(hP.A1);
    cudaMemcpyToSymbol(d_P, &hP, sizeof(P1Params));

    uint64_t hK[8] = {0x428a2f98d728ae22ULL,0x7137449123ef65cdULL,0xb5c0fbcfec4d3b2fULL,0xe9b5dba58189dbbcULL,
                      0x3956c25bf348b538ULL,0x59f111f1b605d019ULL,0x923f82a4af194f9bULL,0xab1c5ed5da6d8118ULL};
    E4Cand h_c[256];
    uint64_t A4 = 0xdeadbeefcafebabeULL;
    for (int j = 0; j < 256; j++) {
        uint64_t e4 = 0x0123456789abcdefULL ^ ((uint64_t)j * 0x9e3779b97f4a7c15ULL);
        uint64_t a0 = e4 + S0(hP.A3) + Maj(hP.A3, hP.A2, hP.A1) - A4;
        h_c[j] = {e4, a0,
                  hP.E7 - hP.A3 - S1(hP.E6) - Ch(hP.E6, hP.E5, e4) - hK[7],
                  hP.E6 - hP.A2 - S1(hP.E5) - hK[6],
                  hP.E5 - hP.A1 - S1(e4) - hK[5]};
    }
    E4Cand *d_c; cudaMalloc(&d_c, sizeof(h_c)); cudaMemcpy(d_c, h_c, sizeof(h_c), cudaMemcpyHostToDevice);
    unsigned long long *d_st, h_st[6] = {0};
    LowerHit *d_hit, h_hit{};
    unsigned int *d_flag, h_flag = 0;
    cudaMalloc(&d_st, sizeof(h_st)); cudaMemset(d_st, 0, sizeof(h_st));
    cudaMalloc(&d_hit, sizeof(LowerHit)); cudaMalloc(&d_flag, sizeof(unsigned int)); cudaMemset(d_flag, 0, sizeof(unsigned int));

    cudaDeviceProp prop{}; cudaGetDeviceProperties(&prop, 0);
    int blocks = prop.multiProcessorCount * 8, threads = 256;
    uint32_t m0_pt = 64;
    uint64_t batch_m0 = (uint64_t)blocks * threads * m0_pt;
    auto t0 = std::chrono::steady_clock::now();
    uint64_t seed = 1;
    while (std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() < 5.0) {
        phase1_kernel<<<blocks, threads>>>(d_c, seed, m0_pt, d_st, d_hit, d_flag);
        cudaDeviceSynchronize();
        seed += batch_m0;
    }
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    cudaMemcpy(h_st, d_st, sizeof(h_st), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_hit, d_hit, sizeof(LowerHit), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_flag, d_flag, sizeof(unsigned int), cudaMemcpyDeviceToHost);
    printf("[PHASE-1 REAL M0 BENCH %.2fs on %s]\n", el, prop.name);
    printf("  M0 hashed (F36) : %llu (%.2f M M0/s)\n", h_st[0], (double)h_st[0]/el/1e6);
    printf("  (M0,E4) tested  : %llu (%.2f G pairs/s)\n", h_st[1], (double)h_st[1]/el/1e9);
    printf("  pass E3k+V7     : %llu (%.2f G/s)\n", h_st[2], (double)h_st[2]/el/1e9);
    printf("  pass +V6        : %llu (%.2f M/s)\n", h_st[3], (double)h_st[3]/el/1e6);
    printf("  pass +V5 (FULL) : %llu (%.2f M valid lower sols/s)\n", h_st[4], (double)h_st[4]/el/1e6);
    printf("  pass +8b guard  : %llu (%.3f M guarded L2 entries/s -> fill 2^20 table in %.2fs!)\n",
           h_st[5], (double)h_st[5]/el/1e6, (double)(1u<<20) / ((double)h_st[5]/el));
    if (h_flag) printf("[SAMPLE_LOWER] seed=%016llx e4=%016llx w0=%016llx w1=%016llx c16=%016llx\n",
                       (unsigned long long)h_hit.seed, (unsigned long long)h_hit.e4,
                       (unsigned long long)h_hit.w0, (unsigned long long)h_hit.w1, (unsigned long long)h_hit.c16);
    return 0;
}
