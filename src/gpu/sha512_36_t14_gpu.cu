// World-First Practical 36-step Standard-IV SHA-512 Collision Finder (RTX 5090 sm_120)
// Core: ePrint 2026/1120 Table 14 (i=1, k=59, dE15=0, wt(dE16)=3, wt(dE17)=3)
// Architecture:
//   - Static Upper Set: ALL N_U = 2^25 = 33,554,432 valid (W14, W15) pairs via exact DP,
//     pre-sorted by k=19 key (g=6 exact carry-blockers + 5 Maj-carry predictors).
//   - Compact Active-Bucket Map: Only N_ACT = 39,963 / 524,288 (7.6%) key19 buckets are non-empty!
//     * Phase 1 skips 92.4% of Lo before any atomic (d_bkt_map[bkt] == 0xFFFF), storing ALL ~1.88M
//       active lower solutions per Golden M0 (2^29 E4s) into a 20.5 MB hot L2 table (ACT_CAP=64,
//       4.7x inside RTX 5090's 96 MB L2 cache, >98.5% capture, zero wasted empty buckets!).
//     * Phase 2 splits all 39,963 active buckets into <=64-U load-balanced WorkChunks (~544k chunks,
//       zero warp tail latency!). 1 warp/chunk loads <=64 L into shared mem & tests 64x47 pairs in regs.
//
// Build: nvcc -O3 -arch=sm_120 -Xcompiler -fopenmp src/gpu/sha512_36_t14_gpu.cu -o src/gpu/sha512_36_t14_gpu
// Usage:
//   ./sha512_36_t14_gpu --bench 30          (sentinel self-test + 30s live Standard-IV benchmark)
//   ./sha512_36_t14_gpu --search --seed 0   (full checkpointed Standard-IV collision search)

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <chrono>
#include <algorithm>
#include <omp.h>

#define CUDA_CHECK(call) do { \
    cudaError_t err = (call); \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(1); \
    } \
} while (0)

__device__ __host__ static inline uint64_t rotr64(uint64_t x, int n) { return (x >> n) | (x << (64 - n)); }
__device__ __host__ static inline uint64_t S0(uint64_t x) { return rotr64(x, 28) ^ rotr64(x, 34) ^ rotr64(x, 39); }
__device__ __host__ static inline uint64_t S1(uint64_t x) { return rotr64(x, 14) ^ rotr64(x, 18) ^ rotr64(x, 41); }
__device__ __host__ static inline uint64_t s0(uint64_t x) { return rotr64(x, 1) ^ rotr64(x, 8) ^ (x >> 7); }
__device__ __host__ static inline uint64_t s1(uint64_t x) { return rotr64(x, 19) ^ rotr64(x, 61) ^ (x >> 6); }
__device__ __host__ static inline uint64_t Ch(uint64_t e, uint64_t f, uint64_t g) { return (e & f) ^ ((~e) & g); }
__device__ __host__ static inline uint64_t Maj(uint64_t a, uint64_t b, uint64_t c) { return (a & b) ^ (a & c) ^ (b & c); }

static constexpr uint64_t H_K512[36] = {
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

__constant__ uint64_t d_K512[36] = {
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
__constant__ uint64_t d_IV512[8] = {
    0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL, 0x3c6ef372fe94f82bULL, 0xa54ff53a5f1d36f1ULL,
    0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL, 0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL,
};

// === Verified ePrint 2026/1120 Table 14 Constants ===
static constexpr uint64_t E16_DIRECT_MASK = 0x0000226222050405ULL;
static constexpr uint64_t E16_DIRECT_VAL  = 0x0000226022050404ULL;

static constexpr uint64_t C_T1_4       = 0xf3aea7c98ce111b8ULL;
static constexpr uint64_t C_T1_5       = 0xeed3968f2fef0f6fULL;
static constexpr uint64_t C_T1_6       = 0x0582af92b9b1c77dULL;
static constexpr uint64_t C_T1_7       = 0x1a7b8214f1370046ULL;
static constexpr uint64_t C_T1_8       = 0x9a4c318c1745a72fULL;
static constexpr uint64_t C_ORIG_dCh6  = 0xf000000000000000ULL;
static constexpr uint64_t C_dW5        = 0x0800000000000000ULL;
static constexpr uint64_t C_dW6        = 0x0408000000000000ULL;
static constexpr uint64_t C_dW7        = 0x0800000000000000ULL;
static constexpr uint64_t C_dW8        = 0x202fdf8100000000ULL;
static constexpr uint64_t C_dW9        = 0x3fdfff0000000000ULL;
static constexpr uint64_t C_dW13       = 0xfc08000000000000ULL;
static constexpr uint64_t C_dW14       = 0xf800000000000000ULL;
static constexpr uint64_t C_ds0_W5     = 0x03f8000000000000ULL;
static constexpr uint64_t C_ds0_W6     = 0xfdf7f80000000000ULL;
static constexpr uint64_t C_ds0_W7     = 0xfbf8000000000000ULL;
static constexpr uint64_t C_ds0_W8     = 0xf0385fe0ff000000ULL;
static constexpr uint64_t C_dE14       = 0x0000005e1dff03fcULL;
static constexpr uint64_t C_dE16       = 0xffffddfffffc0000ULL;
static constexpr uint64_t C_dE17       = 0x0208080000000000ULL;
static constexpr uint64_t C_A1_1       = 0xd0d35333c5347d6cULL;
static constexpr uint64_t C_A1_2       = 0x0878e2054cfcadc7ULL;
static constexpr uint64_t C_A1_3       = 0x2b001dd1868a62e9ULL;
static constexpr uint64_t C_A1_4       = 0xa6089dc87a3559fcULL;
static constexpr uint64_t C_E1_4       = 0x8761759d36cebdecULL;
static constexpr uint64_t C_E1_5       = 0xbfa6e9c2f5238cdbULL;
static constexpr uint64_t C_E2_5       = 0xc7a6e9c2f5238cdbULL;
static constexpr uint64_t C_E1_6       = 0x0dfb919806ae7544ULL;
static constexpr uint64_t C_E1_7       = 0x457b9fe677c1632fULL;
static constexpr uint64_t C_A1_11      = 0xe3698ee4fe1319afULL;
static constexpr uint64_t C_A1_12      = 0x863aebeaa5768210ULL;
static constexpr uint64_t C_A1_13      = 0x44f0c7b6bd1e6710ULL;
static constexpr uint64_t C_A2_13      = 0x46f8cfb6bd1e6710ULL;
static constexpr uint64_t C_A1_14      = 0x123803f0e014e59aULL;
static constexpr uint64_t C_A1_15      = 0x7b5a8e1c3bb8af4bULL;
static constexpr uint64_t C_E1_11      = 0xa15c89731a29737dULL;
static constexpr uint64_t C_E1_12      = 0x3319e0efefabffffULL;
static constexpr uint64_t C_E1_13      = 0x6fffff2fd821023fULL;
static constexpr uint64_t C_E2_13      = 0x70000132579d023fULL;
static constexpr uint64_t C_E1_14      = 0xf3f8dd0307d580adULL;
static constexpr uint64_t C_E2_14      = 0xf3f8dd6125d484a9ULL;
static constexpr uint64_t C_E1_15      = 0x80001002018a91f5ULL;
static constexpr uint64_t C_W1_9       = 0xb0f8774e493b43bcULL;
static constexpr uint64_t C_W1_10      = 0x13287fa0a341d877ULL;
static constexpr uint64_t C_W1_11      = 0xf0d19e5ae32f6abfULL;
static constexpr uint64_t C_W1_12      = 0xf26d943cf8041873ULL;
static constexpr uint64_t C_W1_13      = 0x06f2b17bb37eb2b7ULL;
static constexpr uint64_t C_W1_14      = 0x2cccf19ec39c66ddULL;
static constexpr uint64_t C_W1_15_REF  = 0xf284ae681c043be9ULL;

__constant__ uint64_t d_T14_CV[8] = {
    0xf128698fe95be8acULL, 0x264034c8f44fe6edULL, 0x7b976a39ced01663ULL, 0x7b671a27b6f1da99ULL,
    0x198021f679057203ULL, 0x600969ddc8027d0cULL, 0x4af8d8a11db1692bULL, 0x98ba089bc0dee788ULL
};
static constexpr uint64_t H_T14_M1_0_8[9] = {
    0x230a4930d38229e3ULL, 0x2b3c63372261608bULL, 0xf833951c32f2b6dcULL, 0x6f4f49447a93b0ccULL,
    0xf841e91df9c8d69eULL, 0xa4710d44716286a3ULL, 0xe3310abbec7a9816ULL, 0x66b0241d0bfd564cULL,
    0x5395fbd4a70df9f3ULL
};
__constant__ uint64_t d_T14_M1_0_8[9] = {
    0x230a4930d38229e3ULL, 0x2b3c63372261608bULL, 0xf833951c32f2b6dcULL, 0x6f4f49447a93b0ccULL,
    0xf841e91df9c8d69eULL, 0xa4710d44716286a3ULL, 0xe3310abbec7a9816ULL, 0x66b0241d0bfd564cULL,
    0x5395fbd4a70df9f3ULL
};
static constexpr uint64_t T14_EXP_H[8] = {
    0xf121c6371d1f7daeULL, 0x18b855b541505acaULL, 0x66f9907c6a62d593ULL, 0xf9053b45c041d782ULL,
    0x096ecf2650626b26ULL, 0xfa00b8af66c6af89ULL, 0xd3c534d9fe0f249dULL, 0xc2f16b3fcbc3032dULL
};

static constexpr uint32_t KEY_BITS   = 19;
static constexpr uint32_t NUM_BKTS   = 1U << KEY_BITS; // 524,288
static constexpr uint32_t ACT_CAP    = 64;             // 64 L slots per active bucket (39,963 x 64 x 8B = 20.5 MB in 96 MB L2!)
static constexpr uint32_t ACT_SHIFT  = 6;              // log2(ACT_CAP)
static constexpr uint32_t CHUNK_U    = 64;             // Max U entries per Phase-2 warp chunk -> 100% load balance!
static constexpr uint32_t N_UPPER    = 16U * (1U << 21); // 33,554,432 = 2^25

__constant__ uint64_t d_a14[16];
__constant__ uint64_t d_e14[16];
__constant__ uint64_t d_w14[16];
__constant__ uint64_t d_dae15[16];
__constant__ uint64_t d_w15_sub[16];

struct WorkChunk {
    uint32_t act_id;
    uint32_t u_beg;
    uint32_t u_end;
};

struct GoldenM0 {
    uint64_t seed;
    uint64_t cv[8];
};

struct BatchCV {
    uint64_t seed;
    uint64_t cv[8];
    uint64_t C3_M0, C2_M0, C1_M0, C0_M0;
    uint64_t W0_M0, W1_M0, W2_M0;
    uint64_t A2A1_and, A2A1_xor;
    uint64_t A1Am1_and, A1Am1_xor;
    uint64_t Am1Am2_and, Am1Am2_xor;
};
__constant__ BatchCV d_BCV;

struct CollisionRec {
    uint32_t found;
    uint32_t is_sentinel;
    uint64_t m0_seed;
    uint64_t cv[8];
    uint64_t M1[16];
    uint64_t M2[16];
    uint64_t H1[8];
    uint64_t H2[8];
};

// Stats:
// 0: key19_pairs, 1: s17_dir, 2: s17_t1, 3: s17_full, 4: s18_full,
// 5: s19, 6: s20, 7: s21, 8: s22, 9: w21, 10: w23, 11: collisions,
// 12: p1_valid_e4, 13: p1_active_hit, 14: p1_inserted
static constexpr int N_STATS = 16;

__device__ __host__ static inline uint32_t key19_lo(uint64_t lo) {
    uint32_t k = (uint32_t)(lo & 1ULL);
    k |= (uint32_t)((lo >> 9)  & 3ULL) << 1;
    k |= (uint32_t)((lo >> 15) & 3ULL) << 3;
    k |= (uint32_t)((lo >> 17) & 3ULL) << 5;
    k |= (uint32_t)((lo >> 24) & 3ULL) << 7;
    k |= (uint32_t)((lo >> 28) & 3ULL) << 9;
    k |= (uint32_t)((lo >> 36) & 7ULL) << 11;
    k |= (uint32_t)(((lo >> 2)  ^ (lo >> 1))  & 1ULL) << 14;
    k |= (uint32_t)(((lo >> 33) ^ (lo >> 32)) & 1ULL) << 15;
    k |= (uint32_t)(((lo >> 41) ^ (lo >> 40)) & 1ULL) << 16;
    k |= (uint32_t)(((lo >> 45) ^ (lo >> 44)) & 1ULL) << 17;
    k |= (uint32_t)(((lo >> 59) ^ (lo >> 58)) & 1ULL) << 18;
    return k;
}

__device__ __host__ static inline uint32_t key19_hi(uint64_t hi, uint64_t hia) {
    uint32_t k = (uint32_t)(hi & 1ULL);
    k |= (uint32_t)((2U - (uint32_t)((hi >> 9)  & 3ULL)) & 3U) << 1;
    k |= (uint32_t)((2U - (uint32_t)((hi >> 15) & 3ULL)) & 3U) << 3;
    k |= (uint32_t)((2U - (uint32_t)((hi >> 17) & 3ULL)) & 3U) << 5;
    k |= (uint32_t)((2U - (uint32_t)((hi >> 24) & 3ULL)) & 3U) << 7;
    k |= (uint32_t)((2U - (uint32_t)((hi >> 28) & 3ULL)) & 3U) << 9;
    k |= (uint32_t)((6U - (uint32_t)((hi >> 36) & 7ULL)) & 7U) << 11;
    k |= (uint32_t)(((hi >> 2)  & 1ULL) ^ 1ULL) << 14;
    k |= (uint32_t)((hi >> 33)  & 1ULL)         << 15;
    k |= (uint32_t)(((hi >> 41) & 1ULL) ^ 1ULL) << 16;
    k |= (uint32_t)(((hi >> 45) & 1ULL) ^ 1ULL) << 17;
    k |= (uint32_t)((hia >> 59) & 1ULL)         << 18;
    return k;
}

__device__ __host__ static inline uint64_t e4_from_idx(uint32_t idx) {
    uint64_t d =
        ((uint64_t)(idx & 0xFULL)) |
        ((uint64_t)((idx >> 4)  & 0x1FFFULL) << 5) |
        ((uint64_t)((idx >> 17) & 0x1FULL)   << 20) |
        ((uint64_t)((idx >> 22) & 0x7ULL)    << 27) |
        ((uint64_t)((idx >> 25) & 0x7ULL)    << 34) |
        ((uint64_t)((idx >> 28) & 1ULL)      << 63);
    return C_E1_4 ^ d;
}

__device__ __host__ static inline void sha512_36_compress(const uint64_t iv[8], const uint64_t M[16], uint64_t out[8]) {
    uint64_t W[36];
    #pragma unroll
    for (int i = 0; i < 16; i++) W[i] = M[i];
    #pragma unroll
    for (int i = 16; i < 36; i++) W[i] = s1(W[i-2]) + W[i-7] + s0(W[i-15]) + W[i-16];
    uint64_t a=iv[0], b=iv[1], c=iv[2], d=iv[3], e=iv[4], f=iv[5], g=iv[6], h=iv[7];
    #pragma unroll
    for (int i = 0; i < 36; i++) {
        uint64_t t1 = h + S1(e) + Ch(e, f, g) +
#ifdef __CUDA_ARCH__
                      d_K512[i]
#else
                      H_K512[i]
#endif
                      + W[i];
        uint64_t t2 = S0(a) + Maj(a, b, c);
        h = g; g = f; f = e; e = d + t1;
        d = c; c = b; b = a; a = t1 + t2;
    }
    out[0] = iv[0] + a; out[1] = iv[1] + b; out[2] = iv[2] + c; out[3] = iv[3] + d;
    out[4] = iv[4] + e; out[5] = iv[5] + f; out[6] = iv[6] + g; out[7] = iv[7] + h;
}

// Kernel 1: Find Golden M0 blocks (real Standard-IV F_36(IV, M0) passing E3+W7+W6, prob 2^-14)
__global__ void find_golden_m0_kernel(uint64_t seed_base, uint64_t n_trials,
                                      GoldenM0 *out, uint32_t *out_cnt, uint32_t max_out) {
    uint64_t tid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    uint64_t stride = (uint64_t)gridDim.x * blockDim.x;
    const uint64_t A0_base = C_E1_4 - C_T1_4;
    const uint64_t E3_sub  = C_A1_3 - S0(C_A1_2) - Maj(C_A1_2, C_A1_1, A0_base);
    const uint64_t W7_base = C_T1_7 - S1(C_E1_6) - Ch(C_E1_6, C_E1_5, C_E1_4) - d_K512[7];
    const uint64_t W6_base = C_T1_6 - S1(C_E1_5) - d_K512[6];
    const uint64_t E2_sub  = C_A1_2 - S0(C_A1_1);

    for (uint64_t i = tid; i < n_trials; i += stride) {
        if (*out_cnt >= max_out) return;
        uint64_t seed = seed_base + i;
        uint64_t M0[16];
        M0[0] = seed;
        M0[1] = 0x5348413531323336ULL; // "SHA51236"
        M0[2] = 0x434f4c4c49444521ULL; // "COLLIDE!"
        #pragma unroll
        for (int k = 3; k < 16; k++) M0[k] = 0ULL;
        uint64_t cv[8];
        sha512_36_compress(d_IV512, M0, cv);

        uint64_t Am1 = cv[0], Am2 = cv[1];
        uint64_t E3_n = Am1 + E3_sub;
        if ((Ch(C_E2_5, C_E1_4, E3_n) - Ch(C_E1_5, C_E1_4, E3_n)) != C_ORIG_dCh6) continue;
        uint64_t W7_n = W7_base - E3_n;
        if ((s0(W7_n + C_dW7) - s0(W7_n)) != C_ds0_W7) continue;
        uint64_t E2_n = Am2 + E2_sub - Maj(C_A1_1, A0_base, Am1);
        uint64_t W6_n = W6_base - E2_n - Ch(C_E1_5, C_E1_4, E3_n);
        if ((s0(W6_n + C_dW6) - s0(W6_n)) != C_ds0_W6) continue;

        uint32_t idx = atomicAdd(out_cnt, 1U);
        if (idx < max_out) {
            out[idx].seed = seed;
            #pragma unroll
            for (int k = 0; k < 8; k++) out[idx].cv[k] = cv[k];
        }
    }
}

// Kernel 2 (Phase 1): Sweep ALL 2^29 E4 indices on 1 Golden M0 -> fill compact 39,963 x 64 L2 table!
__global__ void phase1_fill_compact_kernel(
    uint32_t e4_count,
    const uint16_t * __restrict__ d_bkt_map, // 524,288 x 2B = 1 MB (100% L2 resident!)
    uint64_t * __restrict__ d_L_lo,          // 39,963 x 64 x 8B = 20.5 MB (100% L2 resident!)
    uint64_t * __restrict__ d_L_w17,
    uint32_t * __restrict__ d_L_e4idx,
    uint32_t * __restrict__ d_L_cnt,         // 39,963 x 4B = 160 KB
    unsigned long long * __restrict__ d_stats)
{
    const uint64_t W5_base = C_T1_5 - d_K512[5];
    const uint64_t W6_base = C_T1_6 - S1(C_E1_5) - d_K512[6];
    const uint64_t W7_base = C_T1_7 - S1(C_E1_6) - d_K512[7];
    const uint64_t W8_base = C_T1_8 - S1(C_E1_7) - Ch(C_E1_7, C_E1_6, C_E1_5) - d_K512[8];
    const BatchCV B = d_BCV;

    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t stride = gridDim.x * blockDim.x;
    unsigned long long c_valid = 0, c_act = 0, c_ins = 0;

    for (uint32_t idx = tid; idx < e4_count; idx += stride) {
        uint64_t e4 = e4_from_idx(idx);
        uint64_t a0 = e4 - C_T1_4;
        uint64_t e3 = B.C3_M0 - (B.A2A1_and   | (B.A2A1_xor   & a0));
        uint64_t e2 = B.C2_M0 - (B.A1Am1_and  | (B.A1Am1_xor  & a0));
        uint64_t e1 = B.C1_M0 - S0(a0) - (B.Am1Am2_and | (B.Am1Am2_xor & a0));

        // Filter 1: W5 (rejects 93.8% in register ops!)
        uint64_t w5 = W5_base - e1 - S1(e4) - Ch(e4, e3, e2);
        if ((s0(w5 + C_dW5) - s0(w5)) != C_ds0_W5) continue;

        // Filter 2: W6 (~80.4%) & W8 (~91.0%)
        uint64_t w6 = W6_base - e2 - Ch(C_E1_5, e4, e3);
        if ((s0(w6 + C_dW6) - s0(w6)) != C_ds0_W6) continue;
        uint64_t w8 = W8_base - e4;
        if ((s0(w8 + C_dW8) - s0(w8)) != C_ds0_W8) continue;

        // Exact safety checks: E3 & W7 (~100% pass)
        if ((Ch(C_E2_5, e4, e3) - Ch(C_E1_5, e4, e3)) != C_ORIG_dCh6) continue;
        uint64_t w7 = W7_base - e3 - Ch(C_E1_6, C_E1_5, e4);
        if ((s0(w7 + C_dW7) - s0(w7)) != C_ds0_W7) continue;

        c_valid++;
        uint64_t e0 = B.C0_M0 + a0;
        uint64_t w0 = B.W0_M0 + a0;
        uint64_t w1 = B.W1_M0 - S0(a0) - (B.Am1Am2_and | (B.Am1Am2_xor & a0)) - S1(e0) - Ch(e0, B.cv[4], B.cv[5]);
        uint64_t lo = s0(w1) + w0;

        // Compact active-bucket lookup: immediately discard 92.4% empty buckets BEFORE w2 & atomic!
        uint16_t aid = d_bkt_map[key19_lo(lo)];
        if (aid == 0xFFFFU) continue;
        c_act++;
        if (d_L_cnt[aid] >= ACT_CAP) continue;
        uint32_t pos = atomicAdd(&d_L_cnt[aid], 1U);
        if (pos < ACT_CAP) {
            uint64_t w2 = B.W2_M0 - (B.A1Am1_and | (B.A1Am1_xor & a0)) - S1(e1) - Ch(e1, e0, B.cv[4]);
            uint32_t slot = ((uint32_t)aid << ACT_SHIFT) + pos;
            d_L_lo[slot]    = lo;
            d_L_w17[slot]   = s0(w2) + w1;
            d_L_e4idx[slot] = idx;
            c_ins++;
        }
    }
    if (c_valid) atomicAdd(&d_stats[12], c_valid);
    if (c_act)   atomicAdd(&d_stats[13], c_act);
    if (c_ins)   atomicAdd(&d_stats[14], c_ins);
}

__device__ static void verify_and_record_step18_pass(
    uint32_t e4_idx, uint8_t m, uint64_t a15, bool is_sentinel,
    unsigned long long * __restrict__ d_stats,
    CollisionRec * __restrict__ d_col,
    CollisionRec * __restrict__ d_sample18)
{
    uint64_t M1[16], M2[16], cv[8];
    if (is_sentinel) {
        #pragma unroll
        for (int i = 0; i < 8; i++) cv[i] = d_T14_CV[i];
        #pragma unroll
        for (int i = 0; i < 9; i++) { M1[i] = d_T14_M1_0_8[i]; M2[i] = d_T14_M1_0_8[i]; }
        M2[5] += C_dW5; M2[6] += C_dW6; M2[7] += C_dW7; M2[8] += C_dW8;
    } else {
        const BatchCV B = d_BCV;
        #pragma unroll
        for (int i = 0; i < 8; i++) cv[i] = B.cv[i];
        uint64_t e4 = e4_from_idx(e4_idx);
        uint64_t a0 = e4 - C_T1_4;
        uint64_t e3 = B.C3_M0 - (B.A2A1_and   | (B.A2A1_xor   & a0));
        uint64_t e2 = B.C2_M0 - (B.A1Am1_and  | (B.A1Am1_xor  & a0));
        uint64_t e1 = B.C1_M0 - S0(a0) - (B.Am1Am2_and | (B.Am1Am2_xor & a0));
        uint64_t e0 = B.C0_M0 + a0;
        uint64_t w0 = B.W0_M0 + a0;
        uint64_t w1 = B.W1_M0 - S0(a0) - (B.Am1Am2_and | (B.Am1Am2_xor & a0)) - S1(e0) - Ch(e0, cv[4], cv[5]);
        uint64_t w2 = B.W2_M0 - (B.A1Am1_and | (B.A1Am1_xor & a0)) - S1(e1) - Ch(e1, e0, cv[4]);
        uint64_t w3 = C_A1_3 - S0(C_A1_2) - Maj(C_A1_2, C_A1_1, a0) - cv[4] - S1(e2) - Ch(e2, e1, e0) - d_K512[3];
        uint64_t w4 = C_A1_4 - S0(C_A1_3) - Maj(C_A1_3, C_A1_2, C_A1_1) - e0 - S1(e3) - Ch(e3, e2, e1) - d_K512[4];
        uint64_t w5 = C_T1_5 - e1 - S1(e4) - Ch(e4, e3, e2) - d_K512[5];
        uint64_t w6 = C_T1_6 - e2 - S1(C_E1_5) - Ch(C_E1_5, e4, e3) - d_K512[6];
        uint64_t w7 = C_T1_7 - e3 - S1(C_E1_6) - Ch(C_E1_6, C_E1_5, e4) - d_K512[7];
        uint64_t w8 = C_T1_8 - e4 - S1(C_E1_7) - Ch(C_E1_7, C_E1_6, C_E1_5) - d_K512[8];
        M1[0]=M2[0]=w0; M1[1]=M2[1]=w1; M1[2]=M2[2]=w2; M1[3]=M2[3]=w3; M1[4]=M2[4]=w4;
        M1[5]=w5; M2[5]=w5+C_dW5; M1[6]=w6; M2[6]=w6+C_dW6; M1[7]=w7; M2[7]=w7+C_dW7; M1[8]=w8; M2[8]=w8+C_dW8;
    }
    M1[9]  = C_W1_9;  M2[9]  = C_W1_9  + C_dW9;
    M1[10] = C_W1_10; M2[10] = C_W1_10;
    M1[11] = C_W1_11; M2[11] = C_W1_11;
    M1[12] = C_W1_12; M2[12] = C_W1_12;
    M1[13] = C_W1_13; M2[13] = C_W1_13 + C_dW13;
    M1[14] = d_w14[m]; M2[14] = d_w14[m] + C_dW14;
    uint64_t w15 = a15 - d_w15_sub[m];
    M1[15] = w15; M2[15] = w15;

    if (!is_sentinel && d_sample18 != nullptr && atomicCAS(&d_sample18->found, 0U, 1U) == 0U) {
        d_sample18->is_sentinel = 0;
        d_sample18->m0_seed = d_BCV.seed;
        for (int i=0;i<8;i++) d_sample18->cv[i] = cv[i];
        for (int i=0;i<16;i++) { d_sample18->M1[i] = M1[i]; d_sample18->M2[i] = M2[i]; }
    }

    uint64_t W1[36], W2[36];
    for (int i = 0; i < 16; i++) { W1[i] = M1[i]; W2[i] = M2[i]; }
    for (int i = 16; i < 36; i++) {
        W1[i] = s1(W1[i-2]) + W1[i-7] + s0(W1[i-15]) + W1[i-16];
        W2[i] = s1(W2[i-2]) + W2[i-7] + s0(W2[i-15]) + W2[i-16];
    }
    uint64_t a1=cv[0],b1=cv[1],c1=cv[2],d1=cv[3],e1=cv[4],f1=cv[5],g1=cv[6],h1=cv[7];
    uint64_t a2=cv[0],b2=cv[1],c2=cv[2],d2=cv[3],e2=cv[4],f2=cv[5],g2=cv[6],h2=cv[7];
    for (int t = 0; t < 36; t++) {
        uint64_t t1_1 = h1 + S1(e1) + Ch(e1, f1, g1) + d_K512[t] + W1[t];
        uint64_t t2_1 = S0(a1) + Maj(a1, b1, c1);
        h1=g1; g1=f1; f1=e1; e1=d1+t1_1; d1=c1; c1=b1; b1=a1; a1=t1_1+t2_1;

        uint64_t t1_2 = h2 + S1(e2) + Ch(e2, f2, g2) + d_K512[t] + W2[t];
        uint64_t t2_2 = S0(a2) + Maj(a2, b2, c2);
        h2=g2; g2=f2; f2=e2; e2=d2+t1_2; d2=c2; c2=b2; b2=a2; a2=t1_2+t2_2;

        if (t == 18) {
            if (a1 != a2 || e1 != e2) return;
        } else if (t == 19) {
            if (a1 != a2 || (e2 - e1) != 0xf800000000000000ULL) return;
            atomicAdd(&d_stats[5], 1ULL);
        } else if (t == 20) {
            if (a1 != a2 || e1 != e2) return;
            atomicAdd(&d_stats[6], 1ULL);
        } else if (t == 21) {
            if (a1 != a2 || e1 != e2) return;
            atomicAdd(&d_stats[7], 1ULL);
        } else if (t == 22) {
            if (a1 != a2 || e1 != e2) return;
            atomicAdd(&d_stats[8], 1ULL);
            if ((W2[23] - W1[23]) != 0x0800000000000000ULL) return;
            atomicAdd(&d_stats[9], 1ULL);
            if (W1[25] != W2[25]) return;
            atomicAdd(&d_stats[10], 1ULL);
        }
    }
    uint64_t H1[8] = {cv[0]+a1, cv[1]+b1, cv[2]+c1, cv[3]+d1, cv[4]+e1, cv[5]+f1, cv[6]+g1, cv[7]+h1};
    uint64_t H2[8] = {cv[0]+a2, cv[1]+b2, cv[2]+c2, cv[3]+d2, cv[4]+e2, cv[5]+f2, cv[6]+g2, cv[7]+h2};
    bool eq = true;
    for (int i = 0; i < 8; i++) if (H1[i] != H2[i]) eq = false;
    if (eq) {
        atomicAdd(&d_stats[11], 1ULL);
        if (atomicCAS(&d_col->found, 0U, 1U) == 0U) {
            d_col->is_sentinel = is_sentinel ? 1U : 0U;
            d_col->m0_seed = is_sentinel ? 0ULL : d_BCV.seed;
            for (int i = 0; i < 8; i++) { d_col->cv[i] = cv[i]; d_col->H1[i] = H1[i]; d_col->H2[i] = H2[i]; }
            for (int i = 0; i < 16; i++) { d_col->M1[i] = M1[i]; d_col->M2[i] = M2[i]; }
        }
    }
}

// Kernel 3 (Phase 2): 100% load-balanced chunked shared-memory tiled MitM!
// 4 warps/block, 1 warp per WorkChunk (<=64 U entries x <=64 L entries in shared mem!)
__global__ void __launch_bounds__(128) phase2_chunked_kernel(
    uint32_t n_chunks,
    const WorkChunk * __restrict__ d_chunks,
    const uint64_t  * __restrict__ d_L_lo,
    const uint64_t  * __restrict__ d_L_w17,
    const uint32_t  * __restrict__ d_L_e4idx,
    const uint32_t  * __restrict__ d_L_cnt,
    const uint64_t  * __restrict__ d_U_hie16,
    const uint64_t  * __restrict__ d_U_hia16,
    const uint64_t  * __restrict__ d_U_hiw17,
    const uint64_t  * __restrict__ d_U_a15,
    const uint8_t   * __restrict__ d_U_m,
    uint32_t sentinel_aid,
    unsigned long long * __restrict__ d_stats,
    CollisionRec * __restrict__ d_col,
    CollisionRec * __restrict__ d_sample18)
{
    __shared__ uint64_t sh_lo[4][ACT_CAP]; // 4 x 64 x 8B = 2 KB/block!
    const int warp_id = threadIdx.x >> 5;
    const int lane    = threadIdx.x & 31;
    const uint32_t cid = blockIdx.x * 4 + warp_id;
    if (cid >= n_chunks) return;

    const WorkChunk ck = d_chunks[cid];
    uint32_t nL = d_L_cnt[ck.act_id];
    if (nL == 0) return;
    if (nL > ACT_CAP) nL = ACT_CAP;

    const uint32_t l_base = ck.act_id << ACT_SHIFT;
    if (lane < (int)nL)      sh_lo[warp_id][lane]      = d_L_lo[l_base + lane];
    if (lane + 32 < (int)nL) sh_lo[warp_id][lane + 32] = d_L_lo[l_base + lane + 32];
    __syncwarp();

    unsigned long long c_pairs = 0, c_dir = 0, c_t1 = 0, c_s17 = 0, c_s18 = 0;
    const uint64_t E17_CONST = C_E1_13 + d_K512[17];

    for (uint32_t u = ck.u_beg + lane; u < ck.u_end; u += 32) {
        uint64_t hi_e16 = d_U_hie16[u];
        c_pairs += nL;
        for (uint32_t j = 0; j < nL; j++) {
            uint64_t lo = sh_lo[warp_id][j];
            uint64_t e16 = lo + hi_e16;
            if ((e16 & E16_DIRECT_MASK) != E16_DIRECT_VAL) continue;
            c_dir++;
            if ((((e16 >> 8) ^ (e16 >> 49)) & 1ULL) | (((e16 >> 55) ^ (e16 >> 59)) & 1ULL)) continue;
            c_t1++;
            uint64_t a16 = lo + d_U_hia16[u];
            if ((a16 >> 59) & 1ULL) continue;
            c_s17++;

            // Step 18 evaluation (only 5.93% of pairs!)
            uint8_t  m   = d_U_m[u];
            uint64_t a15 = d_U_a15[u];
            uint64_t e15 = a15 + d_dae15[m];
            uint64_t e14 = d_e14[m];
            uint64_t w17 = d_L_w17[l_base + j] + d_U_hiw17[u];
            uint64_t t1_17 = E17_CONST + S1(e16) + Ch(e16, e15, e14) + w17;
            uint64_t e17 = C_A1_13 + t1_17;
            uint64_t e17p = e17 + C_dE17;
            uint64_t e16p = e16 + C_dE16;
            uint64_t dt1_18 = C_dE14 + (S1(e17p) - S1(e17)) + (Ch(e17p, e16p, e15) - Ch(e17, e16, e15));
            if (dt1_18 != 0ULL) continue;

            uint64_t a17 = t1_17 + S0(a16) + Maj(a16, a15, d_a14[m]);
            if ((a17 >> 59) & 1ULL) continue;
            c_s18++;

            bool is_sent = (ck.act_id == sentinel_aid && j == 0);
            uint32_t e4_idx = d_L_e4idx[l_base + j];
            verify_and_record_step18_pass(e4_idx, m, a15, is_sent, d_stats, d_col, d_sample18);
        }
    }

    #pragma unroll
    for (int off = 16; off > 0; off >>= 1) {
        c_pairs += __shfl_down_sync(0xffffffffU, c_pairs, off);
        c_dir   += __shfl_down_sync(0xffffffffU, c_dir,   off);
        c_t1    += __shfl_down_sync(0xffffffffU, c_t1,    off);
        c_s17   += __shfl_down_sync(0xffffffffU, c_s17,   off);
        c_s18   += __shfl_down_sync(0xffffffffU, c_s18,   off);
    }
    if (lane == 0) {
        if (c_pairs) atomicAdd(&d_stats[0], c_pairs);
        if (c_dir)   atomicAdd(&d_stats[1], c_dir);
        if (c_t1)    atomicAdd(&d_stats[2], c_t1);
        if (c_s17)   atomicAdd(&d_stats[3], c_s17);
        if (c_s18)   atomicAdd(&d_stats[4], c_s18);
    }
}

struct UpperHost {
    std::vector<uint16_t>  bkt_map; // NUM_BKTS -> active_id or 0xFFFF
    std::vector<WorkChunk> chunks;
    std::vector<uint64_t>  hie16, hia16, hiw17, a15;
    std::vector<uint8_t>   m_arr;
    uint32_t n_act = 0;
    uint32_t sentinel_aid = UINT32_MAX;
    uint32_t sentinel_u   = UINT32_MAX;
    uint64_t h_a14[16], h_e14[16], h_w14[16], h_dae15[16], h_w15_sub[16];
};

static UpperHost build_upper_2_25() {
    auto t0 = std::chrono::steady_clock::now();
    UpperHost U;
    const int a14_bits[4] = {28, 40, 62, 63};

    uint64_t e15_fix_mask = (1ULL<<18) | (1ULL<<41) | (1ULL<<45);
    for (int b = 0; b < 64; b++) {
        int c0 = (int)((Ch(0ULL, C_E2_14, C_E2_13) >> b) & 1ULL) - (int)((Ch(0ULL, C_E1_14, C_E1_13) >> b) & 1ULL);
        int c1 = (int)((Ch(~0ULL, C_E2_14, C_E2_13) >> b) & 1ULL) - (int)((Ch(~0ULL, C_E1_14, C_E1_13) >> b) & 1ULL);
        if (c0 != c1) e15_fix_mask |= (1ULL << b);
    }
    const int tgt_p0 = (int)(((C_A1_15 >> 1) ^ (C_A1_15 >> 6)) & 1ULL);
    const int tgt_p1 = (int)(((C_A1_15 >> 0) ^ (C_A1_15 >> 53)) & 1ULL);
    const int tgt_p2 = (int)(((C_A1_15 >> 48) ^ (C_A1_15 >> 54)) & 1ULL);
    const int tgt_p  = (tgt_p0) | (tgt_p1 << 1) | (tgt_p2 << 2);

    for (int m = 0; m < 16; m++) {
        uint64_t a14 = C_A1_14;
        for (int i = 0; i < 4; i++) if ((m >> i) & 1) a14 ^= (1ULL << a14_bits[i]);
        uint64_t delta = a14 - C_A1_14;
        uint64_t e14 = C_E1_14 + delta;
        uint64_t w14 = C_W1_14 + delta;
        uint64_t t2_15 = S0(a14) + Maj(a14, C_A1_13, C_A1_12);
        U.h_a14[m] = a14; U.h_e14[m] = e14; U.h_w14[m] = w14;
        U.h_dae15[m] = C_A1_11 - t2_15;
        U.h_w15_sub[m] = t2_15 + C_E1_11 + S1(e14) + Ch(e14, C_E1_13, C_E1_12) + H_K512[15];
    }

    std::vector<std::vector<uint64_t>> a15_per_m(16);
    #pragma omp parallel for schedule(static)
    for (int m = 0; m < 16; m++) {
        uint64_t a14 = U.h_a14[m], dae15 = U.h_dae15[m];
        uint64_t am = (1ULL << 59);
        for (int b = 0; b < 64; b++) {
            if (b == 59) continue;
            int m0 = (int)((Maj(0ULL, a14, C_A2_13) >> b) & 1ULL) - (int)((Maj(0ULL, a14, C_A1_13) >> b) & 1ULL);
            int m1 = (int)((Maj(~0ULL, a14, C_A2_13) >> b) & 1ULL) - (int)((Maj(~0ULL, a14, C_A1_13) >> b) & 1ULL);
            if (m0 != m1) am |= (1ULL << b);
        }
        uint8_t reach[65][16];
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
        std::vector<uint64_t> sols;
        sols.reserve(1U << 21);
        struct Frame { int b; int st; uint64_t val; };
        std::vector<Frame> stk; stk.reserve(128);
        stk.push_back({0, 0, 0ULL});
        while (!stk.empty()) {
            Frame fr = stk.back(); stk.pop_back();
            if (fr.b == 64) { sols.push_back(fr.val); continue; }
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
                if (reach[b + 1][nst]) stk.push_back({b + 1, nst, fr.val | ((uint64_t)ab << b)});
            }
        }
        a15_per_m[m] = std::move(sols);
    }

    std::vector<uint32_t> bkt_cnt(NUM_BKTS, 0);
    std::vector<uint32_t> key_tmp(N_UPPER);
    #pragma omp parallel for schedule(static)
    for (int m = 0; m < 16; m++) {
        uint64_t a14 = U.h_a14[m], e14 = U.h_e14[m], w14 = U.h_w14[m], dae15 = U.h_dae15[m];
        uint64_t base_e16 = C_A1_12 + C_E1_12 + H_K512[16] + s1(w14) + C_W1_9;
        size_t off_m = (size_t)m << 21;
        for (size_t i = 0; i < (1U << 21); i++) {
            uint64_t a15 = a15_per_m[m][i];
            uint64_t e15 = a15 + dae15;
            uint64_t hie16 = base_e16 + S1(e15) + Ch(e15, e14, C_E1_13);
            uint64_t hia16 = hie16 - C_A1_12 + S0(a15) + Maj(a15, a14, C_A1_13);
            key_tmp[off_m + i] = key19_hi(hie16, hia16);
        }
    }
    for (uint32_t i = 0; i < N_UPPER; i++) bkt_cnt[key_tmp[i]]++;

    // Build compact active-bucket map (only non-empty buckets!) + WorkChunks (<=64 U each!)
    U.bkt_map.assign(NUM_BKTS, 0xFFFFU);
    std::vector<uint32_t> act_off;
    act_off.reserve(65536);
    act_off.push_back(0);
    for (uint32_t b = 0; b < NUM_BKTS; b++) {
        if (bkt_cnt[b] == 0) continue;
        if (U.n_act >= 65535U) { fprintf(stderr, "FATAL: n_act exceeds uint16_t!\n"); exit(1); }
        uint32_t aid = U.n_act++;
        U.bkt_map[b] = (uint16_t)aid;
        uint32_t beg = act_off.back();
        uint32_t end = beg + bkt_cnt[b];
        act_off.push_back(end);
        for (uint32_t u = beg; u < end; u += CHUNK_U) {
            U.chunks.push_back({aid, u, std::min(u + CHUNK_U, end)});
        }
    }

    U.hie16.resize(N_UPPER); U.hia16.resize(N_UPPER);
    U.hiw17.resize(N_UPPER); U.a15.resize(N_UPPER); U.m_arr.resize(N_UPPER);
    std::vector<uint32_t> cur = act_off;
    for (int m = 0; m < 16; m++) {
        uint64_t a14 = U.h_a14[m], e14 = U.h_e14[m], w14 = U.h_w14[m], dae15 = U.h_dae15[m], w15s = U.h_w15_sub[m];
        uint64_t base_e16 = C_A1_12 + C_E1_12 + H_K512[16] + s1(w14) + C_W1_9;
        size_t off_m = (size_t)m << 21;
        for (size_t i = 0; i < (1U << 21); i++) {
            uint64_t a15 = a15_per_m[m][i];
            uint64_t e15 = a15 + dae15;
            uint64_t hie16 = base_e16 + S1(e15) + Ch(e15, e14, C_E1_13);
            uint64_t hia16 = hie16 - C_A1_12 + S0(a15) + Maj(a15, a14, C_A1_13);
            uint64_t w15 = a15 - w15s;
            uint16_t aid = U.bkt_map[key_tmp[off_m + i]];
            uint32_t pos = cur[aid]++;
            U.hie16[pos] = hie16; U.hia16[pos] = hia16;
            U.hiw17[pos] = s1(w15) + C_W1_10;
            U.a15[pos]   = a15;   U.m_arr[pos] = (uint8_t)m;
            if (m == 0 && a15 == C_A1_15) {
                U.sentinel_aid = aid; U.sentinel_u = pos;
                if (w15 != C_W1_15_REF) { fprintf(stderr, "FATAL: sentinel w15 mismatch!\n"); exit(1); }
            }
        }
    }
    double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    printf("[INIT] Built N_U=%u (2^25.000) in %.2fs | active_bkts=%u/%u (%.1f%%) | L2_hot=%.1f MB | chunks=%zu (<=%u U/warp) | sentinel=(aid=%u, u=%u)\n",
           N_UPPER, dt, U.n_act, NUM_BKTS, 100.0 * U.n_act / NUM_BKTS,
           (double)U.n_act * ACT_CAP * 8.0 / (1024.0 * 1024.0), U.chunks.size(), CHUNK_U,
           U.sentinel_aid, U.sentinel_u);
    if (U.sentinel_aid == UINT32_MAX) { fprintf(stderr, "FATAL: Table 14 sentinel missing!\n"); exit(1); }
    return U;
}

static BatchCV make_batch_cv(const GoldenM0 &gm) {
    BatchCV B;
    B.seed = gm.seed;
    for (int i = 0; i < 8; i++) B.cv[i] = gm.cv[i];
    uint64_t Am1=gm.cv[0], Am2=gm.cv[1], Am3=gm.cv[2], Am4=gm.cv[3];
    uint64_t Em1=gm.cv[4], Em2=gm.cv[5], Em3=gm.cv[6], Em4=gm.cv[7];
    B.C3_M0 = Am1 + C_A1_3 - S0(C_A1_2);
    B.C2_M0 = Am2 + C_A1_2 - S0(C_A1_1);
    B.C1_M0 = Am3 + C_A1_1;
    B.C0_M0 = Am4 - S0(Am1) - Maj(Am1, Am2, Am3);
    B.W0_M0 = 0ULL - S0(Am1) - Maj(Am1, Am2, Am3) - Em4 - S1(Em1) - Ch(Em1, Em2, Em3) - H_K512[0];
    B.W1_M0 = C_A1_1 - Em3 - H_K512[1];
    B.W2_M0 = C_A1_2 - S0(C_A1_1) - Em2 - H_K512[2];
    B.A2A1_and   = C_A1_2 & C_A1_1; B.A2A1_xor   = C_A1_2 ^ C_A1_1;
    B.A1Am1_and  = C_A1_1 & Am1;    B.A1Am1_xor  = C_A1_1 ^ Am1;
    B.Am1Am2_and = Am1    & Am2;    B.Am1Am2_xor = Am1    ^ Am2;
    return B;
}

static void print_collision(const CollisionRec &rec, const char *label) {
    printf("\n==================== %s ====================\n", label);
    printf("is_sentinel=%u  m0_seed=%llu\n", rec.is_sentinel, (unsigned long long)rec.m0_seed);
    printf("CV :"); for (int i=0;i<8;i++)  printf(" %016llx", (unsigned long long)rec.cv[i]); printf("\n");
    printf("M1 :"); for (int i=0;i<16;i++) printf(" %016llx", (unsigned long long)rec.M1[i]); printf("\n");
    printf("M2 :"); for (int i=0;i<16;i++) printf(" %016llx", (unsigned long long)rec.M2[i]); printf("\n");
    printf("H1 :"); for (int i=0;i<8;i++)  printf(" %016llx", (unsigned long long)rec.H1[i]); printf("\n");
    printf("H2 :"); for (int i=0;i<8;i++)  printf(" %016llx", (unsigned long long)rec.H2[i]); printf("\n");
    printf("============================================================\n");
    fflush(stdout);
}

int main(int argc, char **argv) {
    double bench_sec = 30.0;
    bool search_mode = false;
    uint64_t seed_start = 0;
    bool cpu_init_only = false;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--bench") && i + 1 < argc) bench_sec = atof(argv[++i]);
        else if (!strcmp(argv[i], "--search")) search_mode = true;
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc) seed_start = strtoull(argv[++i], nullptr, 10);
        else if (!strcmp(argv[i], "--cpu-init-only")) cpu_init_only = true;
    }

    UpperHost U = build_upper_2_25();
    if (cpu_init_only) {
        printf("[CPU-INIT-ONLY] Done cleanly without touching GPU.\n");
        return 0;
    }

    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    printf("[GPU] %s | SMs=%d | L2=%.0f MB | VRAM=%.1f GB\n",
           prop.name, prop.multiProcessorCount, prop.l2CacheSize/(1024.0*1024.0), prop.totalGlobalMem/(1024.0*1024.0*1024.0));

    CUDA_CHECK(cudaMemcpyToSymbol(d_a14,     U.h_a14,     sizeof(U.h_a14)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_e14,     U.h_e14,     sizeof(U.h_e14)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_w14,     U.h_w14,     sizeof(U.h_w14)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_dae15,   U.h_dae15,   sizeof(U.h_dae15)));
    CUDA_CHECK(cudaMemcpyToSymbol(d_w15_sub, U.h_w15_sub, sizeof(U.h_w15_sub)));

    uint16_t  *d_bkt_map; WorkChunk *d_chunks;
    uint64_t  *d_U_hie16, *d_U_hia16, *d_U_hiw17, *d_U_a15; uint8_t *d_U_m;
    CUDA_CHECK(cudaMalloc(&d_bkt_map, NUM_BKTS * sizeof(uint16_t)));
    CUDA_CHECK(cudaMalloc(&d_chunks,  U.chunks.size() * sizeof(WorkChunk)));
    CUDA_CHECK(cudaMalloc(&d_U_hie16, N_UPPER * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&d_U_hia16, N_UPPER * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&d_U_hiw17, N_UPPER * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&d_U_a15,   N_UPPER * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&d_U_m,     N_UPPER * sizeof(uint8_t)));

    CUDA_CHECK(cudaMemcpy(d_bkt_map, U.bkt_map.data(), NUM_BKTS * sizeof(uint16_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_chunks,  U.chunks.data(),  U.chunks.size() * sizeof(WorkChunk), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_U_hie16, U.hie16.data(),   N_UPPER * sizeof(uint64_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_U_hia16, U.hia16.data(),   N_UPPER * sizeof(uint64_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_U_hiw17, U.hiw17.data(),   N_UPPER * sizeof(uint64_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_U_a15,   U.a15.data(),     N_UPPER * sizeof(uint64_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_U_m,     U.m_arr.data(),   N_UPPER * sizeof(uint8_t),  cudaMemcpyHostToDevice));

    size_t l_slots = (size_t)U.n_act * ACT_CAP;
    uint64_t *d_L_lo, *d_L_w17; uint32_t *d_L_e4idx, *d_L_cnt;
    CUDA_CHECK(cudaMalloc(&d_L_lo,    l_slots * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&d_L_w17,   l_slots * sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&d_L_e4idx, l_slots * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_L_cnt,   U.n_act * sizeof(uint32_t)));

    unsigned long long *d_stats;
    CollisionRec *d_col, *d_sample18;
    CUDA_CHECK(cudaMalloc(&d_stats,    N_STATS * sizeof(unsigned long long)));
    CUDA_CHECK(cudaMalloc(&d_col,      sizeof(CollisionRec)));
    CUDA_CHECK(cudaMalloc(&d_sample18, sizeof(CollisionRec)));

    uint32_t n_chunks = (uint32_t)U.chunks.size();
    uint32_t p2_blocks = (n_chunks + 3U) / 4U;

    // =========================================================================
    // TEST 1: SENTINEL SELF-TEST (Table 14 injected -> verify GPU captures exact 36-step collision!)
    // =========================================================================
    {
        CUDA_CHECK(cudaMemset(d_L_cnt, 0, U.n_act * sizeof(uint32_t)));
        CUDA_CHECK(cudaMemset(d_stats, 0, N_STATS * sizeof(unsigned long long)));
        CUDA_CHECK(cudaMemset(d_col,   0, sizeof(CollisionRec)));
        uint64_t lo_t14    = s0(H_T14_M1_0_8[1]) + H_T14_M1_0_8[0];
        uint64_t low17_t14 = s0(H_T14_M1_0_8[2]) + H_T14_M1_0_8[1];
        uint32_t one = 1U;
        uint32_t slot = U.sentinel_aid << ACT_SHIFT;
        CUDA_CHECK(cudaMemcpy(d_L_lo  + slot, &lo_t14,    sizeof(uint64_t), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_L_w17 + slot, &low17_t14, sizeof(uint64_t), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_L_cnt + U.sentinel_aid, &one, sizeof(uint32_t), cudaMemcpyHostToDevice));

        phase2_chunked_kernel<<<p2_blocks, 128>>>(
            n_chunks, d_chunks, d_L_lo, d_L_w17, d_L_e4idx, d_L_cnt,
            d_U_hie16, d_U_hia16, d_U_hiw17, d_U_a15, d_U_m,
            U.sentinel_aid, d_stats, d_col, nullptr);
        CUDA_CHECK(cudaDeviceSynchronize());

        CollisionRec h_col;
        CUDA_CHECK(cudaMemcpy(&h_col, d_col, sizeof(CollisionRec), cudaMemcpyDeviceToHost));
        bool hash_ok = (h_col.found == 1);
        for (int i = 0; i < 8; i++) if (h_col.H1[i] != T14_EXP_H[i] || h_col.H2[i] != T14_EXP_H[i]) hash_ok = false;
        printf("[SENTINEL SELF-TEST] %s (found=%u, H1[0]=%016llx vs exp=%016llx)\n",
               hash_ok ? "PASS 100%!!" : "FAIL!!", h_col.found,
               (unsigned long long)h_col.H1[0], (unsigned long long)T14_EXP_H[0]);
        if (!hash_ok) { print_collision(h_col, "SENTINEL FAIL"); return 1; }
    }

    // =========================================================================
    // TEST 2 / PRODUCTION: LIVE STANDARD-IV 36-STEP SHA-512 SEARCH & BENCHMARK
    // =========================================================================
    constexpr uint32_t GM_POOL = 256;
    GoldenM0 *d_gm; uint32_t *d_gm_cnt;
    CUDA_CHECK(cudaMalloc(&d_gm,     GM_POOL * sizeof(GoldenM0)));
    CUDA_CHECK(cudaMalloc(&d_gm_cnt, sizeof(uint32_t)));
    std::vector<GoldenM0> h_gm(GM_POOL);

    CUDA_CHECK(cudaMemset(d_stats,    0, N_STATS * sizeof(unsigned long long)));
    CUDA_CHECK(cudaMemset(d_col,      0, sizeof(CollisionRec)));
    CUDA_CHECK(cudaMemset(d_sample18, 0, sizeof(CollisionRec)));

    uint64_t next_seed = seed_start;
    uint64_t total_batches = 0;
    double t_gm_total = 0.0, t_p1_total = 0.0, t_p2_total = 0.0;
    auto t_wall0 = std::chrono::steady_clock::now();

    cudaEvent_t ev0, ev1, ev2;
    CUDA_CHECK(cudaEventCreate(&ev0));
    CUDA_CHECK(cudaEventCreate(&ev1));
    CUDA_CHECK(cudaEventCreate(&ev2));

    bool done = false;
    while (!done) {
        auto tg0 = std::chrono::steady_clock::now();
        uint32_t gm_found = 0;
        while (gm_found == 0) {
            CUDA_CHECK(cudaMemset(d_gm_cnt, 0, sizeof(uint32_t)));
            uint64_t n_try = 1ULL << 23;
            find_golden_m0_kernel<<<4096, 256>>>(next_seed, n_try, d_gm, d_gm_cnt, GM_POOL);
            CUDA_CHECK(cudaDeviceSynchronize());
            CUDA_CHECK(cudaMemcpy(&gm_found, d_gm_cnt, sizeof(uint32_t), cudaMemcpyDeviceToHost));
            next_seed += n_try;
        }
        if (gm_found > GM_POOL) gm_found = GM_POOL;
        CUDA_CHECK(cudaMemcpy(h_gm.data(), d_gm, gm_found * sizeof(GoldenM0), cudaMemcpyDeviceToHost));
        t_gm_total += std::chrono::duration<double>(std::chrono::steady_clock::now() - tg0).count();

        // Each Golden M0 sweeps ALL 2^29 E4 indices in ONE clean batch (~1.88M active L -> 20.5 MB L2 table!)
        for (uint32_t gi = 0; gi < gm_found && !done; gi++) {
            BatchCV bcv = make_batch_cv(h_gm[gi]);
            CUDA_CHECK(cudaMemcpyToSymbol(d_BCV, &bcv, sizeof(BatchCV)));

            CUDA_CHECK(cudaEventRecord(ev0));
            CUDA_CHECK(cudaMemset(d_L_cnt, 0, U.n_act * sizeof(uint32_t)));
            phase1_fill_compact_kernel<<<4096, 256>>>(
                1U << 29, d_bkt_map, d_L_lo, d_L_w17, d_L_e4idx, d_L_cnt, d_stats);
            CUDA_CHECK(cudaEventRecord(ev1));
            phase2_chunked_kernel<<<p2_blocks, 128>>>(
                n_chunks, d_chunks, d_L_lo, d_L_w17, d_L_e4idx, d_L_cnt,
                d_U_hie16, d_U_hia16, d_U_hiw17, d_U_a15, d_U_m,
                UINT32_MAX, d_stats, d_col, d_sample18);
            CUDA_CHECK(cudaEventRecord(ev2));
            CUDA_CHECK(cudaEventSynchronize(ev2));

            float ms1 = 0, ms2 = 0;
            CUDA_CHECK(cudaEventElapsedTime(&ms1, ev0, ev1));
            CUDA_CHECK(cudaEventElapsedTime(&ms2, ev1, ev2));
            t_p1_total += ms1 * 1e-3;
            t_p2_total += ms2 * 1e-3;
            total_batches++;

            if ((total_batches & 15ULL) == 0ULL || !search_mode) {
                uint32_t col_found = 0;
                CUDA_CHECK(cudaMemcpy(&col_found, &d_col->found, sizeof(uint32_t), cudaMemcpyDeviceToHost));
                if (col_found) {
                    CollisionRec h_col;
                    CUDA_CHECK(cudaMemcpy(&h_col, d_col, sizeof(CollisionRec), cudaMemcpyDeviceToHost));
                    print_collision(h_col, "WORLD-FIRST 36-STEP STANDARD-IV SHA-512 COLLISION!!");
                    done = true;
                    break;
                }
            }
            double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_wall0).count();
            if (!search_mode && wall >= bench_sec) { done = true; break; }
            if (search_mode && (total_batches % 500 == 0)) {
                unsigned long long st[N_STATS];
                CUDA_CHECK(cudaMemcpy(st, d_stats, sizeof(st), cudaMemcpyDeviceToHost));
                double s18_rate = (double)st[4] / wall;
                double eta_hr = (s18_rate > 0) ? (pow(2.0, 31.18) / s18_rate / 3600.0) : 9999.0;
                printf("[PROG %.1fs | %llu GM0 | seed=%llu] pairs=%.2fG/s s17=%.2fG/s s18=%.1fk/s | tail(19/20/21/22/W21/W23)=%llu/%llu/%llu/%llu/%llu/%llu | ETA=%.2fh\n",
                       wall, (unsigned long long)total_batches, (unsigned long long)next_seed,
                       (double)st[0]/wall/1e9, (double)st[3]/wall/1e9, s18_rate/1e3,
                       st[5], st[6], st[7], st[8], st[9], st[10], eta_hr);
                fflush(stdout);
            }
        }
    }

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_wall0).count();
    unsigned long long st[N_STATS];
    CUDA_CHECK(cudaMemcpy(st, d_stats, sizeof(st), cudaMemcpyDeviceToHost));
    printf("\n=== LIVE STANDARD-IV 36-STEP SHA-512 BENCHMARK SUMMARY (%.2fs, %llu Golden M0s) ===\n",
           wall, (unsigned long long)total_batches);
    printf("  Timing breakdown : GM0=%.2f%% (%.3fms/b), Phase1=%.2f%% (%.2fms/b), Phase2=%.2f%% (%.2fms/b), Total=%.2fms/b\n",
           100.0*t_gm_total/wall, 1e3*t_gm_total/total_batches,
           100.0*t_p1_total/wall, 1e3*t_p1_total/total_batches,
           100.0*t_p2_total/wall, 1e3*t_p2_total/total_batches,
           1e3*wall/total_batches);
    printf("  Phase 1 compact  : valid_E4=%.2fM/b (%.2f%%), active_hit=%.2fM/b (%.2f%%), inserted=%.2fM/b (capture=%.1f%%)\n",
           (double)st[12]/total_batches/1e6, 100.0*(double)st[12]/(total_batches*(double)(1U<<29)),
           (double)st[13]/total_batches/1e6, 100.0*(double)st[13]/std::max(st[12],1ULL),
           (double)st[14]/total_batches/1e6, 100.0*(double)st[14]/std::max(st[13],1ULL));
    printf("  Phase 2 rates    : %.2f G pairs/s | %.2f G dir/s (%.2f%%) | %.2f G s17_t1/s (%.2f%%) | %.2f G full_s17/s (%.2f%%)\n",
           (double)st[0]/wall/1e9, (double)st[1]/wall/1e9, 100.0*st[1]/std::max(st[0],1ULL),
           (double)st[2]/wall/1e9, 100.0*st[2]/std::max(st[0],1ULL),
           (double)st[3]/wall/1e9, 100.0*st[3]/std::max(st[0],1ULL));
    double s18_rate = (double)st[4] / wall;
    printf("  Step 18 full     : %llu (%.2f k/s, P(s18|s17)=2^%.2f)\n",
           st[4], s18_rate/1e3, (st[4] && st[3]) ? log2((double)st[4]/st[3]) : -99.0);
    printf("  Tail progression : s19=%llu (2^%.2f) | s20=%llu (2^%.2f) | s21=%llu | s22=%llu | W21=%llu | W23=%llu | COL=%llu\n",
           st[5], (st[5]&&st[4])?log2((double)st[5]/st[4]):-99.0,
           st[6], (st[6]&&st[5])?log2((double)st[6]/st[5]):-99.0,
           st[7], st[8], st[9], st[10], st[11]);
    if (s18_rate > 0) {
        double eta_hr = pow(2.0, 31.18) / s18_rate / 3600.0;
        printf("  >>> MEASURED EXPECTED TIME TO 36-STEP STANDARD-IV COLLISION: %.2f HOURS (%.1f days) <<<\n",
               eta_hr, eta_hr / 24.0);
    }

    CollisionRec h_s18;
    CUDA_CHECK(cudaMemcpy(&h_s18, d_sample18, sizeof(CollisionRec), cudaMemcpyDeviceToHost));
    if (h_s18.found) {
        printf("\n[SAMPLE_STEP18_JSON] {\"seed\":%llu,\"cv\":[", (unsigned long long)h_s18.m0_seed);
        for (int i=0;i<8;i++) printf("%s\"0x%016llx\"", i?",":"", (unsigned long long)h_s18.cv[i]);
        printf("],\"M1\":[");
        for (int i=0;i<16;i++) printf("%s\"0x%016llx\"", i?",":"", (unsigned long long)h_s18.M1[i]);
        printf("],\"M2\":[");
        for (int i=0;i<16;i++) printf("%s\"0x%016llx\"", i?",":"", (unsigned long long)h_s18.M2[i]);
        printf("]}\n");
    }
    return 0;
}
