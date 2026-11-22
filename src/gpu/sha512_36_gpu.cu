// 36-step Standard-IV SHA-512 i=1 Stage-3 GPU Solver & Benchmark (RTX 5090 sm_120)
// Split-Array 23-bit Guarded L2-Cache Table (64 MB hot c16 in 96 MB L2, cold w17_base read only on s17==0!)
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>

__device__ __host__ static inline uint64_t rotr64(uint64_t x, int n) { return (x >> n) | (x << (64 - n)); }
__device__ __host__ static inline uint64_t S0(uint64_t x) { return rotr64(x, 28) ^ rotr64(x, 34) ^ rotr64(x, 39); }
__device__ __host__ static inline uint64_t S1(uint64_t x) { return rotr64(x, 14) ^ rotr64(x, 18) ^ rotr64(x, 41); }
__device__ __host__ static inline uint64_t s0(uint64_t x) { return rotr64(x, 1) ^ rotr64(x, 8) ^ (x >> 7); }
__device__ __host__ static inline uint64_t s1(uint64_t x) { return rotr64(x, 19) ^ rotr64(x, 61) ^ (x >> 6); }
__device__ __host__ static inline uint64_t Ch(uint64_t e, uint64_t f, uint64_t g) { return (e & f) ^ ((~e) & g); }
__device__ __host__ static inline uint64_t Maj(uint64_t a, uint64_t b, uint64_t c) { return (a & b) ^ (a & c) ^ (b & c); }

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

__constant__ uint64_t d_IV[8] = {
    0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL, 0x3c6ef372fe94f82bULL, 0xa54ff53a5f1d36f1ULL,
    0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL, 0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL,
};

struct CoreParams {
    int k_bit;
    uint64_t A0[17], A1[17], E0[17], E1[17];
    uint64_t dA13, dE13, dE14, dA15, dE16;
    uint64_t rhs15;
    uint64_t p16_base;
    uint64_t w15_base;
    uint64_t e15_base;
    uint64_t a15_base;
    uint64_t free_bits[32];
    int n_free;
    uint32_t tgt23;
};

__constant__ CoreParams d_C;

// 23-bit key: 20 cluster bits + 3 S1(E16) linear parity bits!
__device__ __host__ static inline uint32_t pack23(uint64_t x, int k_bit) {
    if (k_bit == 63) {
        uint32_t c20 = (uint32_t)(((x >> 6) & 7ULL) | (((x >> 20) & 7ULL) << 3) |
                                  (((x >> 29) & 127ULL) << 6) | (((x >> 43) & 127ULL) << 13));
        uint32_t p0 = (uint32_t)(((x >> 4) ^ (x >> 41)) & 1ULL);
        uint32_t p1 = (uint32_t)(((x >> 12) ^ (x >> 53)) & 1ULL);
        uint32_t p2 = (uint32_t)(((x >> 59) ^ (x >> 63)) & 1ULL);
        return c20 | (p0 << 20) | (p1 << 21) | (p2 << 22);
    } else {
        uint32_t c20 = (uint32_t)(((x >> 5) & 7ULL) | (((x >> 14) & 127ULL) << 3) |
                                  (((x >> 28) & 127ULL) << 10) | (((x >> 55) & 7ULL) << 17));
        // k=48 S1(E16) active at {12, 20, 30}: b=12->{26,30,53}, b=20->{34,38,61}, b=30->{44,48,7}
        uint32_t p0 = (uint32_t)(((x >> 26) ^ (x >> 53)) & 1ULL);
        uint32_t p1 = (uint32_t)(((x >> 38) ^ (x >> 61)) & 1ULL);
        uint32_t p2 = (uint32_t)(((x >> 44) ^ (x >> 48)) & 1ULL);
        return c20 | (p0 << 20) | (p1 << 21) | (p2 << 22);
    }
}

__device__ __host__ static inline uint32_t sub23(uint32_t tgt, uint32_t c, int k_bit) {
    uint32_t p = (tgt ^ c) & (7u << 20);
    if (k_bit == 63) {
        uint32_t k0 = (tgt - c) & 7u;
        uint32_t k1 = ((tgt >> 3) - (c >> 3)) & 7u;
        uint32_t k2 = ((tgt >> 6) - (c >> 6)) & 127u;
        uint32_t k3 = ((tgt >> 13) - (c >> 13)) & 127u;
        return k0 | (k1 << 3) | (k2 << 6) | (k3 << 13) | p;
    } else {
        uint32_t k0 = (tgt - c) & 7u;
        uint32_t k1 = ((tgt >> 3) - (c >> 3)) & 127u;
        uint32_t k2 = ((tgt >> 10) - (c >> 10)) & 127u;
        uint32_t k3 = ((tgt >> 17) - (c >> 17)) & 7u;
        return k0 | (k1 << 3) | (k2 << 10) | (k3 << 17) | p;
    }
}

struct Hit { uint64_t e15, c16, w17_base; uint32_t bucket, depth; };

__global__ void search_pi_23_kernel(
    const uint64_t * __restrict__ tbl_c16,  // 64 MB HOT IN 96 MB L2 CACHE!
    const uint64_t * __restrict__ tbl_w17,  // 64 MB COLD (read only on s17==0!)
    uint64_t base_idx,
    uint64_t iters_per_thread,
    unsigned long long *stats,
    Hit *hits,
    uint32_t max_hits
) {
    uint64_t tid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    uint64_t start = base_idx + tid * iters_per_thread;

    const int k_bit = d_C.k_bit;
    const uint64_t e15_base  = d_C.e15_base;
    const uint64_t rhs15     = d_C.rhs15;
    const uint64_t p16_base  = d_C.p16_base;
    const uint64_t w15_base  = d_C.w15_base;
    const uint64_t e14_0     = d_C.E0[14], e14_1 = d_C.E1[14];
    const uint64_t e13_0     = d_C.E0[13];
    const uint64_t a14_0     = d_C.A0[14], a13_0 = d_C.A0[13], a12_0 = d_C.A0[12];
    const uint64_t dA13      = d_C.dA13,   dE13  = d_C.dE13,   dE14  = d_C.dE14, dE16 = d_C.dE16;
    const uint64_t e17_base  = a13_0 + e13_0 + d_K512[17];
    const uint64_t a16_const = (uint64_t)0 - a12_0;
    const uint64_t a14_k     = (a14_0 >> k_bit) & 1ULL;

    int nf = d_C.n_free;
    uint64_t g = start ^ (start >> 1);
    uint64_t e15 = e15_base;
    for (int b = 0; b < nf; b++) if ((g >> b) & 1ULL) e15 ^= d_C.free_bits[b];

    unsigned long long c_s17 = 0, c_s18 = 0;
    for (uint64_t it = 0; it < iters_per_thread; it++) {
        uint64_t p16 = p16_base + S1(e15) + Ch(e15, e14_0, e13_0);
        uint32_t key = pack23(p16, k_bit);
        uint64_t c16 = tbl_c16[key];
        uint64_t e16 = c16 + p16;
        uint64_t e16p = e16 + dE16;
        uint64_t s17 = dE13 + (S1(e16p) - S1(e16)) + (Ch(e16p, e15, e14_1) - Ch(e16, e15, e14_0));
        if (s17 == 0ULL) {
            c_s17++;
            uint64_t a15 = e15 - rhs15;
            uint64_t a16 = e16 + a16_const + S0(a15) + Maj(a15, a14_0, a13_0);
            if (((a16 >> k_bit) & 1ULL) == a14_k) {
                uint64_t w15 = e15 - w15_base;
                uint64_t w17 = tbl_w17[key] + s1(w15);
                uint64_t e17 = e17_base + S1(e16) + Ch(e16, e15, e14_0) + w17;
                uint64_t e17p = e17 + dA13;
                uint64_t t1_18 = dE14 + (S1(e17p) - S1(e17)) + (Ch(e17p, e16p, e15) - Ch(e17, e16, e15));
                if (t1_18 == 0ULL) {
                    uint64_t a17 = e17 - a13_0 + S0(a16) + Maj(a16, a15, a14_0);
                    if ((((a17 ^ a16) >> k_bit) & 1ULL) == 0ULL) {
                        c_s18++;
                        uint32_t slot = atomicAdd((unsigned int*)&stats[4], 1u);
                        if (slot < max_hits) hits[slot] = {e15, c16, tbl_w17[key], key, 18};
                    }
                }
            }
        }
        int bit = __ffsll((long long)(start + it + 1)) - 1;
        if (bit < nf) e15 ^= d_C.free_bits[bit];
    }
    atomicAdd(&stats[0], iters_per_thread);
    atomicAdd(&stats[1], c_s17);
    atomicAdd(&stats[2], c_s18);
}

__global__ void fill_bench23_kernel(uint64_t *tbl_c16, uint64_t *tbl_w17, uint32_t tgt23, int k_bit) {
    uint32_t key = blockIdx.x * blockDim.x + threadIdx.x;
    if (key >= (1u << 23)) return;
    uint32_t c23 = sub23(tgt23, key, k_bit);
    uint64_t r = (uint64_t)key * 6364136223846793005ULL + 1442695040888963407ULL;
    r ^= r >> 33; r *= 0xff51afd7ed558ccdULL; r ^= r >> 33;
    uint64_t c16 = r;
    if (k_bit == 63) {
        // 3-bit guard below each of the 4 clusters + 1-bit guard below {4,12,41,53,59,63}!
        uint64_t gmask = (7ULL<<3)|(7ULL<<17)|(7ULL<<26)|(7ULL<<40) | (1ULL<<11)|(1ULL<<52)|(1ULL<<58)|(1ULL<<62);
        uint64_t cmask = (7ULL<<6)|(7ULL<<20)|(127ULL<<29)|(127ULL<<43) | (1ULL<<4)|(1ULL<<41)|(1ULL<<12)|(1ULL<<53)|(1ULL<<59)|(1ULL<<63);
        c16 = r & ~(gmask | cmask);
        c16 |= ((uint64_t)(c23 & 7u) << 6) | ((uint64_t)((c23 >> 3) & 7u) << 20) |
               ((uint64_t)((c23 >> 6) & 127u) << 29) | ((uint64_t)((c23 >> 13) & 127u) << 43);
        // Set parity bits directly on {4, 12, 59} (with {41, 53, 63} = 0):
        if ((c23 >> 20) & 1u) c16 |= (1ULL << 4);
        if ((c23 >> 21) & 1u) c16 |= (1ULL << 12);
        if ((c23 >> 22) & 1u) c16 |= (1ULL << 59);
    }
    tbl_c16[key] = c16;
    tbl_w17[key] = r ^ 0x9e3779b97f4a7c15ULL;
}

int main(int argc, char **argv) {
    if (argc < 2) return 1;
    double run_secs = (argc >= 3) ? atof(argv[2]) : 5.0;
    CoreParams hC{};
    FILE *f = fopen(argv[1], "r");
    if (!f) return 1;
    if (fscanf(f, "%d", &hC.k_bit) != 1) return 1;
    for (int t = 1; t <= 16; t++)
        if (fscanf(f, "%llx %llx %llx %llx", (unsigned long long*)&hC.A0[t], (unsigned long long*)&hC.A1[t],
                   (unsigned long long*)&hC.E0[t], (unsigned long long*)&hC.E1[t]) != 4) return 1;
    int nf_raw = 0;
    if (fscanf(f, "%llx %llx %llx %llx %llx %d", (unsigned long long*)&hC.rhs15, (unsigned long long*)&hC.p16_base,
               (unsigned long long*)&hC.w15_base, (unsigned long long*)&hC.e15_base,
               (unsigned long long*)&hC.a15_base, &nf_raw) != 6) return 1;
    hC.n_free = (nf_raw > 32) ? 32 : nf_raw;
    for (int i = 0; i < nf_raw; i++) {
        unsigned long long v = 0;
        if (fscanf(f, "%llx", &v) != 1) return 1;
        if (i < 32) hC.free_bits[i] = v;
    }
    fclose(f);
    hC.dA13 = hC.A1[13] - hC.A0[13]; hC.dE13 = hC.E1[13] - hC.E0[13];
    hC.dE14 = hC.E1[14] - hC.E0[14]; hC.dA15 = hC.A1[15] - hC.A0[15]; hC.dE16 = hC.E1[16] - hC.E0[16];
    hC.tgt23 = pack23(hC.E0[16], hC.k_bit);

    cudaDeviceProp prop{}; cudaGetDeviceProperties(&prop, 0);
    cudaMemcpyToSymbol(d_C, &hC, sizeof(CoreParams));

    uint64_t *d_c16 = nullptr, *d_w17 = nullptr;
    size_t sz = (size_t)(1u << 23) * sizeof(uint64_t);
    cudaMalloc(&d_c16, sz); cudaMalloc(&d_w17, sz);
    fill_bench23_kernel<<<(1u << 23) / 256, 256>>>(d_c16, d_w17, hC.tgt23, hC.k_bit);
    cudaDeviceSynchronize();

    unsigned long long *d_stats, h_stats[5] = {0};
    Hit *d_hits;
    cudaMalloc(&d_stats, 5 * sizeof(unsigned long long));
    cudaMalloc(&d_hits, 4096 * sizeof(Hit));
    cudaMemset(d_stats, 0, 5 * sizeof(unsigned long long));

    int blocks = prop.multiProcessorCount * 16, threads = 256;
    uint64_t iters = 4096, batch = (uint64_t)blocks * threads * iters;
    auto t0 = std::chrono::steady_clock::now();
    uint64_t base_idx = 1;
    int launches = 0;
    while (true) {
        search_pi_23_kernel<<<blocks, threads>>>(d_c16, d_w17, base_idx, iters, d_stats, d_hits, 4096);
        cudaDeviceSynchronize();
        base_idx += batch; launches++;
        if (std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() >= run_secs) break;
    }
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    cudaMemcpy(h_stats, d_stats, 5 * sizeof(unsigned long long), cudaMemcpyDeviceToHost);

    double pi_rate  = (double)h_stats[0] / el;
    double s17_rate = (double)h_stats[1] / el;
    double s18_rate = (double)h_stats[2] / el;
    printf("[23-BIT SPLIT-L2 BENCH %.2fs, %d launches | 64 MB hot c16 in %.0f MB L2]\n", el, launches, prop.l2CacheSize/(1024.0*1024.0));
    printf("  pi tested   : %llu (%.2f G pi/s)\n", h_stats[0], pi_rate/1e9);
    printf("  s17==0 pass : %llu (%.2f G s17/s, hit/pi=%.3f = 2^%.2f)\n",
           h_stats[1], s17_rate/1e9, (double)h_stats[1]/h_stats[0], log2((double)h_stats[1]/h_stats[0]));
    printf("  Step18 pass : %llu (%.3f M s18/s, P(s18|s17)=2^%.2f)\n",
           h_stats[2], s18_rate/1e6, h_stats[2] ? log2((double)h_stats[2]/h_stats[1]) : -99.0);
    if (s18_rate > 0) {
        double eta_k63 = pow(2.0, 34.77) / s18_rate;
        double eta_k48 = pow(2.0, 33.60) / s18_rate;
        printf("  >>> PROJECTED 1x RTX 5090 TIME TO 36-STEP COLLISION: k63=%.1f min (%.2f hr) | k48=%.1f min (%.2f hr) <<<\n",
               eta_k63/60.0, eta_k63/3600.0, eta_k48/60.0, eta_k48/3600.0);
    }
    cudaFree(d_c16); cudaFree(d_w17); cudaFree(d_stats); cudaFree(d_hits);
    return 0;
}
