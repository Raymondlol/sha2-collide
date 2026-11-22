// Throughput of the r-step SHA-512 compression function on one GPU (r = 36 and r = 80).
// Each thread compresses many distinct runtime-varied blocks from a runtime CV; all 8 output
// words are folded into a checksum so nothing can be optimized away.
// Correctness self-test first: 36-step compress(IV, M0 of the collision) == recorded CV, and
// 80-step SHA-512("abc") == FIPS test vector.
// Build: nvcc -O3 -arch=sm_120 bench_compress.cu -o bench_compress
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <cuda_runtime.h>

#define CK(x) do { cudaError_t e = (x); if (e != cudaSuccess) { printf("CUDA %s @%d\n", cudaGetErrorString(e), __LINE__); return 1; } } while (0)

__constant__ uint64_t K[80] = {
0x428a2f98d728ae22ULL,0x7137449123ef65cdULL,0xb5c0fbcfec4d3b2fULL,0xe9b5dba58189dbbcULL,0x3956c25bf348b538ULL,0x59f111f1b605d019ULL,0x923f82a4af194f9bULL,0xab1c5ed5da6d8118ULL,
0xd807aa98a3030242ULL,0x12835b0145706fbeULL,0x243185be4ee4b28cULL,0x550c7dc3d5ffb4e2ULL,0x72be5d74f27b896fULL,0x80deb1fe3b1696b1ULL,0x9bdc06a725c71235ULL,0xc19bf174cf692694ULL,
0xe49b69c19ef14ad2ULL,0xefbe4786384f25e3ULL,0x0fc19dc68b8cd5b5ULL,0x240ca1cc77ac9c65ULL,0x2de92c6f592b0275ULL,0x4a7484aa6ea6e483ULL,0x5cb0a9dcbd41fbd4ULL,0x76f988da831153b5ULL,
0x983e5152ee66dfabULL,0xa831c66d2db43210ULL,0xb00327c898fb213fULL,0xbf597fc7beef0ee4ULL,0xc6e00bf33da88fc2ULL,0xd5a79147930aa725ULL,0x06ca6351e003826fULL,0x142929670a0e6e70ULL,
0x27b70a8546d22ffcULL,0x2e1b21385c26c926ULL,0x4d2c6dfc5ac42aedULL,0x53380d139d95b3dfULL,0x650a73548baf63deULL,0x766a0abb3c77b2a8ULL,0x81c2c92e47edaee6ULL,0x92722c851482353bULL,
0xa2bfe8a14cf10364ULL,0xa81a664bbc423001ULL,0xc24b8b70d0f89791ULL,0xc76c51a30654be30ULL,0xd192e819d6ef5218ULL,0xd69906245565a910ULL,0xf40e35855771202aULL,0x106aa07032bbd1b8ULL,
0x19a4c116b8d2d0c8ULL,0x1e376c085141ab53ULL,0x2748774cdf8eeb99ULL,0x34b0bcb5e19b48a8ULL,0x391c0cb3c5c95a63ULL,0x4ed8aa4ae3418acbULL,0x5b9cca4f7763e373ULL,0x682e6ff3d6b2b8a3ULL,
0x748f82ee5defb2fcULL,0x78a5636f43172f60ULL,0x84c87814a1f0ab72ULL,0x8cc702081a6439ecULL,0x90befffa23631e28ULL,0xa4506cebde82bde9ULL,0xbef9a3f7b2c67915ULL,0xc67178f2e372532bULL,
0xca273eceea26619cULL,0xd186b8c721c0c207ULL,0xeada7dd6cde0eb1eULL,0xf57d4f7fee6ed178ULL,0x06f067aa72176fbaULL,0x0a637dc5a2c898a6ULL,0x113f9804bef90daeULL,0x1b710b35131c471bULL,
0x28db77f523047d84ULL,0x32caab7b40c72493ULL,0x3c9ebe0a15c9bebcULL,0x431d67c49c100d4cULL,0x4cc5d4becb3e42b6ULL,0x597f299cfc657e2aULL,0x5fcb6fab3ad6faecULL,0x6c44198c4a475817ULL };

__device__ __forceinline__ uint64_t rotr(uint64_t x, int n) { return (x >> n) | (x << (64 - n)); }
#define S0(x) (rotr(x,28) ^ rotr(x,34) ^ rotr(x,39))
#define S1(x) (rotr(x,14) ^ rotr(x,18) ^ rotr(x,41))
#define s0(x) (rotr(x,1) ^ rotr(x,8) ^ ((x) >> 7))
#define s1(x) (rotr(x,19) ^ rotr(x,61) ^ ((x) >> 6))

template <int R>
__device__ __forceinline__ void compress(const uint64_t cv[8], uint64_t W[16], uint64_t out[8]) {
    uint64_t a=cv[0],b=cv[1],c=cv[2],d=cv[3],e=cv[4],f=cv[5],g=cv[6],h=cv[7];
    #pragma unroll
    for (int t = 0; t < R; t++) {
        uint64_t w;
        if (t < 16) w = W[t];
        else { w = s1(W[(t-2)&15]) + W[(t-7)&15] + s0(W[(t-15)&15]) + W[t&15]; W[t&15] = w; }
        uint64_t t1 = h + S1(e) + ((e & f) ^ (~e & g)) + K[t] + w;
        uint64_t t2 = S0(a) + ((a & b) ^ (a & c) ^ (b & c));
        h = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
    }
    out[0]=cv[0]+a; out[1]=cv[1]+b; out[2]=cv[2]+c; out[3]=cv[3]+d;
    out[4]=cv[4]+e; out[5]=cv[5]+f; out[6]=cv[6]+g; out[7]=cv[7]+h;
}

struct Params { uint64_t cv[8]; uint64_t base[16]; };

template <int R>
__global__ void bench(Params P, uint32_t iters, uint64_t *sink) {
    uint64_t tid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    uint64_t acc = 0;
    for (uint32_t it = 0; it < iters; it++) {
        uint64_t x = (tid << 32) ^ it;
        uint64_t W[16], out[8];
        #pragma unroll
        for (int i = 0; i < 16; i++) W[i] = P.base[i] + x;
        compress<R>(P.cv, W, out);
        acc ^= out[0] ^ out[1] ^ out[2] ^ out[3] ^ out[4] ^ out[5] ^ out[6] ^ out[7];
    }
    sink[tid] = acc;
}

template <int R>
__global__ void selftest(const uint64_t *cv, const uint64_t *M, uint64_t *out) {
    uint64_t W[16], c[8], o[8];
    for (int i = 0; i < 16; i++) W[i] = M[i];
    for (int i = 0; i < 8; i++) c[i] = cv[i];
    compress<R>(c, W, o);
    for (int i = 0; i < 8; i++) out[i] = o[i];
}

static const uint64_t IV[8] = {0x6a09e667f3bcc908ULL,0xbb67ae8584caa73bULL,0x3c6ef372fe94f82bULL,0xa54ff53a5f1d36f1ULL,
                               0x510e527fade682d1ULL,0x9b05688c2b3e6c1fULL,0x1f83d9abfb41bd6bULL,0x5be0cd19137e2179ULL};

template <int R>
static int run_selftest(const uint64_t *M, const uint64_t *expect, const char *name) {
    uint64_t *dcv, *dM, *dout, out[8];
    CK(cudaMalloc(&dcv, 64)); CK(cudaMalloc(&dM, 128)); CK(cudaMalloc(&dout, 64));
    CK(cudaMemcpy(dcv, IV, 64, cudaMemcpyHostToDevice)); CK(cudaMemcpy(dM, M, 128, cudaMemcpyHostToDevice));
    selftest<R><<<1,1>>>(dcv, dM, dout); CK(cudaDeviceSynchronize());
    CK(cudaMemcpy(out, dout, 64, cudaMemcpyDeviceToHost));
    int ok = 1; for (int i = 0; i < 8; i++) if (out[i] != expect[i]) ok = 0;
    printf("[selftest] %-40s %s\n", name, ok ? "OK" : "FAIL");
    cudaFree(dcv); cudaFree(dM); cudaFree(dout);
    return ok;
}

template <int R>
static double run_bench(int blocks, int threads, uint32_t iters, int reps) {
    Params P;
    for (int i = 0; i < 8; i++) P.cv[i] = IV[i] ^ 0x0123456789abcdefULL * (i + 1);
    for (int i = 0; i < 16; i++) P.base[i] = 0x9e3779b97f4a7c15ULL * (i + 7);
    uint64_t *sink; size_t n = (size_t)blocks * threads;
    cudaMalloc(&sink, n * 8);
    bench<R><<<blocks, threads>>>(P, 16, sink); cudaDeviceSynchronize();   // warm-up
    cudaEvent_t a, b; cudaEventCreate(&a); cudaEventCreate(&b);
    double best = 0;
    for (int r = 0; r < reps; r++) {
        cudaEventRecord(a);
        bench<R><<<blocks, threads>>>(P, iters, sink);
        cudaEventRecord(b); cudaEventSynchronize(b);
        float ms; cudaEventElapsedTime(&ms, a, b);
        double rate = (double)n * iters / (ms / 1e3);
        if (rate > best) best = rate;
        printf("  R=%d rep %d: %.1f ms, %.3e compressions/s\n", R, r, ms, rate);
    }
    cudaFree(sink);
    return best;
}

int main() {
    cudaDeviceProp p; cudaGetDeviceProperties(&p, 0);
    printf("[gpu] %s, %d SMs\n", p.name, p.multiProcessorCount);
    // 36-step compress(IV, M0 of the found collision) must equal the recorded CV
    uint64_t M0[16] = {0x000001e1b34abbd2ULL, 0x5348413531323336ULL, 0x434f4c4c49444521ULL};
    uint64_t CV[8] = {0xf1e42e9f77e4cca4ULL,0xa785f886ea02658dULL,0xcf52ba6f18e44e11ULL,0x644ba4c897d35d12ULL,
                      0xc4cd7917344f6760ULL,0x3f554360df46b76aULL,0x1a8281d8e746e746ULL,0x0ef673a2b4c69e41ULL};
    // 80-step SHA-512("abc") single padded block (FIPS 180-4 test vector)
    uint64_t ABC[16] = {0x6162638000000000ULL}; ABC[15] = 24;
    uint64_t ABCH[8] = {0xddaf35a193617abaULL,0xcc417349ae204131ULL,0x12e6fa4e89a97ea2ULL,0x0a9eeee64b55d39aULL,
                        0x2192992a274fc1a8ULL,0x36ba3c23a3feebbdULL,0x454d4423643ce80eULL,0x2a9ac94fa54ca49fULL};
    if (!run_selftest<36>(M0, CV, "36-step compress(IV, M0) == recorded CV")) return 1;
    if (!run_selftest<80>(ABC, ABCH, "80-step SHA-512(\"abc\") == FIPS vector")) return 1;

    int blocks = p.multiProcessorCount * 32, threads = 256;
    uint32_t iters = 4096;
    double r36 = run_bench<36>(blocks, threads, iters, 5);
    double r80 = run_bench<80>(blocks, threads, iters / 2, 5);
    printf("[result] 36-step: %.3e compressions/s (2^%.2f /s)\n", r36, log2(r36));
    printf("[result] 80-step: %.3e compressions/s (2^%.2f /s)  ratio 36/80 = %.2f\n", r80, log2(r80), r36 / r80);
    return 0;
}
