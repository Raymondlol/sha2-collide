"""Independent reference implementation of SHA-512 (full 80-step & arbitrary R-step)
with mandatory self-test against Python hashlib.sha512."""
import hashlib, struct

M64 = (1 << 64) - 1

IV512 = [
    0x6a09e667f3bcc908, 0xbb67ae8584caa73b, 0x3c6ef372fe94f82b, 0xa54ff53a5f1d36f1,
    0x510e527fade682d1, 0x9b05688c2b3e6c1f, 0x1f83d9abfb41bd6b, 0x5be0cd19137e2179,
]

K512 = [
    0x428a2f98d728ae22, 0x7137449123ef65cd, 0xb5c0fbcfec4d3b2f, 0xe9b5dba58189dbbc,
    0x3956c25bf348b538, 0x59f111f1b605d019, 0x923f82a4af194f9b, 0xab1c5ed5da6d8118,
    0xd807aa98a3030242, 0x12835b0145706fbe, 0x243185be4ee4b28c, 0x550c7dc3d5ffb4e2,
    0x72be5d74f27b896f, 0x80deb1fe3b1696b1, 0x9bdc06a725c71235, 0xc19bf174cf692694,
    0xe49b69c19ef14ad2, 0xefbe4786384f25e3, 0x0fc19dc68b8cd5b5, 0x240ca1cc77ac9c65,
    0x2de92c6f592b0275, 0x4a7484aa6ea6e483, 0x5cb0a9dcbd41fbd4, 0x76f988da831153b5,
    0x983e5152ee66dfab, 0xa831c66d2db43210, 0xb00327c898fb213f, 0xbf597fc7beef0ee4,
    0xc6e00bf33da88fc2, 0xd5a79147930aa725, 0x06ca6351e003826f, 0x142929670a0e6e70,
    0x27b70a8546d22ffc, 0x2e1b21385c26c926, 0x4d2c6dfc5ac42aed, 0x53380d139d95b3df,
    0x650a73548baf63de, 0x766a0abb3c77b2a8, 0x81c2c92e47edaee6, 0x92722c851482353b,
    0xa2bfe8a14cf10364, 0xa81a664bbc423001, 0xc24b8b70d0f89791, 0xc76c51a30654be30,
    0xd192e819d6ef5218, 0xd69906245565a910, 0xf40e35855771202a, 0x106aa07032bbd1b8,
    0x19a4c116b8d2d0c8, 0x1e376c085141ab53, 0x2748774cdf8eeb99, 0x34b0bcb5e19b48a8,
    0x391c0cb3c5c95a63, 0x4ed8aa4ae3418acb, 0x5b9cca4f7763e373, 0x682e6ff3d6b2b8a3,
    0x748f82ee5defb2fc, 0x78a5636f43172f60, 0x84c87814a1f0ab72, 0x8cc702081a6439ec,
    0x90befffa23631e28, 0xa4506cebde82bde9, 0xbef9a3f7b2c67915, 0xc67178f2e372532b,
    0xca273eceea26619c, 0xd186b8c721c0c207, 0xeada7dd6cde0eb1e, 0xf57d4f7fee6ed178,
    0x06f067aa72176fba, 0x0a637dc5a2c898a6, 0x113f9804bef90dae, 0x1b710b35131c471b,
    0x28db77f523047d84, 0x32caab7b40c72493, 0x3c9ebe0a15c9bebc, 0x431d67c49c100d4c,
    0x4cc5d4becb3e42b6, 0x597f299cfc657e2a, 0x5fcb6fab3ad6faec, 0x6c44198c4a475817,
]

def rotr(x, n): return ((x >> n) | (x << (64 - n))) & M64
def s0(x): return rotr(x, 1) ^ rotr(x, 8) ^ (x >> 7)
def s1(x): return rotr(x, 19) ^ rotr(x, 61) ^ (x >> 6)
def S0(x): return rotr(x, 28) ^ rotr(x, 34) ^ rotr(x, 39)
def S1(x): return rotr(x, 14) ^ rotr(x, 18) ^ rotr(x, 41)
def Ch(e, f, g): return (e & f) ^ ((~e & M64) & g)
def Maj(a, b, c): return (a & b) ^ (a & c) ^ (b & c)

def expand_w(W16, steps=80):
    W = list(W16)
    for t in range(16, steps):
        W.append((W[t-16] + s0(W[t-15]) + W[t-7] + s1(W[t-2])) & M64)
    return W

def compress_trace(H, W16, steps=80):
    """Returns (H_out, A, E, W) where A[t], E[t] are dicts indexed by t in -4..steps-1.
       Convention: A[-1]=H[0]..A[-4]=H[3], E[-1]=H[4]..E[-4]=H[7];
       Round t produces A[t], E[t] using W[t], K512[t]."""
    W = expand_w(W16, steps)
    A = {-1: H[0], -2: H[1], -3: H[2], -4: H[3]}
    E = {-1: H[4], -2: H[5], -3: H[6], -4: H[7]}
    for t in range(steps):
        T1 = (E[t-4] + S1(E[t-1]) + Ch(E[t-1], E[t-2], E[t-3]) + K512[t] + W[t]) & M64
        T2 = (S0(A[t-1]) + Maj(A[t-1], A[t-2], A[t-3])) & M64
        E[t] = (A[t-4] + T1) & M64
        A[t] = (T1 + T2) & M64
    H_out = [
        (H[0] + A[steps-1]) & M64, (H[1] + A[steps-2]) & M64,
        (H[2] + A[steps-3]) & M64, (H[3] + A[steps-4]) & M64,
        (H[4] + E[steps-1]) & M64, (H[5] + E[steps-2]) & M64,
        (H[6] + E[steps-3]) & M64, (H[7] + E[steps-4]) & M64,
    ]
    return H_out, A, E, W

def sha512_reduced(msg_bytes, steps=80):
    bit_len = len(msg_bytes) * 8
    padded = msg_bytes + b'\x80'
    while (len(padded) % 128) != 112:
        padded += b'\x00'
    padded += struct.pack('>QQ', bit_len >> 64, bit_len & M64)
    H = list(IV512)
    for i in range(0, len(padded), 128):
        W16 = struct.unpack('>16Q', padded[i:i+128])
        H, _, _, _ = compress_trace(H, W16, steps)
    return b''.join(struct.pack('>Q', h) for h in H)

def self_test():
    for msg in [b'', b'abc', b'abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmnhijklmnoijklmnopjklmnopqklmnopqrlmnopqrsmnopqrstnopqrstu', bytes(range(256))]:
        got = sha512_reduced(msg, 80).hex()
        exp = hashlib.sha512(msg).hexdigest()
        assert got == exp, f"Mismatch: {got} != {exp}"
    # Also test alternative-equation identities on a random block:
    import random
    random.seed(12345)
    W16 = [random.getrandbits(64) for _ in range(16)]
    H_out, A, E, W = compress_trace(IV512, W16, 36)
    for t in range(36):
        e_check = (A[t] + A[t-4] - S0(A[t-1]) - Maj(A[t-1], A[t-2], A[t-3])) & M64
        w_check = (E[t] - A[t-4] - E[t-4] - S1(E[t-1]) - Ch(E[t-1], E[t-2], E[t-3]) - K512[t]) & M64
        assert e_check == E[t], f"E[{t}] eq failed"
        assert w_check == W[t], f"W[{t}] eq failed"
    print("[OK] sha512_ref.py passed hashlib.sha512 (80-step) and alternative A/E/W identity checks!")

if __name__ == "__main__":
    self_test()
