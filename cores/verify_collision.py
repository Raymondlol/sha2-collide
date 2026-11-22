#!/usr/bin/env python3
"""Independent verifier for 2-block Standard-IV 36-step SHA-512 collisions & Step-18 samples."""
import json, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src", "ref"))
import sha512_ref as S

K = S.K512
IV = [
    0x6a09e667f3bcc908, 0xbb67ae8584caa73b, 0x3c6ef372fe94f82b, 0xa54ff53a5f1d36f1,
    0x510e527fade682d1, 0x9b05688c2b3e6c1f, 0x1f83d9abfb41bd6b, 0x5be0cd19137e2179,
]
MASK = (1 << 64) - 1
rot = lambda x, n: ((x >> n) | (x << (64 - n))) & MASK
S0  = lambda x: rot(x, 28) ^ rot(x, 34) ^ rot(x, 39)
S1  = lambda x: rot(x, 14) ^ rot(x, 18) ^ rot(x, 41)
s0  = lambda x: rot(x, 1) ^ rot(x, 8) ^ (x >> 7)
s1  = lambda x: rot(x, 19) ^ rot(x, 61) ^ (x >> 6)
Ch  = lambda x, y, z: (x & y) ^ ((~x & MASK) & z)
Maj = lambda x, y, z: (x & y) ^ (x & z) ^ (y & z)

def compress36(iv, M):
    W = list(M)
    for t in range(16, 36):
        W.append((s1(W[t-2]) + W[t-7] + s0(W[t-15]) + W[t-16]) & MASK)
    a, b, c, d, e, f, g, h = iv
    for t in range(36):
        t1 = (h + S1(e) + Ch(e, f, g) + K[t] + W[t]) & MASK
        t2 = (S0(a) + Maj(a, b, c)) & MASK
        h, g, f, e = g, f, e, (d + t1) & MASK
        d, c, b, a = c, b, a, (t1 + t2) & MASK
    return [(iv[i] + [a, b, c, d, e, f, g, h][i]) & MASK for i in range(8)]

def verify_record(seed, M1, M2, expected_cv=None):
    M0 = [seed, 0x5348413531323336, 0x434f4c4c49444521] + [0] * 13
    cv = compress36(IV, M0)
    if expected_cv is not None:
        assert cv == expected_cv, "CV mismatch vs Standard-IV F_36(IV, M0)!"
    H1 = compress36(cv, M1)
    H2 = compress36(cv, M2)
    assert M1 != M2, "Trivial M1 == M2!"
    if H1 == H2:
        print("SUCCESS: VERIFIED 2-BLOCK STANDARD-IV 36-STEP SHA-512 COLLISION!!")
        print("  M0[0] seed =", seed)
        print("  H =", " ".join(f"{x:016x}" for x in H1))
        return True
    print(f"Verified Standard-IV CV match (seed={seed}); not a full 36-step collision.")
    return False

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: verify_collision.py <log_file>")
        sys.exit(1)
    for line in open(sys.argv[1]):
        if line.startswith("[SAMPLE_STEP18_JSON]"):
            rec = json.loads(line.split(" ", 1)[1])
            verify_record(rec["seed"], [int(x,16) for x in rec["M1"]], [int(x,16) for x in rec["M2"]],
                          [int(x,16) for x in rec["cv"]])
