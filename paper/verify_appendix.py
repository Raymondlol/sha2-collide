import json, sys, hashlib, struct
M64 = (1 << 64) - 1
def primes(n):
    ps, x = [], 2
    while len(ps) < n:
        if all(x % p for p in ps if p * p <= x): ps.append(x)
        x += 1
    return ps
def iroot(n, k):
    lo, hi = 0, 1 << (n.bit_length() // k + 1)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        lo, hi = (mid, hi) if mid ** k <= n else (lo, mid - 1)
    return lo
P = primes(80)
K  = [iroot(p << 192, 3) & M64 for p in P]      # frac(cbrt p) * 2^64
IV = [iroot(p << 128, 2) & M64 for p in P[:8]]  # frac(sqrt p) * 2^64
rotr = lambda x, n: ((x >> n) | (x << (64 - n))) & M64
S0 = lambda x: rotr(x, 28) ^ rotr(x, 34) ^ rotr(x, 39)
S1 = lambda x: rotr(x, 14) ^ rotr(x, 18) ^ rotr(x, 41)
s0 = lambda x: rotr(x, 1) ^ rotr(x, 8) ^ (x >> 7)
s1 = lambda x: rotr(x, 19) ^ rotr(x, 61) ^ (x >> 6)
def compress(cv, M, steps):
    W = list(M) + [0] * (steps - 16)
    for i in range(16, steps):
        W[i] = (s1(W[i-2]) + W[i-7] + s0(W[i-15]) + W[i-16]) & M64
    a, b, c, d, e, f, g, h = cv
    for t in range(steps):
        t1 = (h + S1(e) + ((e & f) ^ (~e & g)) + K[t] + W[t]) & M64
        t2 = (S0(a) + ((a & b) ^ (a & c) ^ (b & c))) & M64
        h, g, f, e, d, c, b, a = g, f, e, (d + t1) & M64, c, b, a, (t1 + t2) & M64
    return [(x + y) & M64 for x, y in zip(cv, [a, b, c, d, e, f, g, h])]
def blocks(m):
    L = len(m) * 8
    m += b"\x80" + b"\x00" * ((112 - (len(m) + 1) % 128) % 128) + L.to_bytes(16, "big")
    return [list(struct.unpack(">16Q", m[i:i+128])) for i in range(0, len(m), 128)]
def sha512r(m, steps):
    cv = IV
    for B in blocks(m): cv = compress(cv, B, steps)
    return b"".join(x.to_bytes(8, "big") for x in cv)
for m in [b"", b"abc", bytes(range(256)) * 3]:
    assert sha512r(m, 80) == hashlib.sha512(m).digest()
J = json.load(open(sys.argv[1])); h = lambda L: [int(x, 16) for x in L]
M0, M1, M2 = h(J["M0"]), h(J["M1"]), h(J["M2"])
cv = compress(IV, M0, 36)
assert M1 != M2 and compress(cv, M1, 36) == compress(cv, M2, 36)
a = b"".join(x.to_bytes(8, "big") for x in M0 + M1)
b = b"".join(x.to_bytes(8, "big") for x in M0 + M2)
assert sha512r(a, 36) == sha512r(b, 36) and hashlib.sha512(a).digest() != hashlib.sha512(b).digest()
print("36-step SHA-512 collision verified:", sha512r(a, 36).hex())
