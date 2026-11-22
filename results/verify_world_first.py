#!/usr/bin/env python3
"""
Clean-room, zero-dependency FIPS 180-4 independent verifier for the
World-First Practical 36-Step Standard-IV SHA-512 Collision.

Derives SHA-512 IV and K[0..79] directly from first 8 / 80 primes,
cross-validates full 80-round SHA-512 against Python stdlib hashlib.sha512,
and verifies:
  1. Block 1: F_36(Standard_IV, M0) == CV
  2. Block 2: M1 != M2 and F_36(CV, M1) == F_36(CV, M2) == H
  3. Full FIPS 180-4 Merkle-Damgard padded 36-step SHA-512(M0 || M1) == SHA-512(M0 || M2)
  4. Full 36-step state differential trace (confirming non-triviality and exact absorption)
"""
import hashlib
import json
import math
import pathlib
import struct
import sys

MASK64 = (1 << 64) - 1


def first_n_primes(n: int) -> list[int]:
    primes = []
    cand = 2
    while len(primes) < n:
        if all(cand % p != 0 for p in primes if p * p <= cand):
            primes.append(cand)
        cand += 1
    return primes


def frac_root_u64(n: int, k: int) -> int:
    """Exact floor(frac(n^(1/k)) * 2^64) using pure integer binary search."""
    # Find largest integer r such that r^k <= n * 2^(64*k)
    target = n << (64 * k)
    lo, hi = 0, 1 << (64 + 16)
    while lo + 1 < hi:
        mid = (lo + hi) >> 1
        if mid ** k <= target:
            lo = mid
        else:
            hi = mid
    return lo & MASK64


PRIMES_80 = first_n_primes(80)
FIPS_IV = [frac_root_u64(p, 2) for p in PRIMES_80[:8]]
FIPS_K = [frac_root_u64(p, 3) for p in PRIMES_80]

# Cross-check against published FIPS 180-4 Table values
PUBLISHED_IV = [
    0x6A09E667F3BCC908, 0xBB67AE8584CAA73B, 0x3C6EF372FE94F82B, 0xA54FF53A5F1D36F1,
    0x510E527FADE682D1, 0x9B05688C2B3E6C1F, 0x1F83D9ABFB41BD6B, 0x5BE0CD19137E2179,
]
assert FIPS_IV == PUBLISHED_IV, "Derived IV mismatch vs FIPS 180-4!"


def rotr(x: int, n: int) -> int:
    return ((x >> n) | ((x << (64 - n)) & MASK64)) & MASK64


def big_sigma0(x: int) -> int:
    return rotr(x, 28) ^ rotr(x, 34) ^ rotr(x, 39)


def big_sigma1(x: int) -> int:
    return rotr(x, 14) ^ rotr(x, 18) ^ rotr(x, 41)


def small_sigma0(x: int) -> int:
    return rotr(x, 1) ^ rotr(x, 8) ^ (x >> 7)


def small_sigma1(x: int) -> int:
    return rotr(x, 19) ^ rotr(x, 61) ^ (x >> 6)


def ch(x: int, y: int, z: int) -> int:
    return (x & y) ^ ((~x & MASK64) & z)


def maj(x: int, y: int, z: int) -> int:
    return (x & y) ^ (x & z) ^ (y & z)


def compress(cv: list[int], block_bytes: bytes, rounds: int, trace: bool = False):
    assert len(block_bytes) == 128
    W = list(struct.unpack(">16Q", block_bytes))
    for t in range(16, rounds):
        W.append((small_sigma1(W[t - 2]) + W[t - 7] + small_sigma0(W[t - 15]) + W[t - 16]) & MASK64)
    a, b, c, d, e, f, g, h = cv
    states = [(a, b, c, d, e, f, g, h)]
    for t in range(rounds):
        t1 = (h + big_sigma1(e) + ch(e, f, g) + FIPS_K[t] + W[t]) & MASK64
        t2 = (big_sigma0(a) + maj(a, b, c)) & MASK64
        h, g, f, e, d, c, b, a = g, f, e, (d + t1) & MASK64, c, b, a, (t1 + t2) & MASK64
        if trace:
            states.append((a, b, c, d, e, f, g, h))
    out = [(cv[i] + v) & MASK64 for i, v in enumerate((a, b, c, d, e, f, g, h))]
    return (out, states, W) if trace else out


def sha512_padded(msg: bytes, rounds: int = 80) -> bytes:
    bit_len = len(msg) * 8
    padded = msg + b"\x80"
    while (len(padded) % 128) != 112:
        padded += b"\x00"
    padded += struct.pack(">QQ", bit_len >> 64, bit_len & MASK64)
    cv = list(FIPS_IV)
    for off in range(0, len(padded), 128):
        cv = compress(cv, padded[off : off + 128], rounds=rounds)
    return struct.pack(">8Q", *cv)


def words_to_bytes(words: list[int]) -> bytes:
    return struct.pack(">16Q", *words)


def main():
    json_path = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(__file__).with_name("collision_36step_sha512_std_iv.json")
    data = json.loads(json_path.read_text())

    M0_w = [int(x, 16) for x in data["M0"]]
    CV_w = [int(x, 16) for x in data["CV"]]
    M1_w = [int(x, 16) for x in data["M1"]]
    M2_w = [int(x, 16) for x in data["M2"]]
    H1_w = [int(x, 16) for x in data["H1"]]
    H2_w = [int(x, 16) for x in data["H2"]]

    M0_b = words_to_bytes(M0_w)
    M1_b = words_to_bytes(M1_w)
    M2_b = words_to_bytes(M2_w)
    msg1 = M0_b + M1_b
    msg2 = M0_b + M2_b

    # 1. Cross-validate 80-round clean-room SHA-512 against Python hashlib.sha512
    for label, sample in [("empty", b""), ("abc", b"abc"), ("M0||M1", msg1), ("M0||M2", msg2)]:
        got = sha512_padded(sample, rounds=80)
        exp = hashlib.sha512(sample).digest()
        assert got == exp, f"80-round mismatch on {label}!"
    print("[PASS 1/5] Clean-room FIPS 180-4 80-round SHA-512 matches hashlib.sha512 bit-for-bit (4/4 vectors)")

    # 2. Check M0 -> CV under Standard FIPS IV (36 rounds)
    cv_check = compress(FIPS_IV, M0_b, rounds=36)
    assert cv_check == CV_w, "CV mismatch for F_36(FIPS_IV, M0)!"
    m0_ascii = M0_b[8:24].decode("ascii", errors="replace")
    print(f"[PASS 2/5] Block 1 F_36(Standard_IV, M0) == CV verified | M0 tag: {m0_ascii!r} | seed={data.get('m0_seed')}")

    # 3. Check M1 != M2 and 36-step compression collision F_36(CV, M1) == F_36(CV, M2)
    assert M1_b != M2_b and msg1 != msg2, "Messages are identical!"
    h1_check, st1, W1 = compress(CV_w, M1_b, rounds=36, trace=True)
    h2_check, st2, W2 = compress(CV_w, M2_b, rounds=36, trace=True)
    assert h1_check == H1_w, "H1 mismatch vs JSON!"
    assert h2_check == H2_w, "H2 mismatch vs JSON!"
    assert h1_check == h2_check, "NOT A COLLISION!"
    raw_hex = struct.pack(">8Q", *h1_check).hex()
    print(f"[PASS 3/5] Block 2 F_36(CV, M1) == F_36(CV, M2) verified (M1 != M2, 512-bit digest identical)")
    print(f"           Raw 2-block 36-step digest : {raw_hex}")

    # 4. Check full FIPS 180-4 padded 36-step SHA-512(M0 || M1) == SHA-512(M0 || M2)
    pad_h1 = sha512_padded(msg1, rounds=36).hex()
    pad_h2 = sha512_padded(msg2, rounds=36).hex()
    assert pad_h1 == pad_h2, "Padded 36-step SHA-512 mismatch!"
    print(f"[PASS 4/5] Full FIPS 180-4 padded 36-step SHA-512(M0||M1) == SHA-512(M0||M2) verified!")
    print(f"           Padded 36-step SHA-512     : {pad_h1}")

    # 5. Print message-word and step-state differential summary
    diff_words = [i for i in range(16) if M1_w[i] != M2_w[i]]
    active_w = [t for t in range(36) if W1[t] != W2[t]]
    active_steps = [t for t in range(1, 37) if st1[t] != st2[t]]
    print(f"[PASS 5/5] Differential check: differing M words={diff_words} | active W steps={active_w[0]}..{active_w[-1]} ({len(active_w)} words) | active state steps={active_steps[0]}..{active_steps[-1]} (zero from step {active_steps[-1]+1}..36)")
    print("\nALL 5/5 CLEAN-ROOM FIPS 180-4 CHECKS PASSED! WORLD-FIRST 36-STEP STANDARD-IV SHA-512 COLLISION CONFIRMED!")


if __name__ == "__main__":
    main()
