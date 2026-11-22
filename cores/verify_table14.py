#!/usr/bin/env python3
"""Standalone self-contained verifier for ePrint 2026/1120 Table 14 and our algebraic speedups."""
import json, math, os, random, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src", "ref"))
import sha512_ref as S

K = S.K512
MASK = (1 << 64) - 1
rot = lambda x, n: ((x >> n) | (x << (64 - n))) & MASK
S0  = lambda x: rot(x, 28) ^ rot(x, 34) ^ rot(x, 39)
S1  = lambda x: rot(x, 14) ^ rot(x, 18) ^ rot(x, 41)
s0  = lambda x: rot(x, 1) ^ rot(x, 8) ^ (x >> 7)
s1  = lambda x: rot(x, 19) ^ rot(x, 61) ^ (x >> 6)
Ch  = lambda x, y, z: (x & y) ^ ((~x & MASK) & z)
Maj = lambda x, y, z: (x & y) ^ (x & z) ^ (y & z)

def trace(cv, M, steps=36):
    W = list(M)
    for t in range(16, steps):
        W.append((s1(W[t-2]) + W[t-7] + s0(W[t-15]) + W[t-16]) & MASK)
    a, b, c, d, e, f, g, h = cv
    A = {-1: a, -2: b, -3: c, -4: d}
    E = {-1: e, -2: f, -3: g, -4: h}
    for t in range(steps):
        T1 = (E[t-4] + S1(E[t-1]) + Ch(E[t-1], E[t-2], E[t-3]) + K[t] + W[t]) & MASK
        T2 = (S0(A[t-1]) + Maj(A[t-1], A[t-2], A[t-3])) & MASK
        A[t] = (T1 + T2) & MASK
        E[t] = (A[t-4] + T1) & MASK
    H = [(cv[i] + [A[steps-1], A[steps-2], A[steps-3], A[steps-4],
                   E[steps-1], E[steps-2], E[steps-3], E[steps-4]][i]) & MASK for i in range(8)]
    return W, A, E, H

def main():
    spec = json.load(open(os.path.join(ROOT, "cores", "table14_36step_sha512.json")))
    CV = [int(x, 16) for x in spec["CV"]]
    M1 = [int(x, 16) for x in spec["M1"]]
    M2 = [int(x, 16) for x in spec["M2"]]
    EXP = [int(x, 16) for x in spec["expected_hash"]]

    W1, A1, E1, H1 = trace(CV, M1, 36)
    W2, A2, E2, H2 = trace(CV, M2, 36)
    assert M1 != M2 and H1 == H2 == EXP, "Table 14 36-step SFS collision mismatch!"
    print(f"[PASS 1/3] Table 14 36-step SHA-512 SFS verified: H1 == H2 == {H1[0]:016x}...")

    e15_fix_mask = (1 << 18) | (1 << 41) | (1 << 45)
    for b in range(64):
        c0 = ((Ch(0, E2[14], E2[13]) >> b) & 1) - ((Ch(0, E1[14], E1[13]) >> b) & 1)
        c1 = ((Ch(MASK, E2[14], E2[13]) >> b) & 1) - ((Ch(MASK, E1[14], E1[13]) >> b) & 1)
        if c0 != c1:
            e15_fix_mask |= (1 << b)
    tgt_p = (((A1[15]>>1)^(A1[15]>>6))&1, ((A1[15]>>0)^(A1[15]>>53))&1, ((A1[15]>>48)^(A1[15]>>54))&1)
    a14_bits = spec["verified_freedom"]["upper_W14_free_bits"]
    total_u = 0
    for m in range(1 << len(a14_bits)):
        a14 = A1[14]
        for i, b in enumerate(a14_bits):
            if (m >> i) & 1:
                a14 ^= (1 << b)
        t2_15 = (S0(a14) + Maj(a14, A1[13], A1[12])) & MASK
        dae15 = (A1[11] - t2_15) & MASK
        am = (1 << 59)
        for b in range(64):
            if b == 59:
                continue
            if (((Maj(0, a14, A2[13])>>b)&1) - ((Maj(0, a14, A1[13])>>b)&1)) != \
               (((Maj(MASK, a14, A2[13])>>b)&1) - ((Maj(MASK, a14, A1[13])>>b)&1)):
                am |= (1 << b)
        dp = {(0, 0, 0, 0): 1}
        for b in range(64):
            ndp = {}
            da_b = (dae15 >> b) & 1
            a_r = (A1[15] >> b) & 1 if ((am >> b) & 1) else -1
            e_r = (E1[15] >> b) & 1 if ((e15_fix_mask >> b) & 1) else -1
            for (c, p0, p1, p2), cnt in dp.items():
                for ab in (0, 1):
                    if a_r != -1 and ab != a_r:
                        continue
                    s = ab + da_b + c; eb = s & 1; nc = s >> 1
                    if e_r != -1 and eb != e_r:
                        continue
                    k = (nc, p0 ^ (ab if b in (1, 6) else 0),
                             p1 ^ (ab if b in (0, 53) else 0),
                             p2 ^ (ab if b in (48, 54) else 0))
                    ndp[k] = ndp.get(k, 0) + cnt
            dp = ndp
        total_u += sum(cnt for st, cnt in dp.items() if st[1:] == tgt_p)
    assert total_u == 33554432, f"Expected 2^25, got {total_u}"
    print(f"[PASS 2/3] Exact upper (W14, W15) freedom verified: N_U = {total_u:,} = 2^{math.log2(total_u):.3f}")

    rng = random.Random(2026)
    dW5, dW6, dW7, dW8 = [(W2[t]-W1[t])&MASK for t in (5,6,7,8)]
    t5, t6, t7, t8 = [(s0(W2[t])-s0(W1[t]))&MASK for t in (5,6,7,8)]
    orig_dCh6 = (Ch(E2[5], E2[4], E2[3]) - Ch(E1[5], E1[4], E1[3])) & MASK
    orig_dCh7 = (Ch(E2[6], E2[5], E2[4]) - Ch(E1[6], E1[5], E1[4])) & MASK
    T1_4 = (A1[4] - S0(A1[3]) - Maj(A1[3], A1[2], A1[1])) & MASK
    T1_5 = (A1[5] - S0(A1[4]) - Maj(A1[4], A1[3], A1[2])) & MASK
    T1_6 = (A1[6] - S0(A1[5]) - Maj(A1[5], A1[4], A1[3])) & MASK
    T1_7 = (A1[7] - S0(A1[6]) - Maj(A1[6], A1[5], A1[4])) & MASK
    T1_8 = (A1[8] - S0(A1[7]) - Maj(A1[7], A1[6], A1[5])) & MASK
    A0_base = (E1[4] - T1_4) & MASK

    while True:
        Am1, Am2, Am3, Am4 = [rng.getrandbits(64) for _ in range(4)]
        E3_n = (Am1 + A1[3] - S0(A1[2]) - Maj(A1[2], A1[1], A0_base)) & MASK
        if ((Ch(E2[5], E2[4], E3_n) - Ch(E1[5], E1[4], E3_n)) & MASK) != orig_dCh6: continue
        W7_n = (T1_7 - E3_n - S1(E1[6]) - Ch(E1[6], E1[5], E1[4]) - K[7]) & MASK
        if ((s0((W7_n + dW7) & MASK) - s0(W7_n)) & MASK) != t7: continue
        E2_n = (Am2 + A1[2] - S0(A1[1]) - Maj(A1[1], A0_base, Am1)) & MASK
        W6_n = (T1_6 - E2_n - S1(E1[5]) - Ch(E1[5], E1[4], E3_n) - K[6]) & MASK
        if ((s0((W6_n + dW6) & MASK) - s0(W6_n)) & MASK) != t6: continue
        cv_g = [Am1, Am2, Am3, Am4] + [rng.getrandbits(64) for _ in range(4)]
        break

    e4_bits = spec["verified_freedom"]["lower_E4_free_bits_29"]
    ok_all = 0; N = 2000
    for _ in range(N):
        e4 = E1[4]
        for b in e4_bits:
            if rng.getrandbits(1): e4 ^= (1 << b)
        w8 = (T1_8 - e4 - S1(E1[7]) - Ch(E1[7], E1[6], E1[5]) - K[8]) & MASK
        if ((s0((w8+dW8)&MASK) - s0(w8)) & MASK) != t8: continue
        if ((Ch(E2[6], E2[5], e4) - Ch(E1[6], E1[5], e4)) & MASK) != orig_dCh7: continue
        a0 = (e4 - T1_4) & MASK
        e3 = (cv_g[0] + A1[3] - S0(A1[2]) - Maj(A1[2], A1[1], a0)) & MASK
        e2 = (cv_g[1] + A1[2] - S0(A1[1]) - Maj(A1[1], a0, cv_g[0])) & MASK
        e1 = (cv_g[2] + A1[1] - S0(a0) - Maj(a0, cv_g[0], cv_g[1])) & MASK
        if ((Ch(E2[5], e4, e3) - Ch(E1[5], e4, e3)) & MASK) != orig_dCh6: continue
        w7 = (T1_7 - e3 - S1(E1[6]) - Ch(E1[6], E1[5], e4) - K[7]) & MASK
        if ((s0((w7+dW7)&MASK) - s0(w7)) & MASK) != t7: continue
        w6 = (T1_6 - e2 - S1(E1[5]) - Ch(E1[5], e4, e3) - K[6]) & MASK
        if ((s0((w6+dW6)&MASK) - s0(w6)) & MASK) != t6: continue
        w5 = (T1_5 - e1 - S1(e4) - Ch(e4, e3, e2) - K[5]) & MASK
        if ((s0((w5+dW5)&MASK) - s0(w5)) & MASK) != t5: continue
        ok_all += 1
        if ok_all == 1:
            Am1, Am2, Am3, Am4, Em1, Em2, Em3, Em4 = cv_g
            e0 = (Am4 + a0 - S0(Am1) - Maj(Am1, Am2, Am3)) & MASK
            w0 = (a0 - S0(Am1) - Maj(Am1, Am2, Am3) - Em4 - S1(Em1) - Ch(Em1, Em2, Em3) - K[0]) & MASK
            w1 = (A1[1] - S0(a0) - Maj(a0, Am1, Am2) - Em3 - S1(e0) - Ch(e0, Em1, Em2) - K[1]) & MASK
            w2 = (A1[2] - S0(A1[1]) - Maj(A1[1], a0, Am1) - Em2 - S1(e1) - Ch(e1, e0, Em1) - K[2]) & MASK
            w3 = (A1[3] - S0(A1[2]) - Maj(A1[2], A1[1], a0) - Em1 - S1(e2) - Ch(e2, e1, e0) - K[3]) & MASK
            w4 = (A1[4] - S0(A1[3]) - Maj(A1[3], A1[2], A1[1]) - e0 - S1(e3) - Ch(e3, e2, e1) - K[4]) & MASK
            Mt1 = [w0, w1, w2, w3, w4, w5, w6, w7, w8] + M1[9:]
            Mt2 = [w0, w1, w2, w3, w4, (w5+dW5)&MASK, (w6+dW6)&MASK, (w7+dW7)&MASK, (w8+dW8)&MASK] + M2[9:]
            Wt1, At1, Et1, _ = trace(cv_g, Mt1, 36)
            Wt2, At2, Et2, _ = trace(cv_g, Mt2, 36)
            assert At1[16] == At2[16] and ((Et2[16]-Et1[16])&MASK) == ((E2[16]-E1[16])&MASK)
            assert ((Wt2[21]-Wt1[21])&MASK) == ((W2[21]-W1[21])&MASK) and Wt1[20] == Wt2[20] and Wt1[22] == Wt2[22]
    print(f"[PASS 3/3] Lower golden-CV + 29-bit E4 verified: {ok_all}/{N} ({ok_all/N*100:.2f}%, 2^{29+math.log2(ok_all/N):.2f}/M0)")

if __name__ == "__main__":
    main()
