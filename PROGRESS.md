# Progress Log — `sha2-collide`

All runs on `g916` (AMD Ryzen 9 9950X3D 16C/32T, 61 GiB RAM, 1× NVIDIA RTX 5090 32 GiB GDDR7 / 96 MB L2, CUDA 13.0).

---

## Phase 1 — Initial Audits & RTX 5090 Pipeline Benchmarks (2026-11-21 to 2026-11-22 AM)
*Archived in full in `archive/AUDIT_SUMMARY_2026-11-22.md`.*

- Audited SHA-512 36-step (`i=1`) and 38-step (`i=3`) differential structures (`k=48`, `k=63`).
- Built and benchmarked 12-step core SAT solvers from scratch (Signed-CSA `V2` + `R4` + 4-round-collapsed `C7b`, shrinking CNF `422k -> 190k` clauses; 16T clause-sharing CaDiCaL + Kissat 4.0.4). Spent `32.7 CPU-hr` (~500M Kissat conflicts) and mapped the exact scaling cliff (`i=1` 8-step `7.43s`, 9-step `27.5s–53.7s`), confirming finding a 12-step core from scratch via blind CDCL hits the same wall as Goldzweig (`ePrint 2026/2913`, 52.3k vCPU-hr on 32-bit).
- Built & bit-for-bit verified RTX 5090 GPU building blocks (`src/gpu/`):
  - `bench_phase1_gpu.cu`: real Standard-IV `F_36(IV, M0)` + 256 `E0[4]` cands -> **`1.85 G M0/s`**, **`473.9 G (M0, E4)/s`**.
  - `sha512_36_gpu.cu`: 23-bit split-array L2 table (`64 MB` in `96 MB` L2) -> **`126.09 G pi/s`**, **`27.98 G s17==0 / s`**, **`265k` full Step-18/s**.

## Phase 2 — Literature Core Discovery & Verification (`ePrint 2026/1120` Table 14) (2026-11-22 Noon)

- Inspected `ePrint 2026/1080` (CRYPTO 2026: practical 35-step SHA-512, theoretical `2^{136}` 36-step SHA-512) and `ePrint 2026/1120` (practical 36-step SHA-256 `2^{57}`, theoretical `2^{66.2}` 36-step SHA-512).
- Verified `ePrint 2026/1120` **Table 14 (p. 27)** 36-step SHA-512 SFS collision bit-for-bit against `sha512_ref.py`:
  - `H1 == H2 == 0xf121c6371d1f7dae...` (`PASS`).
  - Already aligned at `i=1` (`t=5..16`, `k=59`), `dE15=0`, `dA15=-2^{59}`, `dA16=0`, `wt(dE16)=3`, `wt(dE17)=3`, `dE19=-2^{59}`.
  - Real SHA-512 modular tail (`N=500,000`):
    - Step 17: `15.0b` (`12` direct bits + `2` XOR pairs on `E16`, plus `A16[59] == A14[59] == 0`)
    - Step 18: `16.0b`
    - Steps 19..35: `7.00 + 6.84 + 1.00 + 1.00 + 11.51 + 3.83 = 31.18b`
    - Total raw tail = **`62.18b`** (vs strict `66.0b` in paper).
- Also verified `ePrint 2026/1120` Table 20 (37-step SHA-512 SFS, `H1 == H2 == 0xe28c6359c9f80b41...`), but its `2^{87}` tail keeps it theoretical.

## Phase 3 — Step A: Rigidity Proof & Pure Algebraic Breakthroughs (2026-11-22 Afternoon)

1. **SAT Rigidity Check (`src/sat/expand_t14_upper.cpp`, `expand_t14_v2.cpp`, `check_t14_freedom.cpp`)**:
   - Proved with CaDiCaL that fixing `t=5..8` concrete + Table 14's `signDiffs` locks `W9..W13` uniquely (`Enum stopped at k=1, UNSAT in 0.000s`, `|xE13| <= 33 UNSAT in 0.000s`) because `E13`'s 20-bit borrow chain (`0x6fffff -> 0x700001`) couples with `dE10..dE12`.
2. **Algebraic Breakthrough #1 — Exact Upper Freedom `N_U = 2^{25.000}`**:
   - Exact bit-by-bit carry DP (`cores/verify_table14.py`) showed `W15` has **`2,097,152 = 2^{21.000}`** valid values per `W14` (5 `A15/S0` conditions redundant with `E15[41..60]=0`) and `W14` has **`16 = 2^4`** valid values (`[28, 40, 62, 63]`) -> **`N_U = 33,554,432 = 2^{25.000}`** (`16×` larger than `2026/1120`'s reported `2^{21}`).
   - Verified `1,600 / 1,600` random samples bit-for-bit against full 36-step SHA-512.
3. **Algebraic Breakthrough #2 — Lower `E4` Low-Bit Independence (`12,000×` Step-2 Speedup)**:
   - Discovered that once a golden `M0` satisfies `E3 + W7 + W6` (`2^{-14}`), varying `E4` across its 29 free bits (`[0..36, 63]`) preserves `E3 = 100.00%`, `W7 = 100.00%`, `W8+Ch7 = 91.0%`, `W6 = 80.4%`, needing only `W5` (`6.19%`).
   - Total pass rate = **`4.62%` (`2^{-4.44}`)** -> **`2^{24.56}` (~24.7M) valid lower solutions per golden `M0`** with uniform 64-bit `Lo = s0(W1) + W0`.
4. **Algebraic Breakthrough #3 — `k = 19` Sorted-Bucket `16 × 64` Shared-Memory Tiling**:
   - Computed exact guard-bit coverage vs speed tradeoff (`g = 0..10`):
     - `g = 6, k = 19` (`524,288` buckets, `16` lower entries/bucket = `128 B` = 1 L2 cache line, `5.93%` full Step-17 hit rate): `~10.25h` @ `100 G/s` (`~6–9h` with `16×64` shared-memory tiling).
     - `g = 7, k = 20` (`8` lower entries/bucket, `7.91%` full Step-17 hit rate): `~9.77h` @ `100 G/s` (`T1 = 2ms`).

## Phase 4 — Next: Step B (RTX 5090 End-to-End Solver)
- [ ] Implement `src/gpu/sha512_36_t14_gpu.cu` with sentinel self-test (injecting Table 14 to verify end-to-end collision capture).
- [ ] Run 30s queued `gpu-run` verification & throughput benchmark on RTX 5090.
- [ ] Launch full checkpointed `gpu-run` collision search for world-first practical 36-step Standard-IV SHA-512 collision.
