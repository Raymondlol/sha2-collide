# Progress Log — `sha2-collide`

> **Note (2026-11-26).** Yuan and Yu published an independent 36-step SHA-512 collision on 2026-11-18 (https://github.com/gnyst001/sha512-36-collision; ePrint 2026/3005), four days before this one. This log is kept as written at the time, apart from this note and the renamed verifier.

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

## Phase 4 — Step B: Real Silicon Verification & v2 Optimization on RTX 5090 (2026-11-22)

- [x] Implemented `src/gpu/sha512_36_t14_gpu.cu` (compact active-bucket L2 table + `<=64-U` load-balanced chunks + real Standard-IV `M0` + 29-bit `E4` + full 36-step verifier).
- [x] **v1 Benchmark (`c4ee760`, 30.00s via `gpu-run`)**:
  - Sentinel self-test: **PASS 100%** (Table 14 injected -> exact 36-step collision captured `H1 == H2 == 0xf121c6371d1f7dae...`).
  - Live Standard-IV sample (`seed=953`) independently bit-for-bit verified against `sha512_ref.py` through Round 18 (`dA16=dA17=dA18=0`, `dE18=0`).
  - Speed: `5.78 ms/GM0` (`2.43 ms` Phase 1 + `3.33 ms` Phase 2), `228.08 G pairs/s` wall (`395.5 G pairs/s` pure Phase 2), `31.33 k s18/s` -> **`21.57 hours` ETA**.
  - **Key diagnosis from live silicon**: measured `P(s18 | s17_full) = 2^{-18.72}` revealed **3 forced Step-18 direct bits on `E16`** (`E16[43]=1, E16[51]=0, E16[57]=0`, 100% necessary to cancel `S1(2^{57})` and avoid `Ch18` carry in Step 18!).
- [x] **v2 Optimization (`aa31b0e`, 20.00s via `gpu-run`)**:
  - Added `[43:1, 51:0, 57:0]` into `E16_DIRECT_MASK_15` + switched to optimal **`NEW15_k19_g3`** key (`61,365` active buckets, `30.0 MB` hot L2 in 96 MB L2, `99.0%` capture) + exact 4-bit `W5` identity in Phase 1 + 16-way `low4` `E4` unroll.
  - Sentinel self-test: **PASS 100%**.
  - **Real Silicon Measured Speed**:
    - Timing: **`4.81 ms/GM0`** (`2.06 ms` Phase 1 [42.8%] + `2.74 ms` Phase 2 [57.0%], `GM0` = `0.02%`).
    - Throughput: **`273.58 G pairs/s` wall-clock** (**`480.0 G pairs/s` pure Phase 2**), **`66.04 k full Step-18/s`** (**`2.11×` speedup over v1**, `P(s18 | pre) = 2^{-15.72}` exact).
    - Live tail progression: `s18=1,320,815 -> s19=18,637 (2^{-6.15}) -> s20=85 (2^{-7.78}) -> s21=37 -> s22=21`.
    - **Measured Expected Time to World-First 36-Step Standard-IV SHA-512 Collision: `10.23 HOURS` (`0.43 days`, median 50% hit = `7.09 hours`)!!**

### Exact 4B-Sample Tail Calibration (`60ddcec`, 2.7s on 32T CPU)
- Row 23 (`W21`): `11.5114 ± 0.0012b` (`1,370,221 / 4e9`); Row 25 (`W23`): `3.8301 ± 0.0001b` (exact `-log2(9/128)`).
- Joint `Row 25 | Row 23` with real constrained `W16/W8/W7`: `3.8243 ± 0.0046b` (zero negative correlation).
- Unconstrained Steps 19–22 (`487.7M` real transitions): `16.1900 ± 0.0179b` (`17.7%` `E17` carries -> `+0.35b` vs fixed-`E17` estimate).
- **Exact total tail after Step 18: `31.5315 ± 0.0179` bits** (`95% CI: [31.4957, 31.5673]`) -> RTX 5090 v2 (`66.04 k s18/s`): **`13.06 h`** expected (`9.05 h` median).

- [ ] Launch full checkpointed `gpu-run` collision search (`./src/gpu/sha512_36_t14_gpu --search --seed 0`).

## Phase 5 — Production Search & WORLD-FIRST COLLISION ACHIEVED! (2026-11-22)

- [x] Deployed `v2` (`ab7c13a`) across 3 Blackwell GPUs with disjoint seed spaces:
  - `server-B` GPU 1 (RTX PRO 6000 Blackwell Server Edition, `--seed 0`): `11,220.4s` (`3.12h`), `2,170,000` GM0s, `686,216,835` s18 (`61.16 k/s`), `9,115` s22, `W21=1`, `W23=0`
  - `server-B` GPU 3 (RTX PRO 6000 Blackwell Server Edition, `--seed 1000000000000`): `11,215.4s` (`3.12h`), `2,170,000` GM0s, `686,137,818` s18 (`61.18 k/s`), `9,131` s22, `W21=3`, `W23=0`
  - `g916` GPU 0 (RTX 5090, `--seed 2000000000000`): **`10,113.01s` (`2.81h`), `2,102,288` GM0s, `665,107,101` s18 (`65.77 k/s`), `8,845` s22, `W21=5`, `W23=1`, `COL=1`!!**
- [x] **Combined Effort at Discovery**: **`9.04 GPU-hours`** (`~3.1 h` wall-clock), `2,037,461,754` Step-18 passes (`65.6%` of mean `3.105B`, right at median `2.15B`), `27,091` Step-22 passes -> `9` `W21` hits (`2^{-11.55}`, matching 4B-sample calibration `2^{-11.5114}` to `0.04` bits!) -> `1` `W23` hit -> **1 REAL COLLISION**.
- [x] **Independent Clean-Room Verification (`results/verify_collision.py`)**: **5/5 PASS**:
  1. Prime-derived FIPS 180-4 80-round SHA-512 matches `hashlib.sha512` bit-for-bit on 4 vectors (`""`, `"abc"`, `M0||M1`, `M0||M2`).
  2. `F_36(Standard_IV, M0) == CV` (`m0_seed = 2068887288786`, `M0` tag `"SHA51236COLLIDE!"`).
  3. `M1 != M2` (7 differing words `[5, 6, 7, 8, 9, 13, 14]`), `F_36(CV, M1) == F_36(CV, M2) == 7f987d4f13fdfebe...`.
  4. Full FIPS 180-4 padded 256-byte message collision: `SHA512_36(M0 || M1) == SHA512_36(M0 || M2) == 9eabd799baf71c41...`.
  5. State differential active steps `6..23`, zero identically from step `24..36`.
- [x] Permanently archived in `results/collision_36step_sha512_std_iv.json` & `results/verify_collision.py`; all 3 GPUs cleanly stopped (`0%` idle).
