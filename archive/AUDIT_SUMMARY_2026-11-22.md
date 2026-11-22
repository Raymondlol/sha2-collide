# SHA-512 (64-bit) Standard-IV Collision Audit & Overnight Campaign Summary (2026-11-22)

Directory: `~/Workspace/train-domains/sha512/`
Hardware: AMD Ryzen 9 9950X3D (16C/32T, 61 GiB RAM) + 1× NVIDIA GeForce RTX 5090 (32 GiB GDDR7, 96 MB L2, `sm_120`). All GPU runs queued strictly via `~/.local/bin/gpu-run`.
Solvers installed in `tools/`: CaDiCaL 3.1.0, Kissat 4.0.4, CryptoMiniSat 5.11.21.

---

## 1. Direction 2 (Completed): Non-4-Word 38-Step SHA-512 `i=3` Family Audit
- Script: `audit9_sha512_i3_38step.py` (7-word `i=3`: `{W7..W11, W15, W16 -> W23, W25}`, `dW26..W37 = 0`, 1-bit `dW16 = 2^k`).
- **Measured expansion + state-tail estimate (Steps 19..27, best `k=4`)**:
  - `d15 = 0xf000000000000008` (`wt=2`), `d23 = 0xffd0000000000004` (`wt=3`), `p25 = 8.14b`, `p27 = 2.81b` (`exp = 10.95b`), estimated Steps 19..27 tail `~53.9` bits (before Step 18 / non-free `W16` overhead).
  - Compared to 4-word `j=4` 38-step (`~72.8b` tail): `i=3` shaves `~18` bits on the tail, **at the cost of a 12-step core with non-free `dW16 != 0`**, putting both 38-step SHA-512 families beyond single-machine core+tail reach right now.

---

## 2. Direction 1: 36-Step Standard-IV SHA-512 (`i=1` & 4-word `j=4`)

### 2.1 New Ultra-Sparse `i=1` Parameters (`{W5..W9, W13, W14 -> W21, W23}`, `dW24..W35 = 0`)
- **`k=48` ultra-sparse**: `d14 = +2^48`, `d5 = d7 = d23 = -2^48`, `d13 = -d6 = -2^47 - 2^40` (`wt=2`), `a = +2^46 + 2^33 + 2^32` (`wt=3`), `d9` (`wt=3`), `d8` (`NAF wt=6`), `|V8| = 2^{47.0}`, `|V9| = 2^{52.0}`, `|V13| = 2^{58.0}`, `|V14| = 2^{60.0}`, tail after Step 18 `~33.6` bits.
- **`k=63` zero-carry MSB**: `d14 = d5 = d7 = d23 = 2^63` (`+2^63 ≡ -2^63 mod 2^64`, **zero carry at MSB**, `dW5 + dW14 == 0` automatic!), `s0(2^63) = {55, 56, 62}` (`wt=2` subset `{55, 62}`), `a = {47, 49, 61}` (3 non-adjacent bits → all 8 sign patterns free), `s1(2^63) = {2, 44, 57}`, `s0(s1(2^63))` only 8 raw bits (`wt=5–6`), `|V5| = |V7| = |V14| = 2^{61.0}`, `|V9| = 2^{53.2}`, `|V8| = 2^{45.0}`, exact MC tail after Step 18 = **`34.77` bits**.
  - Discovered critical Step-18 liveness condition on `k=63` bit 49 (shared by `dA13` and `dE16`): `E0[15][49]` must equal `E1[16][49]` (if `dA13[49]=+1`) or `E0[16][49]` (if `dA13[49]=-1`); encoded directly as 4 clauses in `gen_sha512_C7b_cnf.cpp`.

### 2.2 SAT Formulation & Solver Engineering (`gen_sha512_cnf_v3.cpp`, `gen_sha512_k63_cnf.cpp`, `gen_sha512_C7b_cnf.cpp`, `gen_sha512_j4_cnf.cpp`, `par_cadical_share.cpp`)
- **Strict empirical hierarchy** (measured head-to-head on 64-bit SHA-512 & 32-bit Omer SHA-256):
  - **`V2` (Concrete `Eq1_0` + Balanced Signed-CSA `dEq1`/`dEq2` + exact `CL_XOR_D`/`CL_MAJ_D`/`CL_CH_D`) + `--sat`** >> `V3` (compact 2-exec + diff clauses) >> plain `2ex`, `CryptoMiniSat 5` native XOR Gauss-E, `LLW-IPASIR-GF2`, and diff-first `MNS11` (which hit `1,008/1,008` Stage-2 UNSAT even on 8-step `i1-t9` due to concrete bit/linear over-constraining `~1,350 > 1,024` DOF).
- **Exact Algebraic Reductions (`R4` + 4-Round Collapse `C7b`)**:
  - `R4` eliminates `Eq1_0[6, 7, 15, 16]` losslessly (verified on `k63_t8_Ksat2.out`: exact 64-step DP counts `2^{42.00}` valid upper pairs `(A15, E15)` satisfying `Eq1_0[15], dEq1[16], dEq2[16]`).
  - **`C7b`** additionally collapses `t=6`, `t=7` (for `k=63`, `E0[6][63]=1` makes `dCh7[63]` cancel `dW7=2^63` mod `2^64`, shrinking `dE[7]` to **9 exact bits** and `dA[7]` to **23 exact bits** with zero signed-CSA!), `dEq1[15]`, and `dEq1[16]`.
  - Result: full 12-step core CNF shrank **2.22×** from `422,388` clauses (`k63_R4`) to **`190,485` clauses** (`k63_C7b`, `nv=30,068`, Kissat active vars `28.3k → 13.0k`), raising initial conflict rate **2.3×** (~`10k → 23.5k` conflicts/s initial, `10.3–13.8k/s` sustained over 40 min). Verified 100% exact against `sha512_ref.py` on `200/200` random concrete states and on `k63_C7b_t11` (`1.94s` SAT).
- **16-Thread Native Clause-Sharing CaDiCaL 3.1.0 (`par_cadical_share.cpp`)**:
  - Uses `Learner` (`len <= 4`), `FixedAssignmentListener` (units), `ExternalPropagator` (primary vars `1..4000/6500`), and `Terminator`. Calibrated 14.5% faster than 16-thread no-share on 30s `i1_t8` (`29.26s` vs `33.50s`).

### 2.3 Measured Core-Search Scaling Wall (All SAT Results Independently Verified vs `sha512_ref.py`)
| Instance | Width / Active Core | Best Solver & Measured Time | Notes |
|---|---|---|---|
| 32-bit Omer SHA-256 `t=8..16` | 9 steps | **2.49s** (`Kissat --sat`), 7.17s (CMS5), 8.81s (CaDiCaL) | Baseline calibration |
| 32-bit Omer SHA-256 `t=7..16` | 10 steps | **34.16s** (`Kissat --sat`) | 1 round shy of Omer `t=6..16` (52,300 vCPU-hr) |
| 64-bit `j=4` 4-word `t=11..16` / `t=10..16` | 6 / 7 steps | **0.54s** / **10.73s** (`Kissat --sat`), `t=9..16` >20s | Dense explosion at `t=11` (`wA11=32`) → 1–2 rounds harder than `i=1` |
| 64-bit `i=1 k=48` `t=9..16` (`s17=0`) | 8 steps | **7.43s** (`CaDiCaL --sat`), 8.91s (`Kissat --sat`) | Verified 100% PASS |
| 64-bit `i=1 k=48` `t=8..16` (`s17=0`) | 9 steps | **27.54s** (`CaDiCaL --sat`) | Verified 100% PASS |
| 64-bit `i=1 k=63` `t=8..16` (`s17=0`) | 9 steps | **53.73s** (`Kissat --sat`), 58.71s (`CaDiCaL --sat`) | Verified 100% PASS; freeing masks (`k63v2_t8`) hurt (>60s) |
| 64-bit `i=1 k=63 C7b` `t=11..16` (+collapsed `t=6,7`) | 8+4 collapsed | **1.94s** (`Kissat --sat`) | Exact match at `t=5,6,7,11..16` + `s17=0` |
| **Full 12-step 36-step core (`i1_R4`, `k63_R4`, `k63_C7b`, `k48_C7b`)** | **12 steps (`t=5..16` + `s17=0`)** | **NO SAT in 32.7 CPU-hours** | **Wave 1 (14.0 CPU-hr, 24T)**: 157.5M Kissat conflicts + 44.7k shared CaDiCaL clauses. **Wave 2 (18.7 CPU-hr, 28T, `C7b`)**: 341.9M Kissat conflicts + 10.8k shared CaDiCaL clauses. |

**Conclusion on Core Bottleneck**: Adding the final bridge connecting sparse `t=6,7` to `t=8..16` (`t=7` coupling through `S0/S1/Maj/Ch` + `Eq1_0[8..14]`) hits the exact same exponential CDCL cliff observed by Goldzweig on 32-bit 38-step SHA-256 (`52,300` vCPU-hours on Arm) and explains why CRYPTO 2026 (`ePrint 2026/1080`) kept 36-step SHA-512 `i=1` theoretical.

---

## 3. Verified RTX 5090 Stage-3 GPU Pipeline (Bit-for-Bit Verified vs `sha512_ref.py`)
Once **any** verified 12-step `i=1` core exists, Stage 3 is **100% practical on a single consumer RTX 5090** using a **Reversed 23-bit Guarded Split-Array L2-Cache Table**:
- Both `k=63` and `k=48` have all 12 fixed `E16` bits of `s17==0` inside 4 non-adjacent bit clusters of total width **20 bits** + **3 linear `S1(E16)` parity bits** = **23-bit key** (`64 MB` hot `tbl_c16[1<<23]` fits 100% inside the RTX 5090's **96 MB L2 cache**; cold `tbl_w17` read only on `s17==0`).
- **Phase 1 (`bench_phase1_gpu.cu`, real Standard-IV `F_36(IV, M0)` + 256 `E0[4]` lower cands, `gpu-run`)**:
  - **`1.85 G` real `M0/s`**, **`473.9 G` `(M0, E0[4])` pairs/s**, **`57.85 M` full `V5+V6+V7` valid lower solutions/s**, **`226k` 8-bit-guarded L2 entries/s** → fills a `2^{20}` L2 table in **`4.64 seconds`**. Emitted `SAMPLE_LOWER` verified bit-for-bit via `sha512_ref.py compress_trace(IV -> M0 -> H1 -> W0..W7 -> A0..A3, E0..E7)`.
- **Phase 2 (`sha512_36_gpu.cu`, 23-bit split-array L2 table streaming Gray-code `π=(A15, E15)` over `2^{37}` valid upper pairs, `gpu-run`)**:
  - **`126.09 G π/s`** (`7.9 ps/π`), **`27.98 G` full `s17==0` passes/s** (`hit/π = 22.2%`, a `4,000×` speedup over raw `2^{-14.15}`), **`265k` full Step-18 passes/s** (`P(s18|s17) = 2^{-16.69}`). Emitted `SAMPLE_HIT` and `SELFTEST H1` verified bit-for-bit vs `sha512_ref.py`.
  - **Projected single-RTX-5090 wall time to 1 expected 36-step Standard-IV SHA-512 collision** (given 1 verified core):
    - **`k=48` (`2^{33.6}` tail)**: **`13.6 hours`** (`~6.8 hours` adding 1-bit `A16[k]==A14[k]` filter)
    - **`k=63` (`2^{34.77}` tail)**: **`30.7 hours`** (`~15.3 hours` with `A16[k]` filter)
    - Compare Goldzweig `2026/2913` on 32-bit 38-step SHA-256: `440 GPU-hours` expected on enterprise RTX PRO 6000 (`5.44 GB` DRAM table).

---

## 4. Previous Proven Structural Findings (2026-11-22 Early Audit)
1. **Why 1-bit `dW15` fails 4-word MNS13 (`{W7,W8,W10,W15 -> W23,W24}`) on 64-bit SHA-512**: `rotr(8) ^ (>>7)` cancels interior carry chains; min weight 26, min span 57 → `L8 >= 52` (`104 > 64` bits on `W8` → 0 solutions).
2. **Width-8 exact solver for `s0_512(z) - z == target (mod 2^64)`** (`audit1_sha512_exp.py`): 256 iterations (~100 ns, 0 bytes RAM; avg 1.026 exact preimages).
3. **1024-state exact DP for `Delta s0_512(x; d) == tau` over all `2^64` inputs** (`audit6_sha512_exact.c`, < 1 ms/query).
4. **Tail budgets (Standard-IV SHA-512; current practical record = 35 steps)**:
   - 36-step `i=1` (`k=48`/`k=63`): `33.6` / `34.77` bits after Step 18 → **13.6–30.7 hr on 1× RTX 5090** (core is the sole remaining bottleneck).
   - 36-step 4-word `j=4`: `~37.8` bits tail, but core is 1–2 rounds harder for SAT than `i=1`.
   - 37-step 4-word `j=4`: `~53–56.8` bits tail.
   - 38-step `i=3` (`k=4`) vs 4-word `j=4`: `~53.9` bits (Steps 19..27, non-free `W16`) vs `~72.8` bits.
