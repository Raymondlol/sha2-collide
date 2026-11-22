# sha2-collide

Practical collision attacks on reduced-round SHA-2 — currently targeting the **world-first practical 36-step Standard-IV SHA-512 collision** on a single NVIDIA RTX 5090 (32 GiB GDDR7, 96 MB L2 cache).

## Current Target & State of the Art

| Work | 35-step SHA-512 | 36-step SHA-256 | **36-step SHA-512 (Standard IV)** |
| :--- | :---: | :---: | :---: |
| Li et al., CRYPTO 2026 (`ePrint 2026/1080`) | Practical (`2^{54}`, RTX 5060 Ti, 14d) | Theoretical (`2^{71}`) | Theoretical (`2^{136}`) |
| Goldzweig (`ePrint 2026/2913`) | — | Practical (`2^{54.6}`, 440 GPU-h) | — |
| Li et al., `ePrint 2026/1120` | — | Practical (`2^{57}`, Table 2) | **Theoretical (`2^{66.2}`)** — SFS only (Table 14) |
| **This project (`sha2-collide`)** | — | — | **~`2^{49.6}` after L2 MitM (`31.18b` after Step 18) -> ~6–12 h projected on 1× RTX 5090** |

## Verified Breakthroughs on `ePrint 2026/1120` Table 14 (Zero SAT Needed)

All numbers below have been independently bit-for-bit verified against `src/ref/sha512_ref.py` (`python3 cores/verify_table14.py`):

1. **Exact Upper Freedom `N_U = 2^{25.000}` (`16×` larger than reported `2^{21}`)**
   - Exact bit-by-bit carry DP on Table 14 (`t=14, 15`) finds `2^{21.000}` valid `W15` per `W14` (5 `A15/S0` constraints become redundant under `E15[41..60]=0`) × `16` valid `W14` (`A14` free bits `[28, 40, 62, 63]`) = **`33,554,432 = 2^{25.000}`** valid `(W14, W15)` upper probes.
   - `1,600 / 1,600` randomly sampled pairs passed full 36-step trace verification.

2. **Lower `E4` Low-Bit Independence (`12,000×` Step-2 Acceleration)**
   - Instead of paying `2^{-18}` per random `M0` with fixed `E4`, find a "golden `M0`" passing `E3 + W7 + W6` (`2^{-14}`, ~`113k/s` on RTX 5090).
   - Varying `E4` over its **29 free bits** (`[0..3, 5..17, 20..24, 27..29, 34..36, 63]`) keeps `E3 = 100.00%`, `W7 = 100.00%`, `W8+Ch7 = 91.0%`, `W6 = 80.4%`, leaving only `W5` (`6.19%`).
   - Net pass rate = **`4.62%` (`2^{-4.44}`)** -> **`2^{24.56}` (~24.7M) valid lower solutions per single golden `M0`** without recomputing `SHA-512(M0)`.

3. **Real Modular Tail + Carry-Blocker L2-Resident MitM (`k = 19`, `64 MB` in RTX 5090's `96 MB` L2)**
   - Real 64-bit SHA-512 modular tail (`Steps 17..35`): `15.0b (Step 17) + 16.0b (Step 18) + 31.18b (Steps 19..35)` = **`62.18` bits** (vs strict `66.0` bits in `2026/1120`).
   - Pre-sorting the static `2^{25}` upper `(W14, W15)` by a **`19`-bit key** (`6` carry-blocked direct bits + `6` Majority-carry-predicted direct bits + `1` Majority-predicted `A16[59]`, `524,288` buckets × ~`64` `U`/bucket) and filling `M = 2^{23}` lower `L`s per batch (`64 MB`, `16 × 8B = 128B` = **1 exact L2 cache line per bucket**):
     - Shared-memory tiled `16 L × 64 U` (`1,024` pairs per L2 load) completely removes L2 bandwidth bottlenecks.
     - Projected end-to-end wall-clock time on 1× RTX 5090: **`~6–12 hours`** (exact measured ms/batch to be confirmed in Step B's 30s `gpu-run` benchmark).

## Repository Layout

```text
sha2-collide/
├── README.md                              # Project overview & verified architecture
├── PROGRESS.md                            # Chronological experiment & discovery log
├── cores/
│   ├── table14_36step_sha512.json         # Verified Table 14 core constants & key specs
│   └── verify_table14.py                  # Standalone bit-for-bit verifier (run anytime)
├── src/
│   ├── ref/sha512_ref.py                  # Reference reduced-round SHA-512
│   ├── gpu/bench_phase1_gpu.cu            # Verified RTX 5090 Phase-1 bench (1.85 G M0/s)
│   ├── gpu/sha512_36_gpu.cu               # Verified RTX 5090 23-bit L2 bench (126.09 G pi/s)
│   └── sat/                               # Step-A CaDiCaL upper-core rigidity probes
├── papers/                                # Reference ePrint PDFs (2026/1080, 2026/1120, 2026/2913)
├── archive/AUDIT_SUMMARY_2026-11-22.md    # Earlier overnight SAT/GPU audit summary
├── runs/                                  # Checkpoints & GPU run logs (gitignored)
└── tools -> ../train-domains/sha512/tools # CaDiCaL & Kissat builds (symlink, gitignored)
```

## Quick Verification

```bash
python3 cores/verify_table14.py
```

## Next Step (Step B)

1. Implement `src/gpu/sha512_36_t14_gpu.cu` (static `k=19` sorted `2^{25}` upper `U` array + golden-`M0`/29-bit-`E4` Phase 1 + `16×64` shared-memory tiled Phase 2 + full Steps 18–35 verifier).
2. Run a 30-second queued `gpu-run` self-test (injecting Table 14 as a sentinel collision to prove 100% end-to-end detection) and measure exact `ms/batch` + ETA on RTX 5090.
3. Launch the checkpointed full production run via `GPU_RUN_LABEL=sha512-36-t14 gpu-run ...`.
