# sha2-collide

**WORLD-FIRST PRACTICAL 36-STEP STANDARD-IV SHA-512 COLLISION ACHIEVED!** (2026-11-22)

Achieved in **`2.81 hours`** on a single consumer NVIDIA RTX 5090 (`9.04 GPU-hours` total across 1× RTX 5090 + 2× RTX PRO 6000 Blackwell Server Edition) — beating both the prior practical SHA-512 record (35-step SHA-512, CRYPTO 2026 `ePrint 2026/1080`) and the theoretical `2^{66.2}` 36-step SHA-512 complexity of `ePrint 2026/1120`.

## State of the Art Comparison

| Work | 35-step SHA-512 | 36-step SHA-256 | **36-step SHA-512 (Standard IV)** |
| :--- | :---: | :---: | :---: |
| Li et al., CRYPTO 2026 (`ePrint 2026/1080`) | Practical (`2^{54}`, RTX 5060 Ti, 14d) | Theoretical (`2^{71}`) | Theoretical (`2^{136}`) |
| Goldzweig (`ePrint 2026/2913`) | — | Practical (`2^{54.6}`, 440 GPU-h) | — |
| Li et al., `ePrint 2026/1120` | — | Practical (`2^{57}`, Table 2) | Theoretical (`2^{66.2}`, SFS only, Table 14) |
| **This project (`sha2-collide`)** | — | — | **PRACTICAL & VERIFIED (`9.04 GPU-h` total; `2.81 h` on RTX 5090!)** |

## Verified 36-Step Standard-IV SHA-512 Collision

- **Artifact**: [`results/collision_36step_sha512_std_iv.json`](results/collision_36step_sha512_std_iv.json)
- **Clean-room FIPS 180-4 verifier**: `python3 results/verify_world_first.py` (derives IV and K from primes, cross-checks 80-round SHA-512 against Python `hashlib.sha512`, and verifies both raw 2-block compression and full FIPS 180-4 Merkle-Damgård padded 256-byte messages)
- **Block 1 (`M0`, includes ASCII `"SHA51236COLLIDE!"`, `m0_seed = 2068887288786`)**:
  ```text
  000001e1b34abbd2 5348413531323336 434f4c4c49444521 0000000000000000
  0000000000000000 0000000000000000 0000000000000000 0000000000000000
  0000000000000000 0000000000000000 0000000000000000 0000000000000000
  0000000000000000 0000000000000000 0000000000000000 0000000000000000
  ```
- **Chaining Value `CV = F_36(Standard_IV, M0)`**:
  ```text
  f1e42e9f77e4cca4 a785f886ea02658d cf52ba6f18e44e11 644ba4c897d35d12
  c4cd7917344f6760 3f554360df46b76a 1a8281d8e746e746 0ef673a2b4c69e41
  ```
- **Block 2 (`M1 != M2`, differing in words `5, 6, 7, 8, 9, 13, 14`)**:
  - `M1`:
    ```text
    4f0a777b3e20699b 9084d21a7d15ffdb f3eb37ed7c49a6a9 811b8eea0f813426
    9901a1886b80a163 36982650a4377b49 e2734f1dfc05e77e 65f45efd8b33bb64
    d395fbe4be6d579f b0f8774e493b43bc 13287fa0a341d877 f0d19e5ae32f6abf
    f26d943cf8041873 06f2b17bb37eb2b7 2cccf09ec39c66dd 3284b0682fc4054b
    ```
  - `M2`:
    ```text
    4f0a777b3e20699b 9084d21a7d15ffdb f3eb37ed7c49a6a9 811b8eea0f813426
    9901a1886b80a163 3e982650a4377b49 e67b4f1dfc05e77e 6df45efd8b33bb64
    f3c5db65be6d579f f0d8764e493b43bc 13287fa0a341d877 f0d19e5ae32f6abf
    f26d943cf8041873 02fab17bb37eb2b7 24ccf09ec39c66dd 3284b0682fc4054b
    ```
- **Identical 512-bit 36-step Compression Digest `F_36(CV, M1) == F_36(CV, M2)`**:
  ```text
  7f987d4f13fdfebe99203ca3a1bf59931e4de96de6f5fbe5ba20aadb881fc9cf
  33724a7027ab508b3e25595a172fb3c86b67ab3917239813209bb9516ff380cf
  ```
- **Identical Full FIPS 180-4 Padded 36-Step SHA-512 Digest `SHA512_36(M0 || M1) == SHA512_36(M0 || M2)`**:
  ```text
  9eabd799baf71c41b426c3833361abcd9b9405a24d6ff6fa41f43e6071e437f3
  79413d838b9f4634be19868ac36dca673d3e031474cf1a596f05c4e4546e9e84
  ```

## Key Technical Breakthroughs (Zero SAT Needed)

1. **Exact Upper Freedom `N_U = 2^{25.000}` (`16×` larger than reported `2^{21}`)**: exact bit-by-bit carry DP on `ePrint 2026/1120` Table 14 proves `2^{21}` valid `W15` per `W14` × `16` valid `W14` = `33,554,432` upper probes.
2. **Lower `E4` Low-Bit Independence (`12,000×` Step-2 Acceleration)**: sweeping 29 free bits of `E4` per golden `M0` (`2^{-14}`) yields `~20.7M` valid lower solutions per golden `M0`, accelerated with a proven 100%-exact 4-bit `W5` identity + 16-way `low4` unroll.
3. **Forced Step-18 Direct Bits on `E16` (`[43:1, 51:0, 57:0]`) + `NEW15_k19_g3` L2-Resident MitM**: compacts active buckets to `61,365 / 524,288` (`11.7%`, `30.0 MB` hot L2 table, `99.0%` capture) with `<=64-U` load-balanced chunks (`65.8 k s18/s` on RTX 5090, `61.2 k s18/s` on RTX PRO 6000).
4. **Exact 4B-Sample Tail Calibration**: `31.5315 ± 0.0179` bits after Step 18 (`16.1900b` Steps 19–22 + `11.5114b` Row 23 + `3.8301b` Row 25, zero correlation), confirmed to `0.04` bits by `2.04B` live Step-18 GPU samples.

## Quick Verification

```bash
python3 results/verify_world_first.py
python3 cores/verify_table14.py
```
