# sha2-collide

**A practical collision for 36-step SHA-512 with the standard IV** (found 2026-11-22).

The search ran on one RTX 5090 together with two RTX PRO 6000 (Blackwell). The RTX 5090 found the pair after 2.81 h; the three GPUs spent about 8.9 GPU-hours in total. At the measured RTX 5090 throughput this is about 2^49.2 evaluations of the 36-step compression function.

This is a practical instantiation of an existing attack, not a new one. It uses the 36-step differential characteristic, the message-difference pattern and the memory-efficient two-block framework of Li et al. ([ePrint 2026/1120](https://eprint.iacr.org/2026/1120), [ePrint 2026/1080](https://eprint.iacr.org/2026/1080)). It also reuses the internal-state solution of their published 36-step semi-free-start pair (2026/1120, Table 14). The paper estimates the attack at 2^66.2.

Write-up: [`paper/sha512-36step-collision.pdf`](paper/sha512-36step-collision.pdf), LaTeX source in `paper/`.

## Context: collision attacks on step-reduced SHA-512 (standard IV)

| Steps | Status | Time (compression calls) | Reference |
| :---: | :--- | :--- | :--- |
| 35 | practical | ≈ 2^54 (14 days on an RTX 5060 Ti) | Li, Liu, Wang, Shi, CRYPTO 2026 ([2026/1080](https://eprint.iacr.org/2026/1080)) |
| 36 | theoretical | 2^136 | [2026/1080](https://eprint.iacr.org/2026/1080) |
| 36 | theoretical | 2^66.2 | Li, Zhang, Li, Liu, Qian, Zhu ([2026/1120](https://eprint.iacr.org/2026/1120)) |
| **36** | **practical** | **≈ 2^49.2 measured (≈ 2^49.9 expected)** | **this repository** |
| 37 / 38 / 39 | theoretical | 2^87 / 2^125.4 / 2^178 | [2026/1120](https://eprint.iacr.org/2026/1120) |

For SHA-256, the latest practical standard-IV collisions are 36 steps ([2026/1120](https://eprint.iacr.org/2026/1120)) and 38 steps (Goldzweig, [2026/2913](https://eprint.iacr.org/2026/2913), found in 53 GPU-hours).

## The colliding pair

- Record: [`results/collision_36step_sha512_std_iv.json`](results/collision_36step_sha512_std_iv.json)
- Block 1 (`M0`). Word 0 is the seed `2068887288786`; words 1–2 are the ASCII string `"SHA51236COLLIDE!"`:
  ```text
  000001e1b34abbd2 5348413531323336 434f4c4c49444521 0000000000000000
  0000000000000000 0000000000000000 0000000000000000 0000000000000000
  0000000000000000 0000000000000000 0000000000000000 0000000000000000
  0000000000000000 0000000000000000 0000000000000000 0000000000000000
  ```
- Chaining value `CV = F_36(standard IV, M0)`:
  ```text
  f1e42e9f77e4cca4 a785f886ea02658d cf52ba6f18e44e11 644ba4c897d35d12
  c4cd7917344f6760 3f554360df46b76a 1a8281d8e746e746 0ef673a2b4c69e41
  ```
- Block 2. `M1` and `M2` differ in words 5, 6, 7, 8, 9, 13, 14:
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
- Both blocks give the same compression output, `F_36(CV, M1) == F_36(CV, M2)`:
  ```text
  7f987d4f13fdfebe99203ca3a1bf59931e4de96de6f5fbe5ba20aadb881fc9cf
  33724a7027ab508b3e25595a172fb3c86b67ab3917239813209bb9516ff380cf
  ```
- With standard padding, the 256-byte messages `M0||M1` and `M0||M2` have the same 36-step SHA-512 digest:
  ```text
  9eabd799baf71c41b426c3833361abcd9b9405a24d6ff6fa41f43e6071e437f3
  79413d838b9f4634be19868ac36dca673d3e031474cf1a596f05c4e4546e9e84
  ```

## How the attack was made practical

None of these techniques is new; together they bring the expected cost from 2^66.2 down to about 2^49.9 compression-function equivalents.

1. **Degrees of freedom in (W14, W15).** Four unconditioned bits of A14 give 16 values of W14. For each, a bit-level carry DP yields 2^21 values of W15, for 2^25 upper candidates in total. The paper counts 2^21.
2. **Re-solving (A0, E4) on the fly.** A first block is kept if it meets the 14 conditions on E3, W6 and W7. For each kept block, the 29 unconditioned bits of E4 are swept: about 20.7M of the 2^29 variants (3.85%) meet all lower conditions. This makes the first-block cost negligible; the `2^63` term of the paper's estimate disappears.
3. **Meet-in-the-middle on E16, resident in L2.** E16 = ℓ + h, where ℓ depends only on the lower part and h only on the upper part. Lowers and uppers are matched on a 19-bit key: 61,365 non-empty buckets, a 30 MB hot table that stays in L2, and 99.0% capture. Each matched pair then goes through 15 direct bits of E16 plus the step-17/18 filters. Throughput is 65.8k step-18 survivors/s on the RTX 5090 and 61.2k/s on each RTX PRO 6000.
4. **Tail cost model.** After step 18 the remaining conditions cost 2^31.53: 16.19 bits for steps 19–22, 11.51 bits for δW23, and 3.83 bits for δW25. Over the whole run, the observed rates were 2^-16.20 for steps 19–22, and 9 hits for δW23 against 9.3 expected.

The GPU log names two counters "W21" and "W23". They actually check δW23 = 2^59 and δW25 = 0.

## Independent checks

- `paper/verify_appendix.py`: a from-scratch SHA-512 whose constants are derived from primes. It is tested against `hashlib` for the full 80 steps, then verifies the pair at 36 steps, both raw and with padding.
- `audit/`: a CPU-only audit of the tail model on 2·10^10 generated pairs. P(δW23) measured 2^-11.5104 against 2^-11.5114 predicted. P(δW25 | δW23) measured 0.07041 against 9/128. Steps 24–35 held in 482,706 of 482,706 cases. See `audit/README.md`.
- `bench/`: throughput of the 36-step compression function on the RTX 5090, 2.15·10^10/s, with self-tests.

## Repository layout

```text
src/gpu/      CUDA search (sha512_36_t14_gpu.cu) and earlier benchmarks
src/ref/      reference reduced-round SHA-512 (Python)
src/sat/      SAT probes of the upper-part rigidity
cores/        Table 14 constants, verifiers, tail calibration counter
results/      the colliding pair (JSON) and a verifier
audit/        independent CPU audit of the tail model
bench/        36-step compression throughput benchmark
paper/        preprint (LaTeX + PDF) and appendix verifier
papers/       links to the reference papers (PDFs not redistributed)
archive/      earlier audit notes
PROGRESS.md   chronological log
```

## Quick verification

```bash
python3 results/verify_world_first.py
python3 paper/verify_appendix.py results/collision_36step_sha512_std_iv.json
python3 cores/verify_table14.py
```

## AI-use designation

Predominantly AI-generated work, directed by Junhong Huang. Google Gemini 4 Neon, running in the Antigravity IDE, wrote the code and ran the experiments. Anthropic Claude Opus 5.5 contributed the initial exploration, the independent verification, and the write-up.
