# Independent CPU audit of the tail model

`audit_tail.c` re-measures the probabilistic tail of the 36-step search independently of the GPU search code; the only thing it reuses is the generator constants. It runs on the CPU (OpenMP) and does three things:

1. Derives the per-step modular state differences of the characteristic from the Table 14 pair of ePrint 2026/1120.
2. Generates candidate pairs exactly as the GPU search does: golden first blocks, then random E4 variants, then random (W14, W15) upper candidates, enumerated by the same carry DP.
3. For each step t, checks that message 1's real state plus the characteristic difference before step t yields the characteristic difference after step t.

Build and run:

    gcc -O3 -march=native -fopenmp audit_tail.c -o audit_tail -lm
    OMP_NUM_THREADS=24 ./audit_tail 2e10 4096 2048   # pairs, uppers per lower, golden M0s

`run_2e10.txt` holds the output of the 2e10-pair run, which took about 266 s on 24 threads of a Ryzen 9 9950X3D. Summary:

| Quantity | Model | Measured |
|---|---|---|
| P(dW23 = 2^59) | 2^-11.5114 | 2^-11.5104 +/- 0.0006 |
| P(dW25 = 0 given dW23) | 9/128 = 0.07031 | 0.07041 |
| steps 24..35 hold, given dW23 and dW25 | 1 | 482,706 / 482,706 |
| W5 4-bit identity vs. sigma0 condition | exact | 0 mismatches in 1.28e8 |

Conditioned on each of steps 17 and 19 to 22, P(dW23) stays between 2^-11.50 and 2^-11.56, so dW23 shows no measurable correlation with the state conditions.

The GPU log labels two of its counters "W21" and "W23". They actually check dW23 = 2^59 and dW25 = 0, respectively.
