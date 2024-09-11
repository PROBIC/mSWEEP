# mSWEEP GPU comparisons

The following abundance estimations were performed starting with efaec-1_1.aln.gz (2.9 GBs) and efaec-1_2.aln.gz (2.8 GBs) pseudoalignment files which can be obtained by following this mGEMS [tutorial](https://github.com/PROBIC/mGEMS/blob/master/docs/TUTORIAL.md).

## Table 1: algorithm comparisons across HPC platforms in terms of time, iterations, and memory usage.

**Notes:**
- On Turso, A100 GPUs were used, although other GPUs are also possible to use but have less memory and will most likely run slower.
- On LUMI, older versions of LibTorch and ROCm had to be used, most likely affecting the resulting times.
- Since the emgpu algorthm with the default tolenrace of 1e-6 took all 5000 iterations in this case (rare), some results from running the algorithms with a higher tolerance of 1e-3 are shown. This tolerance still seems to provide nearly identical results but in a faster time (see Table 2 for comparison of results).
- Time was acquired from the time taken to execute [this line](https://github.com/Piketulus/mSWEEP-gpu/blob/4ca2acd510c9dfb5f0fed1d3cc3e383a3a7e8572/src/mSWEEP.cpp#L440).

| **Platform** | **Algorithm**          | **Tolerance** | **Time to Estimate Abundances (seconds)** | **Iterations** | **Max Memory Used (GB)** |
|--------------|------------------------|---------------|-------------------------------------------|----------------|--------------------------|
| Turso        | rcgcpu (8 CPUs)         | 1.00E-06      | 1856                                      | 205            | 22.7                     |
| Turso        | rcgcpu (32 CPUs)        | 1.00E-06      | 634                                       | 215            | 23.3                     |
| Turso        | rcgcpu (80 CPUs)        | 1.00E-06      | 485                                       | 215            | 24.4                     |
| Turso        | rcggpu                  | 1.00E-06      | 43                                        | 220            | 27.9 (on GPU)            |
| Turso        | rcggpu                  | 1.00E-03      | 33                                        | 155            | 27.9 (on GPU)            |
| Turso        | emgpu (double)          | 1.00E-06      | 258                                       | 5000           | 14 (on GPU)              |
| Turso        | emgpu (double)          | 1.00E-03      | 143                                       | 2605           | 14 (on GPU)              |
| Turso        | emgpu (float)           | 1.00E-06      | 19                                        | 335            | 7 (on GPU)               |
| LUMI         | rcggpu                  | 1.00E-06      | 103                                       | 225            | 27.9 (on GPU)            |
| LUMI         | emgpu (double)          | 1.00E-06      | 392                                       | 5000           | 14 (on GPU)              |
| LUMI         | emgpu (float)           | 1.00E-06      | 57                                        | 300            | 7 (on GPU)               |

**note** emgpu has lower numerical precision and the results will differ from the rcg algorithms.
