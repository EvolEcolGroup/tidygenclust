- Sun 15 Nov 21:26: job started (commit 347189c (newspeed_port) with test_gt_adamixture.R from 882ad09 (coverage-speed)); sleeping until 2026-11-16 00:00
- Mon 16 Nov 00:00: waiting for the machine to be idle (5-min load < 1.0)
- Mon 16 Nov 00:10: machine idle (load 0.06 0.04 0.06)
- Mon 16 Nov 00:10: start covrO0_omp4_r1 (load 0.06 0.04 0.06)
- Mon 16 Nov 01:00: done covrO0_omp4_r1 (50 min; load after 3.41 3.97 3.91)
- Mon 16 Nov 01:01: start covrO1_omp4_r1 (load 0.78 2.95 3.55)
- Mon 16 Nov 01:58: done covrO1_omp4_r1 (57 min; load after 3.54 3.99 4.07)
- Mon 16 Nov 02:00: start covrO1_omp1_r1 (load 0.81 2.96 3.69)
- Mon 16 Nov 02:33: (by hand) covrO1_omp1 skipped in rounds 2 and 3: the tests set n_cores explicitly, so OMP_NUM_THREADS does not change them (358% CPU at OMP 1); the time goes to the other set-ups
- Mon 16 Nov 02:54: done covrO1_omp1_r1 (54 min; load after 4.32 4.86 4.99)
- Mon 16 Nov 02:58: start O0_omp4_r1 (load 0.97 2.54 3.97)
- Mon 16 Nov 03:00: (by hand) another session (tidytracks) compiled for 1-2 min on 2 cores just before this time: the run in progress may be slightly slowed; it now runs only nice'd single-core jobs
- Mon 16 Nov 03:12: done O0_omp4_r1 (13 min; load after 5.46 5.70 4.91)
- Mon 16 Nov 03:15: start O0cov_omp4_r1 (load 0.70 3.61 4.23)
- Mon 16 Nov 03:33: (by hand, Andrea's call) rounds 2 and 3 and O2_r1 cancelled: the effects (5 vs 50 min) are far above the noise; the job ends after O0cov_omp4_r1; focused single runs follow
- Mon 16 Nov 03:37: (by hand) during O0cov_omp4_r1, other sessions ran niced R jobs (about 1 core, ~03:16-03:30); indicative only
- Mon 16 Nov 04:06: done O0cov_omp4_r1 (51 min; load after 8.66 9.83 9.42)
- Mon 16 Nov 04:13: finished
- Mon 16 Nov 04:14: start covrO0_skip (load 0.66 2.33 5.73)
- Mon 16 Nov 04:19: done covrO0_skip (exit 0, 4 min; coverage 96.83%; load after 1.21 1.64 4.51)
- Mon 16 Nov 04:20: start covrO1_skip (load 0.73 1.49 4.36)
- Mon 16 Nov 04:25: done covrO1_skip (exit 0, 5 min;                       "-O1 --coverage" "-O1 --coverage -D_GLIBCXX_ASSERTIONS" 
"-O1 --coverage -D_GLIBCXX_ASSERTIONS" "-O1 --coverage -D_GLIBCXX_ASSERTIONS" 
"-O1 --coverage -D_GLIBCXX_ASSERTIONS" "-O1 --coverage -D_GLIBCXX_ASSERTIONS" 
"-O1 --coverage -D_GLIBCXX_ASSERTIONS"                       "-O1 --coverage" 
coverage 96.82%; load after 1.57 1.53 3.58)
- Mon 16 Nov 04:25: all done; other sessions told the machine is free
