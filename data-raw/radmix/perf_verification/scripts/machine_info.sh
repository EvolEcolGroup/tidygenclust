#!/bin/bash
# What the timings depend on: processor, cores, instruction sets, memory,
# frequency scaling, compiler, R, BLAS and OpenMP settings.
#
#   bash machine_info.sh > <results dir>/machine.txt
echo "== date"; date
echo "== host"; uname -a
echo "== CPU"; lscpu
echo "== instruction sets"
for f in sse4_2 avx avx2 fma avx512f avx512dq avx512vl avx512_fp16; do
  if grep -qw "$f" /proc/cpuinfo; then echo "  $f: yes"; else echo "  $f: no"; fi
done
echo "== NUMA"; (numactl --hardware 2>/dev/null || echo "numactl not installed")
echo "== memory"; free -g
echo "== frequency governor"
cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo "unknown"
cat /sys/devices/system/cpu/intel_pstate/no_turbo 2>/dev/null | sed 's/^/intel no_turbo: /'
cat /sys/devices/system/cpu/cpufreq/boost 2>/dev/null | sed 's/^/boost: /'
echo "== load"; uptime
echo "== perf_event_paranoid"; cat /proc/sys/kernel/perf_event_paranoid
echo "== compiler"; $(R CMD config CXX | cut -d' ' -f1) --version | head -1
echo "== R's C++ flags"
for v in CXX CXXFLAGS CXX17FLAGS CXX20FLAGS; do
  echo "  $v: $(R CMD config $v)"
done
echo "== ~/.R/Makevars (not used by the builds of build_variants.sh)"
cat ~/.R/Makevars 2>/dev/null || echo "none"
echo "== OpenMP environment"; env | grep -E '^(OMP|GOMP|OPENBLAS|MKL)_' || echo "none set"
echo "== R"
Rscript -e 'cat(R.version.string, "\n"); print(extSoftVersion()[c("BLAS")]); print(La_library());
  for (p in c("tidypopgen", "bigstatsr", "bigsnpr", "Rcpp", "RcppArmadillo"))
    cat(p, as.character(packageVersion(p)), "\n");
  if (requireNamespace("RhpcBLASctl", quietly = TRUE))
    cat("BLAS threads:", RhpcBLASctl::blas_get_num_procs(), " OpenMP threads:", RhpcBLASctl::omp_get_max_threads(), "\n")'
