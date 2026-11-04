#!/bin/bash
# Full local timing suite across LAPACK/BLAS backends, run sequentially so jobs never
# compete for cores. Usage: ./run_suite.sh [outdir]
# Stages (cheapest first, so partial runs are still useful):
#   1. kernel microbenchmarks, 1 thread and 4 threads (the M1's performance cores)
#   2. matched-accuracy grid (time to tol 1e-6 vs a 1e-10 reference), 1 thread
#   3. bit-compatible ablation grid, 1 thread
# Load average is logged before and after every job; check it before trusting a timing.
out=${1:-results/suite-$(date +%Y%m%d)}
mkdir -p "$out"
export BENCH_SPIECEASI_REF=$(cat .Rlib/SpiecEasi/BENCH_REF 2>/dev/null)
log="$out/suite.log"
backends="ref openblas accelerate"
job() {  # job <label> <backend> <threads> <script> <args...>
  local label=$1 b=$2 t=$3; shift 3
  echo "$(date '+%F %T') START $label $b t$t  load: $(sysctl -n vm.loadavg)" >> "$log"
  BENCH_THREADS=$t ./rbench "$b" "$@" > "$out/$label-$b-t$t.out" 2>&1
  echo "$(date '+%F %T') END   $label $b t$t  status $?  load: $(sysctl -n vm.loadavg)" >> "$log"
}
echo "SpiecEasi $BENCH_SPIECEASI_REF on $(sysctl -n machdep.cpu.brand_string)" >> "$log"
for t in 1 4; do for b in $backends; do
  job kernels "$b" $t R/kernels.R "$out/kernels-local-$b-t$t.csv"
done; done
for b in $backends; do job accuracy "$b" 1 R/accuracy_grid.R 300 1e-6 "$out/accuracy-local-$b-t1.csv"; done
for b in $backends; do job grid "$b" 1 R/run_grid.R 300 3 "$out/grid-local-$b-t1.csv"; done
echo "$(date '+%F %T') DONE" >> "$log"
