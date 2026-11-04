# SpiecEasi-bench

Benchmarks for the SpiecEasi sparse + low-rank (SLR) ADMM solver, kept outside
the package so `R CMD check` stays fast and any SpiecEasi ref can be measured
after the fact.

## What it measures

- `R/kernels.R` — per-iteration kernels (full vs partial eigen L-update, old vs
  new R-update) at p = 100…2000.
- `R/run_grid.R` — every flag configuration in `R/configs.R` on the data grid in
  `R/sim.R` (n < p, n > p, ill-conditioned, hub, penalty mode, tied spectrum).
  Each fit is checked against the same job's baseline, so no timing is reported
  for a result that differs; fits that hit `maxiter` are flagged.
- `R/accuracy_grid.R` — time and ADMM steps to reach a fixed tolerance, with each
  fit's error against a tightly converged reference. Options that change the
  iteration path (Anderson acceleration, path order) are only comparable this way.
- `R/warmstart_experiments.R` — warm starts along the lambda path, across StARS
  subsamples and across an r grid.
- `R/compare_systems.R` — combines results from several machines. The main check
  is that configurations rank the same on every BLAS/LAPACK backend (Kendall's τ),
  using speedups relative to each machine's own baseline rather than raw times.

## Running locally

```sh
export R_LIBS=$PWD/.Rlib          # isolated library; see install notes below
Rscript R/kernels.R results/kernels-local.csv
Rscript R/run_grid.R 300 3 results/grid-local.csv
Rscript R/compare_systems.R results
```

The SpiecEasi under test is whatever is installed in `R_LIBS`
(e.g. `R CMD INSTALL ../SpiecEasi` from the `perf/slr-admm` branch).

### Full suite across backends (macOS)

```sh
.backends/build.sh     # once: shim dylibs standing in for R's LAPACK/BLAS
./run_suite.sh results/suite-YYYYMMDD
Rscript R/summarize_suite.R results/suite-YYYYMMDD
```

`./rbench <ref|openblas|accelerate> script.R ...` runs a script against reference
LAPACK, OpenBLAS's LAPACK or Apple Accelerate without modifying the R install.
Jobs run one at a time, and the load average before and after each is logged in
`suite.log`: timings are only meaningful on an otherwise idle machine.

### Results so far

`results/suite-20261103/` (Apple M1, single-threaded, SpiecEasi `d8f1095`), time to
tol 1e-6 relative to the same backend's baseline:

| configuration | geomean speedup | ADMM steps |
|---|---|---|
| partial eigen L-update + one-eig R-update (identical iterates) | 1.33x | same |
| + Anderson acceleration | 2.07x | 3.6x fewer |
| + sparse-first lambda path | 2.73x (1.85-4.24x) | 4.0x fewer |

Backend choice barely matters at p <= 300 on this machine; at p = 2000 a
4-thread OpenBLAS halves the eigensolver time.

On macOS with anaconda on `PATH`, build igraph with anaconda removed from `PATH`;
otherwise it links anaconda's libxml2 and fails to load.

## CI

`.github/workflows/bench.yml` runs the kernels and the grid (p ≤ 300 by default)
on reference LAPACK, OpenBLAS and MKL (Ubuntu), Accelerate (macOS) and reference
LAPACK (Windows), then runs the cross-system comparison. Trigger it manually with the SpiecEasi
repo and ref to measure; the ref must be pushed to that repo.
