# SpiecEasi-bench

Benchmarks for the SpiecEasi sparse + low-rank (SLR) ADMM solver, kept outside
the package so `R CMD check` stays fast and any SpiecEasi ref can be measured
after the fact. Plan: https://claude.ai/artifact/TT13zUu7iK7MjQH4JooSJC

## What it measures

- `R/kernels.R` — per-iteration kernels (full vs partial eigen L-update, old vs
  new R-update) at p = 100…2000.
- `R/run_grid.R` — every flag configuration in `R/configs.R` on the data grid in
  `R/sim.R` (n < p, n > p, ill-conditioned, hub, penalty mode, tied spectrum).
  Each fit is checked against the same job's baseline, so no timing is reported
  for a result that differs; fits that hit `maxiter` are flagged.
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

On macOS with anaconda on `PATH`, build igraph with anaconda removed from `PATH`;
otherwise it links anaconda's libxml2 and fails to load.

## CI

`.github/workflows/bench.yml` runs the kernels and the grid (p ≤ 300 by default)
on reference LAPACK, OpenBLAS and MKL (Ubuntu), Accelerate (macOS) and reference
LAPACK (Windows), then runs the cross-system comparison. Triggers: manual with a
SpiecEasi ref, and nightly against `devel`. Not pushed yet.
