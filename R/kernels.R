## Kernel microbenchmarks: the per-iteration building blocks of the SLR ADMM.
## Usage: Rscript R/kernels.R [out.csv]
## Kernels missing from the installed SpiecEasi are skipped, so the script runs
## against any ref on the perf/slr-admm branch.

suppressPackageStartupMessages(library(SpiecEasi))
args <- commandArgs(trailingOnly=TRUE)
out  <- if (length(args) >= 1) args[1] else "results/kernels.csv"

if (requireNamespace("RhpcBLASctl", quietly=TRUE) && nzchar(Sys.getenv("BENCH_THREADS")))
  RhpcBLASctl::blas_set_num_threads(as.integer(Sys.getenv("BENCH_THREADS")))

ns <- asNamespace("SpiecEasi")
has <- function(f) exists(f, envir=ns, inherits=FALSE)
kern <- function(f) get(f, envir=ns)

time_it <- function(expr_fun, reps) {
  expr_fun()  # warm-up
  t <- vapply(seq_len(reps), function(i) system.time(expr_fun())[['elapsed']], numeric(1))
  median(t)
}

ps <- as.numeric(strsplit(Sys.getenv("KERNEL_P", "100,300,1000,2000"), ",")[[1]])
res <- list()
for (p in ps) {
  set.seed(p)
  X <- matrix(rnorm(p * p), p); M <- (X + t(X)) / 2
  reps <- if (p >= 1000) 3 else 10
  add <- function(kernel, k, sec) {
    res[[length(res) + 1]] <<- data.frame(p=p, kernel=kernel, k=k, seconds=sec)
    cat(sprintf("p=%5d  %-24s k=%-3s %.4fs\n", p, kernel, k, sec))
  }
  for (k in c(1, 2, 5, 20)) {
    if (has(".eigPSDThresh2"))       add("L full eig",    k, time_it(function() kern(".eigPSDThresh2")(M, k, 0), reps))
    if (has(".eigPSDThreshPartial")) add("L partial eig", k, time_it(function() kern(".eigPSDThreshPartial")(M, k, 0), reps))
  }
  if (has(".Rupdate_old")) add("R update old", NA, time_it(function() kern(".Rupdate_old")(M, 1), reps))
  if (has(".Rupdate_new")) add("R update new", NA, time_it(function() kern(".Rupdate_new")(M, 1), reps))
}
res <- do.call(rbind, res)
res$lapack <- La_library(); res$blas <- extSoftVersion()[['BLAS']]
res$runner <- Sys.getenv("BENCH_RUNNER", "local")
res$spieceasi_ref <- Sys.getenv("BENCH_SPIECEASI_REF", NA)
dir.create(dirname(out), showWarnings=FALSE, recursive=TRUE)
write.csv(res, out, row.names=FALSE)
