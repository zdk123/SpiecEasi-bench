## Time to a fixed accuracy, per configuration.
## Usage: Rscript R/accuracy_grid.R [max_p] [tol] [out.csv]
##
## Configurations that change the iteration path stop at different points for the same
## tol, so speed is only meaningful together with accuracy. Each case gets a reference
## solution (accelerated run at tol 1e-10, cross-checked against plain ADMM in the package
## tests); every configuration is then run at `tol` and reports time, ADMM steps and its
## error against the reference.

source("R/sim.R")

args  <- commandArgs(trailingOnly=TRUE)
max_p <- if (length(args) >= 1) as.numeric(args[1]) else 300
tol   <- if (length(args) >= 2) as.numeric(args[2]) else 1e-6
out   <- if (length(args) >= 3) args[3] else "results/accuracy.csv"

if (requireNamespace("RhpcBLASctl", quietly=TRUE) && nzchar(Sys.getenv("BENCH_THREADS")))
  RhpcBLASctl::blas_set_num_threads(as.integer(Sys.getenv("BENCH_THREADS")))

## Each family is a separate estimator with its own reference solution; compositional
## configurations run on clr input only.
ALL <- list(partial.eig=TRUE, fast.R=TRUE, anderson=10, path.order='sparse_first')
FAMILIES <- list(default=list(), compositional=list(compositional=TRUE))
CONFIGS <- list(
  baseline      = list(family='default', opts=list()),
  exact_kernels = list(family='default', opts=list(partial.eig=TRUE, fast.R=TRUE)),
  anderson      = list(family='default', opts=list(anderson=10)),
  all           = list(family='default', opts=ALL),
  comp_baseline = list(family='compositional', opts=list(compositional=TRUE)),
  comp_anderson = list(family='compositional', opts=list(compositional=TRUE, anderson=10)),
  comp_all      = list(family='compositional', opts=c(list(compositional=TRUE), ALL))
)
MAXITER <- 5e4

rel_err <- function(a, b) norm(as.matrix(a) - as.matrix(b), 'F') / max(1, norm(as.matrix(b), 'F'))

fit <- function(X, row, opts) {
  a <- list(data=X, cor=TRUE, nlambda=5, lambda.min.ratio=1e-2,
            lambda.max=pulsar::getMaxCov(cor(X)), opts=opts)
  if (row$mode == 'rank') a$r <- 2 else a$beta <- .1
  t0 <- proc.time()[['elapsed']]
  f <- do.call(sparseLowRankiCov, a)
  f$elapsed <- proc.time()[['elapsed']] - t0
  f
}

cases <- grid_cases(max_p)
cases <- cases[!duplicated(cases$case) | cases$p > 100, ]
if (nzchar(Sys.getenv("BENCH_CASES"))) cases <- cases[grepl(Sys.getenv("BENCH_CASES"), cases$id), ]  # regex on id
res <- list()
for (ci in seq_len(nrow(cases))) {
  row <- cases[ci, ]
  X <- sim_case(row$p, row$n, kappa=row$kappa, graph=row$graph, tied=row$tied, clr=row$clr)
  fams <- if (isTRUE(row$clr)) names(FAMILIES) else 'default'
  refs <- lapply(setNames(fams, fams), function(fm)
    fit(X, row, c(FAMILIES[[fm]], list(anderson=10, partial.eig=TRUE, fast.R=TRUE,
                                       tol=1e-10, maxiter=1e5))))
  for (cfg in names(CONFIGS)) {
    fm <- CONFIGS[[cfg]]$family
    if (!fm %in% fams) next
    ref <- refs[[fm]]
    f <- fit(X, row, c(CONFIGS[[cfg]]$opts, list(tol=tol, maxiter=MAXITER)))
    err <- max(mapply(rel_err, f$icov, ref$icov), mapply(rel_err, f$resid, ref$resid))
    capped <- sum(f$iter >= MAXITER)
    ref_capped <- sum(ref$iter >= 1e5)
    res[[length(res) + 1]] <- data.frame(id=row$id, case=row$case, p=row$p, n=row$n, clr=row$clr,
      config=cfg, family=fm, tol=tol, seconds=f$elapsed, steps=sum(f$iter),
      lambdas_capped=capped, ref_capped=ref_capped, max_rel_err=err)
    cat(sprintf("%-28s %-14s %8.2fs  steps %7d  capped %d/%d  err %.1e%s\n", row$id, cfg, f$elapsed,
                sum(f$iter), capped, length(f$iter), err, if (ref_capped) "  (reference capped!)" else ""))
  }
}
res <- do.call(rbind, res)
res$lapack <- La_library(); res$blas <- extSoftVersion()[['BLAS']]
res$runner <- Sys.getenv("BENCH_RUNNER", "local")
res$spieceasi_ref <- Sys.getenv("BENCH_SPIECEASI_REF", NA)
dir.create(dirname(out), showWarnings=FALSE, recursive=TRUE)
write.csv(res, out, row.names=FALSE)
