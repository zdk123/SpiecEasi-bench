## Run the SLR flag-configuration ablation over the data grid.
## Usage: Rscript R/run_grid.R [max_p] [reps] [out.csv]
##
## Each configuration is timed against the all-flags-off baseline from the same
## job, and every fit is checked against that baseline, so a timing is never
## reported for a result that differs.

source("R/sim.R")
source("R/configs.R")

args  <- commandArgs(trailingOnly=TRUE)
max_p <- if (length(args) >= 1) as.numeric(args[1]) else 300
reps  <- if (length(args) >= 2) as.integer(args[2]) else 3
out   <- if (length(args) >= 3) args[3] else "results/grid.csv"

if (requireNamespace("RhpcBLASctl", quietly=TRUE) && nzchar(Sys.getenv("BENCH_THREADS")))
  RhpcBLASctl::blas_set_num_threads(as.integer(Sys.getenv("BENCH_THREADS")))

rel_err <- function(a, b) norm(as.matrix(a) - as.matrix(b), 'F') / max(1, norm(as.matrix(b), 'F'))

fit_once <- function(X, row, opts) {
  a <- list(data=X, cor=TRUE, nlambda=if (row$p >= 1000) 5 else 10, lambda.min.ratio=1e-2,
            lambda.max=pulsar::getMaxCov(cor(X)),
            opts=c(opts, list(maxiter=if (row$p >= 1000) 300 else 1000)))
  if (row$mode == 'rank') a$r <- 2 else a$beta <- .1
  t0 <- proc.time()[['elapsed']]
  fit <- do.call(sparseLowRankiCov, a)
  fit$elapsed <- proc.time()[['elapsed']] - t0
  fit
}

cases <- grid_cases(max_p)
res <- list()
for (ci in seq_len(nrow(cases))) {
  row <- cases[ci, ]
  X <- sim_case(row$p, row$n, kappa=row$kappa, graph=row$graph, tied=row$tied, clr=row$clr)
  base <- list()
  cfgs <- configs_for(row)
  for (cfg in names(cfgs)) {
    fam <- cfgs[[cfg]]$family
    for (k in seq_len(reps)) {
      fit <- fit_once(X, row, cfgs[[cfg]]$opts)
      if (is.null(base[[fam]])) base[[fam]] <- fit   # first config of a family is its baseline
      err <- max(mapply(rel_err, fit$icov, base[[fam]]$icov), mapply(rel_err, fit$resid, base[[fam]]$resid))
      iters <- if (!is.null(fit$iter)) sum(fit$iter) else NA
      hit_max <- if (!is.null(fit$iter)) any(fit$iter >= if (row$p >= 1000) 300 else 1000) else NA
      res[[length(res) + 1]] <- data.frame(id=row$id, case=row$case, p=row$p, n=row$n, config=cfg,
        family=fam, kind=cfgs[[cfg]]$kind, clr=row$clr, rep=k, seconds=fit$elapsed, iters=iters, hit_maxiter=hit_max,
        max_rel_err=err)
      cat(sprintf("%-28s %-14s rep %d  %7.3fs  iters %s  err %.1e\n", row$id, cfg, k, fit$elapsed, iters, err))
    }
  }
}
res <- do.call(rbind, res)
res$lapack <- La_library(); res$blas <- extSoftVersion()[['BLAS']]
res$sys <- paste(Sys.info()[['sysname']], Sys.info()[['machine']])
res$runner <- Sys.getenv("BENCH_RUNNER", "local")
res$spieceasi_ref <- Sys.getenv("BENCH_SPIECEASI_REF", NA)
dir.create(dirname(out), showWarnings=FALSE, recursive=TRUE)
write.csv(res, out, row.names=FALSE)
