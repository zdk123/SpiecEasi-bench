## Warm-start experiments for the SLR solver.
## Usage: Rscript R/warmstart_experiments.R [out_prefix]
##
## 1. threshold.warm: apply the L threshold on the first step of each warm-started fit.
## 2. anderson.carry: carry the Anderson damping counter across the lambda path.
##    (1 and 2 are compared at matched accuracy: steps to tol, error vs a tight reference.)
## 3. StARS subsample warm starts: start each subsample's lambda path from the full-data
##    solution (all lambdas, or the first lambda only) instead of cold. Compared on ADMM
##    steps AND on what StARS selects (pulsar's own stars.stability), at the current
##    default tol and at a tight tol.

source("R/sim.R")
args <- commandArgs(trailingOnly=TRUE)
prefix <- if (length(args) >= 1) args[1] else "results/warmstart"
rel_err <- function(a, b) norm(as.matrix(a) - as.matrix(b), 'F') / max(1, norm(as.matrix(b), 'F'))

fit <- function(X, opts, r=2, nlambda=10) {
  t0 <- proc.time()[['elapsed']]
  f <- sparseLowRankiCov(X, cor=TRUE, r=r, nlambda=nlambda, lambda.min.ratio=1e-2,
                         lambda.max=pulsar::getMaxCov(cor(X)), opts=opts)
  f$elapsed <- proc.time()[['elapsed']] - t0
  f
}

## ---- 1 & 2: warm-start details along the path ---------------------------------
ACC <- list(anderson=10, path.order='sparse_first', partial.eig=TRUE, fast.R=TRUE)
cfg12 <- list(
  plain            = list(),
  plain_thrwarm    = list(threshold.warm=TRUE),
  aa               = ACC,
  aa_thrwarm       = c(ACC, list(threshold.warm=TRUE)),
  aa_carry         = c(ACC, list(anderson.carry=TRUE)),
  aa_thrwarm_carry = c(ACC, list(threshold.warm=TRUE, anderson.carry=TRUE)))
cases12 <- list(np20=c(20, 60), pn30=c(30, 15), np100=c(100, 500), pn100=c(100, 50))
res12 <- list()
for (cn in names(cases12)) {
  pn <- cases12[[cn]]; X <- sim_case(pn[1], pn[2])
  ref <- fit(X, c(ACC, list(tol=1e-11, maxiter=1e5)))
  for (k in names(cfg12)) for (tol in c(1e-5, 1e-7)) {
    f <- fit(X, c(cfg12[[k]], list(tol=tol, maxiter=1e5)))
    err <- max(mapply(rel_err, f$icov, ref$icov), mapply(rel_err, f$resid, ref$resid))
    res12[[length(res12) + 1]] <- data.frame(case=cn, config=k, tol=tol, steps=sum(f$iter),
                                             seconds=f$elapsed, err=err)
    cat(sprintf("[1/2] %-6s %-17s tol %.0e  steps %6d  %6.2fs  err %.1e\n", cn, k, tol, sum(f$iter), f$elapsed, err))
  }
}
write.csv(do.call(rbind, res12), paste0(prefix, "-path.csv"), row.names=FALSE)

## ---- 3: StARS subsample warm starts ------------------------------------------
stars_run <- function(X, opts, mode, rep.num=20, thresh=0.05, seed=10010) {
  n <- nrow(X); p <- ncol(X)
  ratio <- if (n > 144) 10 * sqrt(n) / n else 0.8      # pulsar's default subsample ratio
  set.seed(seed)
  idx <- replicate(rep.num, sample(seq_len(n), floor(n * ratio)), simplify=FALSE)
  t0 <- proc.time()[['elapsed']]
  full <- fit(X, c(opts, list(return.states=(mode != 'cold'))))
  init <- if (mode == 'cold') list() else
          list(init.states=full$states, init.first.only=(mode == 'warm_first'))
  subs <- lapply(idx, function(ix) fit(X[ix, ], c(opts, init)))
  secs <- proc.time()[['elapsed']] - t0
  nl <- length(full$lambda)
  premerge <- lapply(seq_len(nl), function(i) lapply(subs, function(s) s$path[[i]]))
  st <- pulsar:::stars.stability(premerge, thresh, rep.num, p)
  list(opt=st$opt.index, summary=st$summary, steps=sum(full$iter) + sum(sapply(subs, function(s) sum(s$iter))),
       seconds=secs)
}

cases3 <- list(np50=c(50, 200), pn50=c(50, 30), np100=c(100, 300))
modes <- c('cold', 'warm_all', 'warm_first')
tols3 <- list(default=list(tol=1e-3), tight=c(ACC, list(tol=1e-7, maxiter=1e5)))
res3 <- list()
for (cn in names(cases3)) {
  pn <- cases3[[cn]]; X <- sim_case(pn[1], pn[2])
  for (tn in names(tols3)) {
    cold <- NULL
    for (md in modes) {
      s <- stars_run(X, tols3[[tn]], md)
      if (md == 'cold') cold <- s
      res3[[length(res3) + 1]] <- data.frame(case=cn, tol=tn, mode=md, steps=s$steps, seconds=s$seconds,
        opt_index=s$opt, max_summary_diff=max(abs(s$summary - cold$summary)))
      cat(sprintf("[3]   %-6s %-8s %-11s steps %7d  %7.1fs  StARS opt %2d  max |dStability| %.3f\n",
                  cn, tn, md, s$steps, s$seconds, s$opt, max(abs(s$summary - cold$summary))))
    }
  }
}
write.csv(do.call(rbind, res3), paste0(prefix, "-stars.csv"), row.names=FALSE)

## ---- 4: warm starts across an r grid ---------------------------------------------
## The user's workflow: fit r = 1..R, pick r by eBIC at each StARS-selected lambda.
## Each r's path starts cold, or from the previous r's per-lambda solution.
## Also records how far each r's fit is from its own tight reference at the current
## default tol, since eBIC differences between ranks are only meaningful if each fit
## is close to its optimum.
cases4 <- list(np50=c(50, 200), pn50=c(50, 30), np100=c(100, 300))
res4 <- list()
for (cn in names(cases4)) {
  pn <- cases4[[cn]]; X <- sim_case(pn[1], pn[2], r=3)
  for (tn in names(tols3)) {
    prev <- NULL
    for (r in 1:5) {
      ref  <- fit(X, c(ACC, list(tol=1e-11, maxiter=1e5)), r=r)
      cold <- fit(X, tols3[[tn]], r=r)
      warm <- if (is.null(prev)) cold else fit(X, c(tols3[[tn]], list(init.states=prev$states)), r=r)
      prev <- fit(X, c(tols3[[tn]], list(return.states=TRUE)), r=r)
      e <- function(f) max(mapply(rel_err, f$icov, ref$icov), mapply(rel_err, f$resid, ref$resid))
      res4[[length(res4) + 1]] <- data.frame(case=cn, tol=tn, r=r, cold_steps=sum(cold$iter),
        warm_steps=sum(warm$iter), cold_err=e(cold), warm_err=e(warm),
        cold_loglik=paste(signif(cold$loglik, 6), collapse=";"),
        ref_loglik=paste(signif(ref$loglik, 6), collapse=";"))
      cat(sprintf("[4]   %-6s %-8s r=%d  steps cold %6d warm %6d  err cold %.1e warm %.1e\n",
                  cn, tn, r, sum(cold$iter), sum(warm$iter), e(cold), e(warm)))
    }
  }
}
write.csv(do.call(rbind, res4), paste0(prefix, "-rgrid.csv"), row.names=FALSE)
