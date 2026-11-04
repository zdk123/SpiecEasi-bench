## Cross-system comparison: is the speed ranking of configurations the same on
## every BLAS/LAPACK backend?
## Usage: Rscript R/compare_systems.R <results dir>
##
## Absolute times are not comparable across shared runners, so each system's
## configurations are ranked by speedup relative to that system's own baseline.
## Covers both the bit-compatible ablation (grid-*.csv) and the matched-accuracy
## runs (accuracy-*.csv).

args <- commandArgs(trailingOnly=TRUE)
dir  <- if (length(args) >= 1) args[1] else "results"
failed <- FALSE

compare <- function(prefix) {
  files <- list.files(dir, pattern=paste0("^", prefix, "-.*\\.csv$"), full.names=TRUE, recursive=TRUE)
  if (!length(files)) { cat("no", prefix, "results\n"); return(invisible()) }
  d <- do.call(rbind, lapply(files, read.csv))
  cat("\n==", prefix, "==\n")

  # speedups are within a family (default vs compositional are different estimators)
  if (is.null(d$family)) d$family <- "default"
  med <- aggregate(seconds ~ runner + id + family + config, data=d, FUN=median)
  base <- med[med$config %in% c("baseline", "comp_baseline"), c("runner", "id", "family", "seconds")]
  names(base)[4] <- "base_seconds"
  med <- merge(med, base)
  med$speedup <- med$base_seconds / med$seconds
  # summarise over the cases every system finished (some run a smaller grid)
  common <- Reduce(intersect, split(med$id, med$runner))
  med <- med[med$id %in% common, ]

  sp <- aggregate(speedup ~ runner + config, data=med, FUN=function(x) exp(mean(log(x))))
  wide <- reshape(sp, idvar="config", timevar="runner", direction="wide")
  names(wide) <- sub("speedup.", "", names(wide), fixed=TRUE)
  cat("Geometric-mean speedup vs each system's own baseline:\n"); print(wide, row.names=FALSE)
  write.csv(wide, file.path(dir, paste0("comparison-", prefix, "-speedups.csv")), row.names=FALSE)

  runners <- unique(med$runner)
  taus <- list()
  for (i in seq_along(runners)) for (j in seq_along(runners)) if (i < j) {
    for (cid in unique(med$id)) {
      a <- med[med$runner == runners[i] & med$id == cid, c("config", "seconds")]
      b <- med[med$runner == runners[j] & med$id == cid, c("config", "seconds")]
      ab <- merge(a, b, by="config")
      if (nrow(ab) >= 3)
        taus[[length(taus) + 1]] <- data.frame(a=runners[i], b=runners[j], id=cid,
          tau=suppressWarnings(cor(ab$seconds.x, ab$seconds.y, method="kendall")))
    }
  }
  if (length(taus)) {
    taus <- do.call(rbind, taus)
    write.csv(taus, file.path(dir, paste0("comparison-", prefix, "-kendall.csv")), row.names=FALSE)
    cat("\nKendall tau of configuration rankings, worst case per system pair (1 = identical order):\n")
    print(aggregate(tau ~ a + b, data=taus, FUN=min), row.names=FALSE)
  }

  # Correctness
  # grid: bit-compatible configs must match the baseline; accuracy: no config may be
  # more than 10x less accurate than the baseline at the same tol
  bad <- if (prefix == "grid") d[d$kind == "exact" & d$max_rel_err > 1e-8, ] else {
    be <- d[d$config %in% c("baseline", "comp_baseline"), c("runner", "id", "family", "max_rel_err")]
    names(be)[4] <- "base_err"
    dd <- merge(d, be)
    dd[dd$max_rel_err > 10 * dd$base_err + 1e-4, ]
  }
  if (nrow(bad)) { cat("\nFits failing their accuracy check:\n"); print(bad); failed <<- TRUE }
}

compare("grid")
compare("accuracy")
if (failed) quit(status=1)
