## Cross-system comparison: is the speed ranking of configurations the same on
## every BLAS/LAPACK backend?
## Usage: Rscript R/compare_systems.R <results dir>
##
## Absolute times are not comparable across shared runners, so each system's
## configurations are ranked by speedup relative to that system's own baseline.

args <- commandArgs(trailingOnly=TRUE)
dir  <- if (length(args) >= 1) args[1] else "results"
files <- list.files(dir, pattern="^grid-.*\\.csv$", full.names=TRUE, recursive=TRUE)
if (!length(files)) stop("no grid-*.csv files under ", dir)
d <- do.call(rbind, lapply(files, read.csv))

med <- aggregate(seconds ~ runner + id + config, data=d, FUN=median)
base <- med[med$config == "baseline", c("runner", "id", "seconds")]
names(base)[3] <- "base_seconds"
med <- merge(med, base)
med$speedup <- med$base_seconds / med$seconds

# Speedup per configuration and system, geometric mean over cases
sp <- aggregate(speedup ~ runner + config, data=med, FUN=function(x) exp(mean(log(x))))
wide <- reshape(sp, idvar="config", timevar="runner", direction="wide")
names(wide) <- sub("^speedup\\.", "", names(wide))

# Kendall's tau between every pair of systems, on each case's configuration ranking
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
taus <- if (length(taus)) do.call(rbind, taus) else data.frame()

# Correctness: every fit must match its baseline
bad <- d[d$kind == "exact" & d$max_rel_err > 1e-8, ]

write.csv(wide, file.path(dir, "comparison-speedups.csv"), row.names=FALSE)
write.csv(taus, file.path(dir, "comparison-kendall.csv"), row.names=FALSE)
cat("Geometric-mean speedup vs each system's own baseline:\n"); print(wide, row.names=FALSE)
if (nrow(taus)) { cat("\nKendall tau of configuration rankings (1 = identical order):\n")
  print(aggregate(tau ~ a + b, data=taus, FUN=min), row.names=FALSE) }
if (nrow(bad)) { cat("\nFits that differ from baseline:\n"); print(bad); quit(status=1) }
