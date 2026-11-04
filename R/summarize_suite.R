## Summarize a run_suite.sh output directory across backends.
## Usage: Rscript R/summarize_suite.R <suite dir>
args <- commandArgs(trailingOnly=TRUE)
dir  <- if (length(args) >= 1) args[1] else "results"
rd <- function(prefix) {
  f <- list.files(dir, pattern=paste0("^", prefix, "-.*\\.csv$"), full.names=TRUE)
  d <- do.call(rbind, lapply(f, read.csv))
  d$backend <- sub("^local-(.*)-t[0-9]+$", "\\1", d$runner)
  d$threads <- as.integer(sub(".*-t", "", d$runner))
  d
}
gm <- function(x) exp(mean(log(x)))
w <- function(d, v, row, col) {
  x <- reshape(d[, c(row, col, v)], idvar=row, timevar=col, direction="wide")
  names(x) <- sub(paste0(v, "."), "", names(x), fixed=TRUE); x
}

k <- rd("kernels")
k$what <- ifelse(is.na(k$k), k$kernel, paste0(k$kernel, " k=", k$k))
k <- k[is.na(k$k) | k$k == 2, ]
cat("== Kernel time (ms), by backend and thread count\n")
k$col <- paste0(k$backend, "-t", k$threads); k$ms <- round(1000 * k$seconds, 1)
print(w(k[order(k$p), ], "ms", c("p", "what"), "col"), row.names=FALSE)

a <- rd("accuracy")
cat("\n== Accuracy grid: total seconds over all cases, per backend\n")
print(w(aggregate(seconds ~ config + backend, a, sum), "seconds", "config", "backend"), row.names=FALSE)
base <- a[a$config %in% c("baseline", "comp_baseline"), c("backend", "id", "family", "seconds", "steps")]
names(base)[4:5] <- c("b_sec", "b_steps")
m <- merge(a, base)
m$speedup <- m$b_sec / m$seconds; m$step_ratio <- m$b_steps / m$steps
cat("\n== Geomean speedup vs same-backend baseline (time to tol 1e-6)\n")
print(w(aggregate(speedup ~ config + backend, m, gm), "speedup", "config", "backend"), row.names=FALSE)
cat("\n== Geomean step reduction (backend-independent check)\n")
print(w(aggregate(step_ratio ~ config + backend, m, gm), "step_ratio", "config", "backend"), row.names=FALSE)
cat("\n== Backend speed: geomean ref-time / backend-time, same config\n")
r <- merge(a[a$backend == "ref", c("id", "config", "seconds")], a[, c("id", "config", "backend", "seconds")], by=c("id", "config"))
r$x <- r$seconds.x / r$seconds.y
print(w(aggregate(x ~ config + backend, r, gm), "x", "config", "backend"), row.names=FALSE)
cat("\n== Steps differ across backends? (max relative step difference per config)\n")
st <- aggregate(steps ~ id + config, a, function(s) (max(s) - min(s)) / min(s))
print(aggregate(steps ~ config, st, max), row.names=FALSE)
cat("\n== Accuracy: max_rel_err vs family reference, worst over cases; capped lambdas\n")
print(aggregate(cbind(max_rel_err, lambdas_capped, ref_capped) ~ config + backend, a, max), row.names=FALSE)

cat("\n== Per case, ref backend: baseline s -> all s (speedup)\n")
rr <- m[m$backend == "ref" & m$config %in% c("all", "comp_all"), c("id", "config", "b_sec", "seconds", "speedup", "b_steps", "steps")]
rr[, 3:5] <- round(rr[, 3:5], 2); print(rr, row.names=FALSE)

g <- rd("grid")
cat("\n== Ablation grid: worst max_rel_err vs family baseline (exact configs), and geomean speedup\n")
gb <- aggregate(seconds ~ backend + id + family + config, g, median)
b2 <- gb[gb$config %in% c("baseline", "comp_baseline"), c("backend", "id", "family", "seconds")]; names(b2)[4] <- "b"
gb <- merge(gb, b2); gb$speedup <- gb$b / gb$seconds
print(merge(aggregate(max_rel_err ~ config, g, max), w(aggregate(speedup ~ config + backend, gb, gm), "speedup", "config", "backend")), row.names=FALSE)
cat("\nfits hitting maxiter in ablation grid:", sum(g$hit_maxiter, na.rm=TRUE), "of", nrow(g), "\n")
