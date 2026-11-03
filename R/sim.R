## Benchmark data grid for the SLR solver.
## Each case returns a data matrix with a real low-rank (latent factor) component.

suppressPackageStartupMessages(library(SpiecEasi))

# Precision matrix with a chosen condition number, built directly so the target
# is met exactly (graph2prec's binary search may not reach large condition numbers).
prec_with_condition <- function(graph, kappa) {
  p <- ncol(graph)
  A <- as.matrix(graph) * 0.3
  diag(A) <- 0
  ev <- eigen(A, symmetric=TRUE)
  lo <- min(ev$values); hi <- max(ev$values)
  # shift-and-scale so eigenvalues span [1, kappa]
  d <- 1 + (ev$values - lo) / (hi - lo) * (kappa - 1)
  ev$vectors %*% diag(d) %*% t(ev$vectors)
}

sim_case <- function(p, n, r=2, kappa=100, graph='band', tied=FALSE, seed=1,
                     clr=FALSE, depth=5000) {
  set.seed(seed)
  g <- make_graph(graph, p, p)
  Prec <- if (kappa <= 100) graph2prec(g, targetCondition=kappa) else prec_with_condition(g, kappa)
  Cor  <- cov2cor(prec2cov(Prec))
  Z <- rmvnorm(n, rep(0, p), Cor)
  U <- matrix(rnorm(n * r), n, r)
  W <- if (tied) {
    # orthogonal loadings of equal norm -> tied top eigenvalues of L
    q <- qr.Q(qr(matrix(rnorm(p * r), p, r))); q * sqrt(p) * .3
  } else matrix(rnorm(p * r, sd=.5), p, r)
  X <- Z + U %*% t(W)
  if (!clr) return(X)
  # compositional counts: log-normal abundances, multinomial sampling at a fixed depth,
  # then SpiecEasi's own clr pipeline (the input spiec.easi passes to the SLR solver)
  mu <- rnorm(p, sd=1)
  P  <- exp(sweep(X, 2, mu, '+'))
  counts <- t(apply(P, 1, function(pr) rmultinom(1, depth, pr / sum(pr))))
  SpiecEasi:::.spiec.easi.norm(counts)
}

grid_cases <- function(max_p=Inf) {
  rows <- rbind(
    data.frame(case='n<p',      p=c(100, 300, 1000), n=NA, kappa=100, graph='band', tied=FALSE, mode='rank'),
    data.frame(case='n>p',      p=c(50, 100, 300),   n=NA, kappa=100, graph='band', tied=FALSE, mode='rank'),
    data.frame(case='n>p 20x',  p=100,               n=2000, kappa=100, graph='band', tied=FALSE, mode='rank'),
    data.frame(case='ill 1e4',  p=c(100, 300),       n=NA, kappa=1e4, graph='band', tied=FALSE, mode='rank'),
    data.frame(case='ill 1e6',  p=c(100, 300),       n=NA, kappa=1e6, graph='band', tied=FALSE, mode='rank'),
    data.frame(case='ill hub',  p=c(100, 300),       n=NA, kappa=1e4, graph='hub',  tied=FALSE, mode='rank'),
    data.frame(case='penalty',  p=c(100, 300),       n=NA, kappa=100, graph='band', tied=FALSE, mode='penalty'),
    data.frame(case='tied',     p=c(100, 300),       n=NA, kappa=100, graph='band', tied=TRUE,  mode='rank')
  )
  rows$n <- ifelse(!is.na(rows$n), rows$n,
            ifelse(rows$case == 'n<p', rows$p / 2,
            ifelse(rows$case == 'n>p', rows$p * 5, rows$p * 2)))
  # ill-conditioned cases also run with n < p
  extra <- rows[rows$case %in% c('ill 1e4', 'ill 1e6'), ]
  extra$n <- extra$p / 2; extra$case <- paste(extra$case, 'n<p')
  rows <- rbind(rows, extra)
  rows$clr <- FALSE
  # clr-transformed count data, the input spiec.easi actually gives the solver
  clr_rows <- data.frame(case=c('clr n<p', 'clr n<p', 'clr n<p', 'clr n>p', 'clr n>p'),
                         p=c(20, 100, 300, 100, 300), n=c(10, 50, 150, 500, 1500),
                         kappa=100, graph='band', tied=FALSE, mode='rank', clr=TRUE)
  rows <- rbind(rows, clr_rows)
  rows$id <- sprintf('%s_p%d_n%d', gsub('[^a-z0-9]+', '-', rows$case), rows$p, rows$n)
  rows[rows$p <= max_p, ]
}
