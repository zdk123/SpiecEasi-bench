# sparseLowRankiCov: loglik uses det() of the sparse part only

**Draft issue for zdk123/SpiecEasi — not filed.**

In `R/SparseLowRankICov.R`, the per-lambda log-likelihood is

```r
R <- icov[[i]] #- resid[[i]]
loglik[[i]] <- log(Matrix::det(R[z,z])) - sum(Matrix::diag(R[z,z] %*% SigmaO[z,z])) - (p-q)
```

Two problems:

1. **Overflow/underflow.** `det()` of a p×p matrix over- or underflows for moderate p
   (product of p eigenvalues), giving `Inf`, `0` and then `log(0) = -Inf`. Use
   `Matrix::determinant(R[z,z], logarithm=TRUE)$modulus`.
2. **Which precision matrix.** `icov` is the sparse part S. The SLR model's precision
   is R = S − L; the commented-out `- resid[[i]]` suggests this was considered. S alone
   need not be positive definite, so its log-determinant can be `NaN`. This loglik
   feeds EBIC model selection, so the choice changes selected networks.

Decision needed before a fix: use log-det of S − L (model precision; selection results
change), or keep S and only fix the numerics.
