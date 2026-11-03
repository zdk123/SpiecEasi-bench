## Flag configurations, applied cumulatively in plan order.
## kind: 'exact' changes must match the baseline to ~1e-8; 'path' changes are
## judged against a tightly converged reference instead.

CONFIGS <- list(
  baseline    = list(kind='exact', opts=list()),
  partial_eig = list(kind='exact', opts=list(partial.eig=TRUE)),
  fast_R      = list(kind='exact', opts=list(partial.eig=TRUE, fast.R=TRUE))
)
## Path-changing options (balance.mu, path.order, anderson) are compared at matched
## accuracy in R/accuracy_grid.R instead.
