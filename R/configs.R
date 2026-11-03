## Flag configurations, applied cumulatively in plan order.
## family: configurations are compared within a family, against that family's baseline.
##   'default'       — the PSD solver
##   'compositional' — compositional=TRUE (R q = 0); a different estimator, so it has its
##                     own baseline and reference. Only run on clr input.
## kind: 'exact' changes must match their family baseline to ~1e-8.

CONFIGS <- list(
  baseline         = list(family='default', kind='exact', opts=list()),
  partial_eig      = list(family='default', kind='exact', opts=list(partial.eig=TRUE)),
  fast_R           = list(family='default', kind='exact', opts=list(partial.eig=TRUE, fast.R=TRUE)),
  comp_baseline    = list(family='compositional', kind='exact', opts=list(compositional=TRUE)),
  comp_fast        = list(family='compositional', kind='exact',
                          opts=list(compositional=TRUE, partial.eig=TRUE, fast.R=TRUE))
)
## Path-changing options (balance.mu, path.order, anderson) are compared at matched
## accuracy in R/accuracy_grid.R instead.

configs_for <- function(row) {
  keep <- vapply(CONFIGS, function(cf) cf$family == 'default' || isTRUE(row$clr), logical(1))
  CONFIGS[keep]
}
