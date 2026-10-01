# Gauss quadrature rules: 1D rules quadrature rules from statmod, multivariate
# rules via tensor products of those.

#' Gauss-Hermite rule for the weight `exp(-t^2)`; weights sum to `sqrt(pi)`
#' @keywords internal
#' @noRd
gauss_hermite <- function(n) {
  rule <- statmod::gauss.quad(n, kind = "hermite")
  list(nodes = rule$nodes, weights = rule$weights)
}


#' Gauss-Jacobi rule on `[0, 1]` for the weight `(1-u)^a u^b`
#'
#' statmod's rule is on `[-1, 1]` for `(1-t)^a (1+t)^b`; `u = (1 + t) / 2`
#' scales the weights by `2^-(a+b+1)`, so they sum to `B(a+1, b+1)`. The nodes
#' lie strictly inside `(0, 1)`.
#' @keywords internal
#' @noRd
gauss_jacobi_01 <- function(n, a, b = 0) {
  rule <- statmod::gauss.quad(n, kind = "jacobi", alpha = a, beta = b)
  list(nodes = (1 + rule$nodes) / 2, weights = rule$weights / 2^(a + b + 1))
}


#' Tensor product of one-dimensional rules
#'
#' Returns `list(nodes, log_w)`, where `nodes` is an (M, d) matrix,
#' and `log_w` is a length-M vector.
#' @keywords internal
#' @noRd
tensor_rule <- function(rules) {
  idx <- as.matrix(expand.grid(lapply(rules, function(r) seq_along(r$nodes))))
  nodes <- matrix(0, nrow = nrow(idx), ncol = length(rules))
  log_w <- numeric(nrow(idx))
  for (i in seq_along(rules)) {
    nodes[, i] <- rules[[i]]$nodes[idx[, i]]
    log_w <- log_w + log(rules[[i]]$weights[idx[, i]])
  }
  list(nodes = nodes, log_w = log_w)
}
