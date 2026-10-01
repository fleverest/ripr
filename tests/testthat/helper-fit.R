# Brute-force references for the fitting tests, through the exported API only.
#
# Each enumerates the family's sample space and evaluates densities with
# `log_density()`, so it shares nothing with the quadrature engine or the step
# layer. A test comparing a verb's output with one of these checks the
# package's own reductions rather than restating them, and survives any
# refactor of the internals that leaves the answers alone. Only finite sample
# spaces are supported, which is every multinomial fixture in the suite.

# `log q(x)` and `log p_W(x)` at every outcome of the state's sample space.
exact_densities <- function(state, W = state@mixing) {
  fam <- state@null@family
  x <- enumerate_space(fam@sample_space)
  list(
    family = fam,
    x = x,
    log_q = log_density(state@alternative, x),
    log_p = log_density(fam(W), x)
  )
}

# KL(Q || P_W), with the `0 log 0 = 0` convention.
exact_kl <- function(state, W = state@mixing) {
  e <- exact_densities(state, W)
  q <- exp(e$log_q)
  sum(ifelse(q > 0, q * (e$log_q - e$log_p), 0))
}

# G(theta) = E_Q[p_theta / P_W] at each row of `theta`: the linear oracle, and
# the gradient of KL in the weights up to sign.
exact_g <- function(state, theta, W = state@mixing) {
  e <- exact_densities(state, W)
  theta <- if (is.matrix(theta)) theta else matrix(theta, nrow = 1L)
  vapply(
    seq_len(nrow(theta)),
    function(i) {
      sum(exp(e$log_q + log_density(e$family(theta[i, ]), e$x) - e$log_p))
    },
    numeric(1)
  )
}

# Is every atom of the state inside the part it is filed under?
in_own_part <- function(state, tol = 1e-5) {
  all(vapply(
    seq_along(state@part),
    function(j) {
      contains(
        parts(state@null@region)[[state@part[j]]],
        atoms(state@mixing)[j, ],
        tol = tol
      )
    },
    logical(1)
  ))
}

# Steps taken by each verb, read off the trace's `phase` column.
verb_counts <- function(state) {
  vapply(
    c(fw = "fw", lb = "lb", em = "em", weight = "weight"),
    \(p) sum(state@trace$phase == p),
    integer(1)
  )
}
