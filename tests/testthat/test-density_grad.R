# Properties of compile_density_grad(), the weighted density gradient behind
# the linear and Li--Barron oracles.
#
# Checked against finite differences, in particular including the a face of the
# simplex, where the density is 0 and the score is Inf.

# sum_x exp(log_m(x)) p_theta(x) by enumeration for finite differences.
weighted_density <- function(fam, x, theta, log_m) {
  sum(exp(log_m + compile_loglik(fam, x)(theta)))
}

test_that("the multinomial density gradient matches finite differences", {
  fam <- multinomial_family(n_trials = 6L, k = 3L)
  x <- enumerate_space(fam@sample_space)
  set.seed(1)
  log_m <- log(stats::runif(nrow(x)))
  dg <- ripr:::compile_density_grad(fam, x)
  h <- 1e-6

  for (theta in list(c(0.3, 0.3, 0.4), c(0.5, 0.5, 0))) {
    got <- dg(theta, log_m)
    # One-sided into the simplex, so the face point stays a valid argument.
    fd <- vapply(
      seq_along(theta),
      function(j) {
        step <- replace(numeric(3), j, h)
        (weighted_density(fam, x, theta + step, log_m) -
          weighted_density(fam, x, theta, log_m)) /
          h
      },
      numeric(1)
    )
    expect_true(all(is.finite(got)))
    expect_equal(as.vector(got), fd, tolerance = 1e-4)
  }
})

test_that("the default density gradient agrees with the multinomial one inside the simplex", {
  fam <- multinomial_family(n_trials = 5L, k = 3L)
  x <- enumerate_space(fam@sample_space)
  log_m <- -seq_len(nrow(x)) / 10
  theta_mat <- rbind(c(0.2, 0.5, 0.3), c(0.6, 0.1, 0.3))
  default <- S7::method(ripr:::compile_density_grad, ripr::parametric_family)
  expect_equal(
    default(fam, x)(theta_mat, log_m),
    ripr:::compile_density_grad(fam, x)(theta_mat, log_m)
  )
})

test_that("the linear oracle's gradient is finite on a face of the simplex", {
  fam <- multinomial_family(n_trials = 8L, k = 3L)
  null <- null_model(
    fam,
    simplex_region(
      vertices = rbind(c(0.5, 0.5, 0), c(0.5, 0, 0.5), c(0, 0.5, 0.5))
    )
  )
  st <- ripr_init(fam(c(0.6, 0.2, 0.2)), null)
  ld <- ripr:::compile_engine(st@engine)
  log_p <- ripr:::log_p_at_nodes(st, ld)
  grad <- ripr:::linear_oracle(st, log_p, ld)$grad(c(0.5, 0.5, 0))
  expect_true(all(is.finite(grad)))
})
