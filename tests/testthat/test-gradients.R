# Gradient checks guarding the torch-removal rewrite: every analytic gradient
# used on the hot path is compared against central finite differences. These
# are the tests the plan calls out as highest-value.

test_that("multinomial score matches finite differences along simplex tangents", {
  n <- 6
  fam <- multinomial_family(n, 3)
  theta <- c(0.4, 0.35, 0.25)
  sc <- score(fam, theta) # (M, 3), rows over the support

  # The multinomial score carries the simplex-tangent `-n` correction, so only
  # tangent directions (summing to zero) agree with a naive coordinate FD. Test
  # the directional derivative of log p for a few outcomes along e_1 - e_2.
  X <- support(fam)
  u <- c(1, -1, 0)
  for (i in c(1L, 5L, nrow(X))) {
    directional <- sum(sc[i, ] * u)
    f <- function(th) log_density(fam, th)[i]
    fd <- (f(theta + 1e-6 * u) - f(theta - 1e-6 * u)) / (2e-6)
    expect_equal(directional, fd, tolerance = 1e-5)
  }
})

test_that("Frank-Wolfe on-face gradient matches finite differences", {
  n <- 6
  alt <- small_alt()
  fam <- multinomial_family(n, 3)
  eng <- exact_engine(fam, as_marginal(alt, fam))
  log_P <- log(rep(1 / nrow(support(fam)), nrow(support(fam))))
  obj <- ripr:::fw_objective(eng, fam, log_P)

  f <- plurality_faces(3)[[1L]]
  V <- f@vertices
  a2v <- ripr:::alpha_from_v
  sj <- ripr:::softmax_jacobian

  value_v <- function(v) obj$value(as.vector(V %*% a2v(v)))
  grad_v <- function(v) {
    a <- a2v(v)
    as.vector(obj$grad_theta(as.vector(V %*% a)) %*% V %*% sj(a))
  }

  set.seed(1)
  v0 <- rnorm(ncol(V) - 1L)
  expect_equal(grad_v(v0), fd_grad(value_v, v0), tolerance = 1e-5)
})

test_that("EM M-step on-face gradient matches finite differences", {
  n <- 6
  alt <- small_alt()
  fam <- multinomial_family(n, 3)
  eng <- exact_engine(fam, as_marginal(alt, fam))
  M <- nrow(support(fam))
  set.seed(2)
  log_r <- log(runif(M)) # arbitrary positive responsibility weights
  obj <- ripr:::em_objective(eng, fam, log_r)

  f <- plurality_faces(3)[[2L]]
  V <- f@vertices
  a2v <- ripr:::alpha_from_v
  sj <- ripr:::softmax_jacobian

  value_v <- function(v) obj$value(as.vector(V %*% a2v(v)))
  grad_v <- function(v) {
    a <- a2v(v)
    as.vector(obj$grad_theta(as.vector(V %*% a)) %*% V %*% sj(a))
  }

  set.seed(3)
  v0 <- rnorm(ncol(V) - 1L)
  expect_equal(grad_v(v0), fd_grad(value_v, v0), tolerance = 1e-5)
})

test_that("state KL weight-gradient matches finite differences", {
  n <- 6
  alt <- small_alt()
  fam <- multinomial_family(n, 3)
  prob <- ripr_problem(fam, plurality_faces(3), as_marginal(alt, fam))
  st <- mixture_state(prob$engine, 3L)
  faces <- prob$null
  for (k in seq_len(2L)) {
    ripr:::state_add_atom(st, init_point(faces[[k]], alt_mean(alt)), k)
  }
  w0 <- c(0.4, 0.6)
  ripr:::state_set_weights(st, w0)

  g_analytic <- state_objective(st, w0)$grad
  # Gradient is wrt the raw (unnormalised) weight vector.
  loss_raw <- function(w) state_objective(st, w)$loss
  expect_equal(g_analytic, fd_grad(loss_raw, w0), tolerance = 1e-5)
})
