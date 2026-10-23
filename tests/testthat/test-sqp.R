# Properties of the general SQP over x >= 0 and its QP kernel.
#
# The QP kernel must agree with quadprog; the SQP must solve a quadratic in
# one iteration and a non-quadratic convex problem with active bounds to a
# KKT residual near rounding, and report a stop short of tol.

test_that("qp_nonneg agrees with quadprog", {
  # Our QP solve is very similar to using `quadprog::solve.QP`, but quadprog
  # requires a p.d. Hessian. Our hessians will often be nearly singular (when
  # atoms are close to one another).
  set.seed(2)
  m <- matrix(stats::rnorm(40), 8L, 5L)
  h <- crossprod(m) + diag(1e-3, 5L)
  c <- stats::rnorm(5L)
  ours <- ripr:::qp_nonneg(h, c, numeric(5L))
  ref <- quadprog::solve.QP(h, -c, diag(5L), numeric(5L))$solution
  expect_equal(ours, pmax(ref, 0), tolerance = 1e-8)
})

test_that("sqp_nonneg solves a non-negative quadratic in one iteration", {
  # A sanity check for sequential quadratic programming when the objective is
  # genuinely quadratic; the "sequential" part is not necessary.
  set.seed(3)
  m <- matrix(stats::rnorm(40), 8L, 5L)
  h <- crossprod(m) + diag(1e-3, 5L)
  c <- stats::rnorm(5L)
  fit <- ripr:::sqp_nonneg(
    rep(1, 5L),
    fn = \(x) sum(x * (h %*% x)) / 2 + sum(c * x),
    gr = \(x) as.vector(h %*% x) + c,
    hess = \(x) h
  )
  expect_identical(fit$iterations, 1L)
  expect_lte(fit$kkt, 1e-12)
  expect_equal(fit$par, ripr:::qp_nonneg(h, c, numeric(5L)), tolerance = 1e-10)
})

test_that("sqp_nonneg reaches the bound-constrained optimum of a convex function", {
  # sum(exp(x) - b x) over x >= 0 is minimised at pmax(log(b), 0).
  b <- c(0.2, 0.9, 1, 2, 5)
  fit <- ripr:::sqp_nonneg(
    rep(1, 5L),
    fn = \(x) sum(exp(x) - b * x),
    gr = \(x) exp(x) - b,
    hess = \(x) diag(exp(x), length(x))
  )
  expect_lte(fit$kkt, 1e-12)
  expect_equal(fit$par, pmax(log(b), 0), tolerance = 1e-12)
  expect_identical(fit$par[1:2], c(0, 0)) # strictly active bounds are exact
  short <- ripr:::sqp_nonneg(
    rep(1, 5L),
    fn = \(x) sum(exp(x) - b * x),
    gr = \(x) exp(x) - b,
    hess = \(x) diag(exp(x), length(x)),
    max_iter = 1L
  )
  expect_gt(short$kkt, 1e-12)
})
