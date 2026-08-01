# Properties of R/numerics.R.
#
# These hold by definition, so a failure here is a bug in a helper and nothing
# downstream can be trusted. Since the log-sum-exp reductions delegate to
# matrixStats, the assertions on them are a contract with a dependency: if
# matrixStats ever changes its -Inf or empty-input conventions, these say so.

# --- NaN handling -------------------------------------------------------------

test_that("nan_to_zero short-circuits without altering clean input", {
  clean <- matrix(rnorm(12), nrow = 3)
  expect_identical(nan_to_zero(clean), clean)

  dirty <- clean
  dirty[2, 2] <- NaN
  expect_equal(nan_to_zero(dirty)[2, 2], 0)
  expect_equal(nan_to_zero(dirty)[-5], clean[-5])
})

test_that("nan_to_zero leaves NA alone", {
  # anyNA() fires for NA too, but only NaN should be rewritten: an NA here means
  # a bug upstream and must not be silently zeroed.
  out <- nan_to_zero(c(1, NA, NaN))
  expect_true(is.na(out[2]))
  expect_false(is.nan(out[2]))
  expect_equal(out[3], 0)
})

test_that("nan_to_neginf rewrites NaN and leaves NA alone", {
  out <- nan_to_neginf(c(1, NA, NaN))
  expect_true(is.na(out[2]))
  expect_false(is.nan(out[2]))
  expect_equal(out[3], -Inf)
})

# --- Column-wise arithmetic ---------------------------------------------------

test_that("add_by_col and div_by_col agree with sweep", {
  m <- matrix(rnorm(24), nrow = 6)
  w <- c(1.5, -0.5, 2, 0.25)
  expect_equal(add_by_col(m, w), sweep(m, 2L, w, "+"))
  expect_equal(div_by_col(m, w), sweep(m, 2L, w, "/"))
})

test_that("recycling reproduces a row-wise sweep", {
  # Why the row-wise case needs no helper: R recycles column-major.
  m <- matrix(rnorm(24), nrow = 6)
  v <- rnorm(6)
  expect_equal(m - v, sweep(m, 1L, v, "-"))
  expect_equal(m + v, sweep(m, 1L, v, "+"))
})

test_that("column-wise helpers handle a single-column matrix", {
  m <- matrix(rnorm(5), ncol = 1L)
  expect_equal(add_by_col(m, 2), m + 2)
  expect_equal(div_by_col(m, 2), m / 2)
})

# --- Log-sum-exp reductions ---------------------------------------------------

test_that("logsumexp_vec matches the naive form when it is safe to compute", {
  v <- c(-1.5, 0.2, 3.7, -8)
  expect_equal(logsumexp_vec(v), log(sum(exp(v))))
})

test_that("logsumexp_vec is -Inf iff every entry is -Inf", {
  expect_equal(logsumexp_vec(rep(-Inf, 4)), -Inf)
  expect_true(is.finite(logsumexp_vec(c(-Inf, -Inf, 2))))
})

test_that("logsumexp_vec survives arguments that would overflow exp()", {
  v <- c(1000, 1001, 999)
  expect_equal(logsumexp_vec(v), 1001 + log(exp(-1) + 1 + exp(-2)))
  expect_true(is.finite(logsumexp_vec(v)))
})

test_that("row_logsumexp matches the naive form where it is safe to compute", {
  m <- matrix(rnorm(60), nrow = 12)
  expect_equal(row_logsumexp(m), log(rowSums(exp(m))))
})

test_that("row and column logsumexp agree under transposition", {
  m <- matrix(rnorm(35), nrow = 5)
  expect_equal(row_logsumexp(m), col_logsumexp(t(m)))
})

test_that("row_logsumexp survives arguments that would overflow exp()", {
  m <- rbind(c(1000, 1001, 999), c(-1000, -1001, -999))
  out <- row_logsumexp(m)
  expect_true(all(is.finite(out)))
  expect_equal(out[1], 1001 + log(exp(-1) + 1 + exp(-2)))
})

test_that("row_logsumexp handles all--Inf rows without NaN", {
  m <- rbind(c(-Inf, -Inf), c(0, log(3)))
  out <- row_logsumexp(m)
  expect_equal(out[1], -Inf)
  expect_equal(out[2], log(4))
  expect_false(anyNA(out))
})

test_that("row_logsumexp of a zero-column matrix is -Inf per row", {
  # log(sum over nothing) = log(0). An empty subnull contributes a zero-column
  # atom matrix, so this path is reached in ordinary use.
  out <- row_logsumexp(matrix(numeric(0), nrow = 3L, ncol = 0L))
  expect_equal(out, rep(-Inf, 3))
  expect_false(anyNA(out))
})

test_that("row_logsumexp returns a plain vector", {
  expect_null(dim(row_logsumexp(matrix(rnorm(12), nrow = 4L))))
  expect_length(row_logsumexp(matrix(rnorm(4), ncol = 1L)), 4L)
})

test_that("logsumexp_weighted computes log(sum w_i exp(v_i))", {
  v <- c(0.3, -2, 1.1)
  w <- c(0.5, 0.2, 0.3)
  expect_equal(logsumexp_weighted(v, log(w)), log(sum(w * exp(v))))
})

test_that("logsumexp_weighted is -Inf when every weight is zero", {
  expect_equal(logsumexp_weighted(c(1, 2, 3), rep(-Inf, 3)), -Inf)
})

# --- -Inf-safe matrix multiplication -----------------------------------------

test_that("matmul_0_ninf treats 0 * -Inf as 0 where %*% gives NaN", {
  a <- matrix(c(0, 1), nrow = 1)
  b <- matrix(c(-Inf, 2), ncol = 1)
  expect_true(is.nan(drop(a %*% b)))
  expect_equal(drop(matmul_0_ninf(a, b)), 2)
})

test_that("matmul_0_ninf still propagates -Inf against a positive weight", {
  a <- matrix(c(1, 1), nrow = 1)
  b <- matrix(c(-Inf, 2), ncol = 1)
  expect_equal(drop(matmul_0_ninf(a, b)), -Inf)
})

test_that("matmul_0_ninf agrees with %*% when no -Inf is present", {
  a <- matrix(abs(rnorm(12)), nrow = 3)
  b <- matrix(rnorm(20), nrow = 4)
  expect_equal(matmul_0_ninf(a, b), a %*% b)
})

test_that("the matmul_0_ninf guard does not change the answer", {
  # The guarded fast path and the always-two-matmul form must agree on both
  # sides of the branch.
  two_matmul <- function(a, b) {
    neg_inf <- is.infinite(b) & b < 0
    b_safe <- b
    b_safe[neg_inf] <- 0
    out <- a %*% b_safe
    out[(a != 0) %*% neg_inf > 0] <- -Inf
    out
  }
  a <- matrix(rpois(30, 3), nrow = 6)
  b_interior <- matrix(rnorm(25), nrow = 5)
  b_boundary <- b_interior
  b_boundary[2L, 3L] <- -Inf
  b_boundary[5L, 1L] <- -Inf

  expect_equal(matmul_0_ninf(a, b_interior), two_matmul(a, b_interior))
  expect_equal(matmul_0_ninf(a, b_boundary), two_matmul(a, b_boundary))
})

# --- Coercion and predicates --------------------------------------------------

test_that("as_outcome_matrix treats a bare vector as a single outcome", {
  expect_equal(dim(as_outcome_matrix(c(3, 7))), c(1L, 2L))
  expect_equal(dim(as_outcome_matrix(matrix(1:6, nrow = 3))), c(3L, 2L))
})

test_that("is_prob_vector accepts the simplex and rejects everything else", {
  expect_true(is_prob_vector(c(0.25, 0.25, 0.5)))
  expect_true(is_prob_vector(c(1, 0, 0)))
  expect_false(is_prob_vector(c(0.5, 0.6)))
  expect_false(is_prob_vector(c(1.5, -0.5)))
  expect_false(is_prob_vector(c(0.5, NA)))
})
