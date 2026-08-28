# Properties of R/random_variable.R.
#
# A random variable maps elements of a sample space to numbers. Two things are
# worth guarding: that it refuses input which is not an element of its own
# sample space, since silently reshaping the wrong thing returns a number rather
# than a complaint; and that arithmetic composes without anyone having to trust
# it, since `Q / P_star / bound` is how an e-variable candidate gets built.
#
# Nothing here claims validity. That is the point -- an object with no place to
# record a claim cannot record a false one.

fixture <- function(n = 8L, k = 3L) {
  fam <- multinomial_family(n_trials = n, k = k)
  theta <- seq(k, 1) / sum(seq(k, 1))
  list(
    family = fam,
    space = fam@sample_space,
    Q = induced_distribution(fam, point_mixing(theta)),
    P = induced_distribution(fam, point_mixing(rep(1 / k, k)))
  )
}

# Some testthat expect_error calls kind of break with S7:
#
# Use this instead of `expect_error()` when it breaks. `expect_error()` labels
# the frames it captures, which forces the `x` argument of the S7 method that
# threw (a promise still under evaluation at that moment). Re-entering it fails,
# and the failure is reported against whichever function was forcing it, hiding
# the real error behind a formatting one.
error_message <- function(expr) {
  tryCatch(
    {
      force(expr)
      NA_character_
    },
    error = conditionMessage
  )
}

# --- Evaluation ---------------------------------------------------------------

test_that("a random variable takes one element or many", {
  # One element as a length-d vector gives one number; n of them as an (n, d)
  # matrix gives n. The wrapped function only ever sees the matrix form.
  f <- fixture()
  first <- random_variable(function(x) x[, 1L], f$space)
  expect_equal(first(c(4, 2, 2)), 4)
  expect_equal(first(rbind(c(4, 2, 2), c(0, 8, 0))), c(4, 0))
  expect_length(first(c(4, 2, 2)), 1L)
})

test_that("the wrapped function is never handed unchecked input", {
  # It is written for an (n, d) matrix and may assume that shape.
  f <- fixture()
  seen <- NULL
  probe <- random_variable(
    function(x) {
      seen <<- x
      rep(0, nrow(x))
    },
    f$space
  )
  invisible(probe(c(4, 2, 2)))
  expect_true(is.matrix(seen))
  expect_identical(dim(seen), c(1L, 3L))
})

test_that("a likelihood is the density, and a ratio is a quotient of two", {
  # There is no ratio class: `Q / P` is ordinary arithmetic on two likelihoods.
  f <- fixture()
  L <- likelihood(f$Q)
  x <- rbind(c(4, 2, 2), c(8, 0, 0))
  expect_equal(L(x), exp(log_density(f$Q, x)))

  R <- likelihood(f$Q) / likelihood(f$P)
  expect_equal(R(x), exp(log_density(f$Q, x) - log_density(f$P, x)))
})

test_that("a ratio integrates to 1 against its own denominator", {
  # `E_P[Q/P] = 1` is the one value known in closed form without integrating
  # anything, so it catches a mis-signed or mis-weighted density.
  f <- fixture()
  R <- likelihood(f$Q) / likelihood(f$P)
  outcomes <- enumerate_space(f$family@sample_space)
  expect_equal(
    sum(exp(log_density(f$P, outcomes)) * R(outcomes)),
    1,
    tolerance = rounding_tol(1)
  )
})

test_that("infinite values are allowed", {
  # A ratio is genuinely infinite where its denominator vanishes and its
  # numerator does not, and that is reachable: a mixture with an atom at a
  # simplex vertex gives zero probability to almost every outcome.
  f <- fixture()
  vertex <- induced_distribution(f$family, point_mixing(c(1, 0, 0)))
  R <- likelihood(f$Q) / likelihood(vertex)
  expect_true(is.infinite(R(c(4, 2, 2))))
  expect_true(R(c(8, 0, 0)) < Inf)
})

# --- Refusing what is not an outcome ------------------------------------------

test_that("input must have the right shape", {
  f <- fixture()
  R <- likelihood(f$Q)
  expect_match(error_message(R(c(4, 4))), "length-3 vector")
  expect_match(error_message(R(c(2, 2, 2, 2))), "length-3 vector")
  expect_match(error_message(R(matrix(1, nrow = 2L, ncol = 4L))), "3 columns")
})

test_that("input must be numeric and present", {
  f <- fixture()
  R <- likelihood(f$Q)
  expect_match(error_message(R(c("a", "b", "c"))), "outcomes must be numeric")
  expect_match(error_message(R(c(4, NA, 4))), "outcomes must not be missing")
  expect_match(error_message(R(c(4, NaN, 4))), "outcomes must not be missing")
})

test_that("count outcomes must be counts summing to the total", {
  # The shape check alone would let through a vector of the right length that
  # is not a point of the sample space at all.
  f <- fixture()
  R <- likelihood(f$Q)
  expect_match(error_message(R(c(4, 2, 1))), "summing to 8", fixed = TRUE)
  expect_match(error_message(R(c(9, -1, 0))), "non-negative whole numbers")
  expect_match(error_message(R(c(4.5, 2, 1.5))), "non-negative whole numbers")
  expect_silent(R(c(8, 0, 0)))
})

test_that("Gaussian outcomes may be anything finite", {
  # The sample space is all of R^d, so only finiteness is added: an infinite
  # outcome has zero density under every parameter, making a ratio there 0 / 0.
  fam <- gaussian_family(dim = 2L)
  L <- likelihood(induced_distribution(fam, point_mixing(c(0.5, 0.5))))
  expect_silent(L(c(-3, 40)))
  expect_match(error_message(L(c(Inf, 0))), "finite")
  expect_match(error_message(L(1)), "length-2 vector")
})

# --- Arithmetic ---------------------------------------------------------------

test_that("each operator does what it says, either way round", {
  f <- fixture()
  X <- likelihood(f$Q)
  Y <- likelihood(f$P)
  x <- rbind(c(4, 2, 2), c(8, 0, 0), c(0, 4, 4))

  expect_equal((X + Y)(x), X(x) + Y(x))
  expect_equal((X - Y)(x), X(x) - Y(x))
  expect_equal((X * Y)(x), X(x) * Y(x))
  expect_equal((X / Y)(x), X(x) / Y(x))

  expect_equal((X + 2)(x), X(x) + 2)
  expect_equal((2 + X)(x), 2 + X(x))
  expect_equal((X - 2)(x), X(x) - 2)
  expect_equal((2 - X)(x), 2 - X(x))
  expect_equal((X * 3)(x), X(x) * 3)
  expect_equal((3 * X)(x), 3 * X(x))
  expect_equal((X / 3)(x), X(x) / 3)
  expect_equal((3 / X)(x), 3 / X(x))
})

test_that("results are random variables, and compose", {
  f <- fixture()
  R <- likelihood(f$Q) / likelihood(f$P)
  x <- rbind(c(4, 2, 2), c(8, 0, 0))
  expect_true(S7_inherits(R / 2, random_variable))
  expect_equal(((R / 2 + 1) * 3)(x), (R(x) / 2 + 1) * 3)
})

test_that("a derived variable still checks its input", {
  # Validation is not something the outermost variable does on everyone's
  # behalf; each operand checks, so any variable is sound on its own.
  f <- fixture()
  R <- likelihood(f$Q) / 2
  expect_error(R(c(4, 2, 1)), "summing to 8")
})

test_that("variables on different sample spaces cannot be combined", {
  # Checked when written, not when evaluated, so the error names the line that
  # made the mistake.
  a <- likelihood(fixture(k = 3L)$Q)
  b <- likelihood(fixture(k = 4L)$Q)
  expect_error(a + b, "different sample spaces")
  expect_error(a / b, "different sample spaces")
})

test_that("only a single number may be combined with a variable", {
  f <- fixture()
  X <- likelihood(f$Q)
  expect_error(X * c(1, 2), "single number")
})

# --- Printing -----------------------------------------------------------------

test_that("a leaf prints its label", {
  f <- fixture()
  Q <- f$Q
  expect_match(rv_expression(likelihood(Q)), "^Q$")
  expect_match(rv_expression(likelihood(Q, label = "alt")), "^alt$")
  expect_match(
    rv_expression(random_variable(function(x) x[, 1L], f$space)),
    "rv"
  )
})

test_that("an unlabelled call records how it was written, warts and all", {
  # `substitute()` captures syntax rather than identity, so a call behind a
  # helper labels itself with the helper's argument name. Cosmetic, and why
  # `label` exists.
  f <- fixture()
  helper <- function(d) likelihood(d)
  expect_match(rv_expression(helper(f$Q)), "^d$")
})

test_that("a label must be a single string or NULL", {
  f <- fixture()
  expect_error(likelihood(f$Q, label = c("a", "b")), "single string")
  expect_error(likelihood(f$Q, label = NA_character_), "single string")
  expect_error(likelihood(f$Q, label = 3), "single string")
})

test_that("arithmetic prints as the expression that built it", {
  f <- fixture()
  Q <- f$Q
  P <- f$P
  R <- likelihood(Q) / likelihood(P)
  expect_equal(rv_expression(R), "Q / P")
  expect_equal(rv_expression(R / 1.5), "Q / P / 1.5")
  expect_equal(rv_expression(R + R), "Q / P + Q / P")
})

test_that("brackets appear only where they change the reading", {
  # `*` and `/` bind tighter than `+` and `-`; neither `-` nor `/` associates,
  # so an equally strong right operand of either needs bracketing.
  f <- fixture()
  Q <- f$Q
  P <- f$P
  R <- likelihood(Q) / likelihood(P)
  expect_equal(rv_expression(2 / R), "2 / (Q / P)")
  expect_equal(rv_expression((R + 1) * 3), "(Q / P + 1) * 3")
  expect_equal(rv_expression(R - R / 2), "Q / P - Q / P / 2")
  expect_equal(rv_expression(R / R / 2), "Q / P / (Q / P) / 2")
})

test_that("the expression tree records the operands", {
  f <- fixture()
  R <- likelihood(f$Q) / 2
  expect_identical(R@op, "/")
  expect_length(R@operands, 2L)
  expect_true(S7_inherits(R@operands[[1L]], random_variable))
  expect_equal(R@operands[[2L]], 2)
  expect_true(is.na(likelihood(f$Q)@op))
})

test_that("a long label is shortened", {
  f <- fixture()
  long <- likelihood(induced_distribution(f$family, point_mixing(c(0.5, 0.3, 0.2))))
  expect_true(nchar(rv_expression(long)) <= 24L)
  expect_match(rv_expression(long), "\\.\\.\\.$")
})

test_that("printing returns the variable invisibly", {
  f <- fixture()
  R <- likelihood(f$Q)
  expect_silent(invisible(capture.output(out <- print(R))))
  expect_true(S7_inherits(out, random_variable))
})

test_that("format() gives the expression, and does not error", {
  family <- multinomial_family(n_trials = 4L, k = 3L)
  Q <- induced_distribution(family, point_mixing(c(0.5, 0.3, 0.2)))
  x <- likelihood(Q, label = "Q")
  p <- likelihood(Q, label = "P*")

  expect_identical(format(x), "Q")
  expect_identical(format(x / p), "Q / P*")
  expect_identical(format(x / 4.27), "Q / 4.27")

  # A scalar, not the deparsed function, and not multi-line.
  expect_type(format(x), "character")
  expect_length(format(x), 1L)

  # `print()` is the expression plus a class banner, so the two must agree.
  expect_true(any(grepl(
    format(x / p),
    capture.output(print(x / p)),
    fixed = TRUE
  )))
})
