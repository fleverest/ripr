#' @include family.R random_variable.R
NULL

# Bounds on `log E_theta[X]` at given parameters: the values the
# difference-of-convex (DC) enclosure interpolates between. A family
# supplies them by implementing `log_expectation_bounds()`.

# --- The family's side --------------------------------------------------------

#' Difference-of-convex bounds on null expectations
#'
#' `log_expectation_bounds(family, x)` yields bounds on a family's
#' \eqn{\log E_\theta[X]}{log E_theta[X]} for given parameters and X. It is a
#' function of parameters theta, one per row, returning a two-column matrix of
#' lower and upper bounds at each row, with equal columns when the value is
#' known exactly.
#'
#' A new family registers this method for each kind of random variable `x` that
#' it can handle; for our Gaussian family, this is just [likelihood_ratio]s.
#' When `x` is of that class but not of a form the family can bound (a ratio
#' of distributions from another family, say), the method returns instead a
#' string explaining why. The default says that nothing is implemented:
#' "no way to evaluate or bound the expectation of a <class of x> at a point is
#' implemented for a <class of family>."
#'
#' The Gaussian family implements it in closed form for likelihood ratios with
#' a single-Gaussian denominator.
#' @param family A [parametric_family].
#' @param x A [random_variable].
#' @return A function, or a string.
#' @seealso [certify()]
#' @examples
#' fam <- gaussian_family(d = 2L)
#' x <- likelihood(fam(c(1, 0))) / likelihood(fam(c(0, 0)))
#' bounds <- log_expectation_bounds(fam, x)
#' bounds(rbind(c(0, 0), c(0.5, 0.5)))
#' @name dc_enclosure
NULL

#' @rdname dc_enclosure
#' @export
log_expectation_bounds <- new_generic(
  "log_expectation_bounds",
  c("family", "x"),
  function(family, x) S7::S7_dispatch()
)

method(log_expectation_bounds, list(parametric_family, random_variable)) <-
  function(family, x) {
    paste0(
      "no way to evaluate or bound the expectation of a `",
      class_name(x),
      "` at a point is implemented for a `",
      class_name(family),
      "`"
    )
  }
