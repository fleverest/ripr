#' @include parameter_space.R family.R
NULL

# --- The null hypothesis ------------------------------------------------------

#' A null hypothesis: a family together with its parameter region
#'
#' \eqn{H_0 = \{P_\theta : \theta \in \bigcup_i \Theta_{0i}\}}{H_0 = {P_theta : theta in union_i Theta_0i}}.
#' The null is a set of *distributions*, so it needs both the model and the
#' geometry; keeping them together means nothing downstream has to carry them as
#' separate arguments that could disagree.
#'
#' @param family A [parametric_family].
#' @param subnulls A non-empty list of [parameter_space] objects, the convex
#'   pieces of the null.
#' @return A `null_model`.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' # The plurality null
#' null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' @export
null_model <- new_class(
  "null_model",
  properties = list(family = parametric_family, subnulls = class_list),
  validator = function(self) {
    if (length(self@subnulls) == 0L) {
      return("`subnulls` must be a non-empty list")
    }
    ok <- vapply(
      self@subnulls,
      \(s) S7_inherits(s, parameter_space),
      logical(1)
    )
    if (!all(ok)) {
      return("every element of `subnulls` must be a `parameter_space`")
    }
    # A region of the wrong dimension is not a subset of the family's parameter
    # space, and comparing it against a parameter would silently recycle rather
    # than complain. Checked here so the error names the call that built it.
    d <- space_dim(self@family@parameter_space)
    dims <- vapply(self@subnulls, space_dim, integer(1))
    if (any(dims != d)) {
      return(paste0(
        "every element of `subnulls` must have dimension ",
        d,
        ", matching the family's parameter space; got ",
        paste(unique(dims[dims != d]), collapse = ", ")
      ))
    }
    NULL
  }
)


#' Number of convex pieces in a null
#' @param null A [null_model].
#' @return Integer.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' n_subnulls(plurality)
#' @export
n_subnulls <- function(null) length(null@subnulls)


#' Does any subnull contain this parameter value?
#' @param null A [null_model].
#' @param theta Parameter vector.
#' @param tol Tolerance.
#' @return `TRUE` or `FALSE`.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' null <- null_model(fam, list(halfspace_region(c(1, -1, 0))))
#' in_null(null, c(0.2, 0.5, 0.3))
#' in_null(null, c(0.6, 0.2, 0.2))
#' @export
in_null <- function(null, theta, tol = 1e-8) {
  any(vapply(null@subnulls, \(s) contains(s, theta, tol), logical(1)))
}
