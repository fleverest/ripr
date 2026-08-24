#' @include region.R family.R
NULL

# --- The null hypothesis ------------------------------------------------------

#' A null hypothesis: a family together with its parameter region
#'
#' \eqn{H_0 = \{P_\theta : \theta \in \bigcup_i \Theta_{0i}\}}{H_0 = {P_theta : theta in union_i Theta_0i}}.
#' The null is a set of *distributions*, so it needs both the model and the
#' geometry; keeping them together means nothing downstream has to carry them as
#' separate arguments that could disagree.
#'
#' Each convex part of the region is itself a null hypothesis. This package
#' calls them parts, because [parts()] is what any [region] answers whether
#' or not it happens to be a null, but the two words mean the same thing here.
#'
#' @param family A [parametric_family].
#' @param region The null's geometry: any [region]. A single [convex_region]
#'   is stored as it comes, since one convex set is already a region; a list of
#'   them becomes the [union_region] of its elements.
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
  properties = list(family = parametric_family, region = region),
  constructor = function(family, region) {
    new_object(S7_object(), family = family, region = as_region(region))
  },
  validator = function(self) {
    # `union_region` already validates that its cells are convex regions of one
    # shared dimension; all that is left is whether it is compatible with
    # `family`.
    d <- space_dim(self@family@parameter_space)
    d_region <- space_dim(self@region)
    if (d_region != d) {
      return(paste0(
        "the null's region must have dimension ",
        d,
        ", matching the family's parameter space; got ",
        d_region
      ))
    }
    NULL
  }
)


#' Does the null contain this parameter value?
#'
#' Membership of \eqn{H_0}{H_0}, which is [contains()] on the null's region:
#' true when any one part holds `theta`.
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
  contains(null@region, theta, tol)
}
