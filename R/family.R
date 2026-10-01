#' @include space.R region.R
NULL

#' Parametric families
#'
#' A `parametric_family` is the model \eqn{p_\theta(x)}{p_theta(x)}: a
#' parameter space \eqn{\Theta}{Theta} (a [convex_region]), a sample [space],
#' and the map \eqn{\theta \mapsto p_\theta}{theta -> p_theta} between them.
#' It provides a log-likelihood compiler, score and sampler.
#'
#' Families are callable: `fam(theta)` is the [distribution]
#' \eqn{p_\theta}{p_theta}, and `fam(W)` for a [distribution] `W` over the
#' parameter space is the [mixture()] \eqn{P_W}{P_W}.
#'
#' This is effectively an abstract class. It has no [compile_loglik()] method,
#' so use a concrete family instead.
#'
#' @param sample_space The [space] that outcomes belong to.
#' @param parameter_space The [convex_region] that parameter lives in. For
#'   instance, the standard simplex for Multinomial proportions.
#' @return A callable `parametric_family`.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' S7::S7_inherits(fam, parametric_family)
#'
#' # The map theta -> p_theta, and its extension to mixing measures.
#' fam(c(0.5, 0.3, 0.2))
#' fam(finite_dist(
#'   atoms = rbind(c(0.6, 0.2, 0.2), c(0.2, 0.6, 0.2)),
#'   weights = c(0.5, 0.5)
#' ))
#'
#' # Properties and other generics are unaffected by being callable.
#' enumerate_space(fam@sample_space)
#' @export
parametric_family <- new_class(
  "parametric_family",
  parent = class_function,
  properties = list(
    sample_space = space,
    parameter_space = convex_region
  ),
  constructor = function(sample_space, parameter_space) {
    new_object(
      at_theta,
      sample_space = sample_space,
      parameter_space = parameter_space
    )
  }
)


#' The map `theta -> p_theta`, shared by every family. Defined at namespace
#' level, not in a constructor, so a family does not serialise a copy of the
#' constructor's frame; `sys.function()` recovers the family being called.
#' @keywords internal
#' @noRd
at_theta <- function(at) {
  mixture(sys.function(), at)
}


#' @rdname parametric_family
#' @usage NULL
method(print, parametric_family) <- function(x, ...) {
  cat("<", class_name(x), ">\n", sep = "")
  cat("  parameters ", space_label(x@parameter_space), "\n", sep = "")
  cat("  outcomes   ", space_label(x@sample_space), "\n", sep = "")
  invisible(x)
}


#' @description `format()` gives the two spaces on one line, without the class
#'   banner `print()` adds.
#' @rdname parametric_family
#' @usage NULL
method(format, parametric_family) <- function(x, ...) {
  sprintf(
    "%s: %s -> %s",
    class_name(x),
    class_name(x@parameter_space),
    class_name(x@sample_space)
  )
}


#' Compile the log-likelihood function for a fixed set of outcomes
#'
#' The density method that a [parametric_family] implements. Compiling lets the
#' family precompute whatever depends on `x`. For the density at one parameter,
#' use `log_density(family(theta), x)`.
#' @param family A [parametric_family].
#' @param x `(M, K)` matrix of outcomes, where `M` is the number of outcomes and
#'   `K` the dimension of the sample space.
#' @return A function of `theta_mat`, a `(C, d)` matrix with one parameter
#'   vector per row, returning the `(M, C)` matrix of log densities at `x`:
#'   one row per outcome, one column per parameter.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' x <- rbind(c(2L, 1L, 1L), c(4L, 0L, 0L))
#' ll <- compile_loglik(fam, x)
#' ll(rbind(c(0.5, 0.3, 0.2), c(0.25, 0.25, 0.5)))
#' @export
compile_loglik <- new_generic("compile_loglik", "family", function(family, x) {
  S7::S7_dispatch()
})

#' `(M, C)` log densities at `x` for each row of `theta_mat`. Recompiles on
#' every call; in a loop, call [compile_loglik()] once instead.
#' @keywords internal
#' @noRd
kernel_loglik_batch <- function(family, theta_mat, x) {
  compile_loglik(family, x)(theta_mat)
}


#' Score `d log P_theta(x) / d theta`
#'
#' Per-outcome contributions in the family's own parameter coordinates, with no
#' constraint projection applied.
#' @param family A [parametric_family].
#' @param theta Parameter vector.
#' @param x `(M, K)` matrix of outcomes.
#' @return `(M, d)` matrix.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' score(fam, c(0.5, 0.3, 0.2), c(2L, 1L, 1L))
#' @export
score <- new_generic("score", "family", function(family, theta, x) {
  S7::S7_dispatch()
})


#' Draw one observation from `P_theta` per parameter
#'
#' One draw per row of `theta_mat`. Repeat rows for repeated draws. Usually
#' reached through `draw(fam(theta), n)`.
#' @param family A [parametric_family].
#' @param theta_mat `(M, d)` matrix with one parameter vector per row; a
#'   length-`d` vector is taken as a single row.
#' @return `(M, K)` numeric matrix, one observation per row, drawn from the
#'   parameter in the matching row of `theta_mat`.
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#'
#' # Five draws from one parameter: repeat it across five rows.
#' kernel_draw(fam, matrix(c(0.5, 0.3, 0.2), nrow = 5L, ncol = 3L, byrow = TRUE))
#'
#' # One draw from each of three different parameters.
#' kernel_draw(fam, rbind(c(0.5, 0.3, 0.2), c(0.2, 0.2, 0.6), c(0.9, 0.05, 0.05)))
#' @seealso [compile_loglik()], the density half of the same pair.
#' @export
kernel_draw <- new_generic(
  "kernel_draw",
  "family",
  function(family, theta_mat) {
    S7::S7_dispatch()
  }
)
