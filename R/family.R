#' @include space.R region.R
NULL

#' Parametric families
#'
#' A `parametric_family` defines the model \eqn{p_\theta(x)}{p_theta(x)}. It has
#' no knowledge of null hypotheses, alternatives, or any optimisation procedure.
#' Families provide a log-likelihood compiler, score functions and a sampler.
#'
#' A family is the pair of a [convex_region] \eqn{\Theta}{Theta} and the map
#' \eqn{\theta \mapsto p_\theta}{theta -> p_theta} into laws on a
#' [space]; the two spaces are what the family carries, and everything
#' else it offers is a way of navigating that map.
#'
#' Families are callable, which is that map written down: `fam(theta)` is the
#' [distribution] \eqn{p_\theta}{p_theta}. A kernel extends canonically from
#' points to distributions, so `fam(W)` for a [distribution] `W`` over the
#' parameter space is the same map and gives the [induced_distribution()]
#' \eqn{P_W}{P_W}.
#'
#' Not marked abstract, because S7 forbids that alongside a `class_function`
#' parent -- abstract classes must have abstract parents. It is one in every
#' other sense: it supplies no [compile_loglik()] method, so constructing it
#' directly gives a family with no kernel, which errors on first use exactly as
#' any other incomplete family does.
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
#'   components = cbind(c(0.6, 0.2, 0.2), c(0.2, 0.6, 0.2)),
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


#' The map `theta -> p_theta`, shared by every family
#'
#' Defined once at namespace level and never inside a constructor. A closure
#' built per family would capture that constructor's frame, and the family's
#' own properties with it: on a family with a 20k-row property that is roughly
#' twice the serialised size, for a copy nothing reads.
#'
#' `sys.function()` recovers the family being called, with its S7 attributes
#' intact, so the closure captures nothing at all.
#' @keywords internal
#' @noRd
at_theta <- function(at) {
  induced_distribution(sys.function(), at)
}


#' @rdname parametric_family
#' @usage NULL
method(print, parametric_family) <- function(x, ...) {
  cat("<", attr(S7_class(x), "name"), ">\n", sep = "")
  cat("  parameters ", space_label(x@parameter_space), "\n", sep = "")
  cat("  outcomes   ", space_label(x@sample_space), "\n", sep = "")
  invisible(x)
}


#' @description `format()` gives the two spaces on one line, without the class
#'   banner `print()` adds.
#'
#' Both are needed rather than inherited: the parent is `class_function`, so
#' the defaults reach `deparse()` and print the shared closure plus an
#' attribute dump.
#' @rdname parametric_family
#' @usage NULL
method(format, parametric_family) <- function(x, ...) {
  sprintf(
    "%s: %s -> %s",
    attr(S7_class(x), "name"),
    attr(S7_class(x@parameter_space), "name"),
    attr(S7_class(x@sample_space), "name")
  )
}


#' Compile the log-likelihood function for a fixed set of outcomes
#'
#' The single density method a family must implement; [kernel_loglik_batch()] and
#' [kernel_loglik()] are thin wrappers over it. Returns a function of `theta_mat`.
#'
#' Compiling lets a family precompute whatever depends on `x` alone, which we
#' require because the our optimiser fixes the outcomes (either enumerating the
#' entire sample space or fixing monte carlo draws at the start) and evaluates
#' likelihoods for many different parameter values. Closing over `x` rather than
#' taking a cache argument means precomputed constants cannot be paired with the
#' wrong outcomes.
#' @param family A [parametric_family].
#' @param x `(M, K)` matrix of outcomes, where `M` is the number of outcomes and
#'   `K` the dimension of the sample space.
#' @return A function of `theta_mat`, a `(d, C)` matrix of parameter columns,
#'   returning the `(M, C)` matrix of log densities at `x`.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' x <- rbind(c(2L, 1L, 1L), c(4L, 0L, 0L))
#' ll <- compile_loglik(fam, x)
#' ll(cbind(c(0.5, 0.3, 0.2), c(0.25, 0.25, 0.5)))
#' @export
compile_loglik <- new_generic("compile_loglik", "family", function(family, x) {
  S7::S7_dispatch()
})

#' Batched log density over parameter columns
#'
#' Recompiles on every call, so in a loop over many `theta` at fixed `x`,
#' use [compile_loglik()] once and call its result instead.
#'
#' @param family A [parametric_family].
#' @param theta_mat `(d, C)` matrix of parameter columns, where `d` is the
#'   dimension of the parameter and `C` the number of columns. The same `x` is
#'   used for every column.
#' @param x `(M, K)` matrix of outcomes.
#' @return `(M, C)` matrix of log densities.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' x <- rbind(c(2L, 1L, 1L), c(4L, 0L, 0L))
#' kernel_loglik_batch(fam, cbind(c(0.5, 0.3, 0.2), c(0.25, 0.25, 0.5)), x)
#' @export
kernel_loglik_batch <- function(family, theta_mat, x) {
  compile_loglik(family, x)(theta_mat)
}


#' Log density `log p_theta(x)`
#' @param family A [parametric_family].
#' @param theta Parameter vector of length `space_dim(family@parameter_space)`.
#' @param x `(M, K)` matrix of outcomes, or a length-`K` vector for one outcome.
#' @return Length-`M` numeric vector.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' kernel_loglik(fam, c(0.5, 0.3, 0.2), c(2L, 1L, 1L))
#' @export
kernel_loglik <- function(family, theta, x) {
  as.vector(kernel_loglik_batch(family, matrix(theta, ncol = 1L), x))
}


#' Score `d log P_theta(x) / d theta`
#'
#' Per-outcome contributions in the family's own parameter coordinates, with no
#' constraint projection applied. Applying the Jacobian of a parametrisation
#' belongs to whatever owns that parametrisation.
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
#' Take one draw per column of `theta_mat`, so the number of draws is the number
#' of parameters. To take repeated draws from a single parameter value you would
#' repeat the parameter value across columns, which is what
#' [induced_distribution()] does for a point mass.
#'
#' Users should sample via `draw(fam(theta), n)`, which does the same thing and
#' routes here.
#' @param family A [parametric_family].
#' @param theta_mat `(d, M)` matrix of parameter columns; a length-`d` vector is
#'   taken as a single column.
#' @return `(M, k)` numeric matrix, one observation per row.
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#'
#' # Five draws from one parameter: repeat it across five columns.
#' kernel_draw(fam, matrix(c(0.5, 0.3, 0.2), nrow = 3L, ncol = 5L))
#'
#' # One draw from each of three different parameters.
#' kernel_draw(fam, cbind(c(0.5, 0.3, 0.2), c(0.2, 0.2, 0.6), c(0.9, 0.05, 0.05)))
#' @seealso [compile_loglik()], the density half of the same pair.
#' @export
kernel_draw <- new_generic(
  "kernel_draw",
  "family",
  function(family, theta_mat) {
    S7::S7_dispatch()
  }
)
