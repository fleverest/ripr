#' @include sample_space.R parameter_space.R
NULL

#' Parametric families
#'
#' A `parametric_family` defines the model \eqn{p_\theta(x)}{p_theta(x)}. It has
#' no knowledge of null hypotheses, alternatives, or any optimisation procedure.
#' Families provide a log-likelihood compiler, score functions and a sampler.
#'
#' A family is the pair of a [parameter_space] \eqn{\Theta}{Theta} and the map
#' \eqn{\theta \mapsto p_\theta}{theta -> p_theta} into laws on a
#' [sample_space]; the two spaces are what the family carries, and everything
#' else it offers is a way of navigating that map.
#' @param sample_space The [sample_space] that outcomes belong to.
#' @param parameter_space The [parameter_space] that parameters belong to.
#' @examples
#' # `parametric_family` is abstract; families subclass it, e.g.
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' S7::S7_inherits(fam, parametric_family)
#' enumerate_space(fam@sample_space)
#' @export
parametric_family <- new_class(
  "parametric_family",
  abstract = TRUE,
  properties = list(
    sample_space = sample_space,
    parameter_space = parameter_space
  )
)


#' Compile the log-likelihood function for a fixed set of outcomes
#'
#' The single density method a family must implement; [log_density_batch()] and
#' [log_density()] are thin wrappers over it. Returns a function of `theta_mat`.
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
#' log_density_batch(fam, cbind(c(0.5, 0.3, 0.2), c(0.25, 0.25, 0.5)), x)
#' @export
log_density_batch <- function(family, theta_mat, x) {
  compile_loglik(family, x)(theta_mat)
}


#' Log density `log p_theta(x)`
#' @param family A [parametric_family].
#' @param theta Parameter vector of length `space_dim(family@parameter_space)`.
#' @param x `(M, K)` matrix of outcomes, or a length-`K` vector for one outcome.
#' @return Length-`M` numeric vector.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' log_density(fam, c(0.5, 0.3, 0.2), c(2L, 1L, 1L))
#' @export
log_density <- function(family, theta, x) {
  as.vector(log_density_batch(family, matrix(theta, ncol = 1L), x))
}


#' Score `d log p_theta(x) / d theta`
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


#' Draw observations from `p_theta`
#' @param family A [parametric_family].
#' @param theta Parameter vector.
#' @param n_obs Number of draws.
#' @return `(n_obs, k)` numeric matrix.
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' draw(fam, c(0.5, 0.3, 0.2), n_obs = 5L)
#' @export
draw <- new_generic("draw", "family", function(family, theta, n_obs) {
  S7::S7_dispatch()
})


#' A reference point for the parameter space.
#'
#' Used as a fallback for initialising the atoms for a RIPr optimisation run,
#' where the alternative does not take the form of a mixture. Defaults to the
#' point of the parameter space closest to the origin, which is the centroid for
#' a simplex and the origin itself for an unconstrained space.
#' @keywords internal
#' @noRd
reference_parameter <- new_generic("reference_parameter", "family", \(family) {
  S7::S7_dispatch()
})


method(reference_parameter, parametric_family) <- function(family) {
  space <- family@parameter_space
  project(space, rep(0, space_dim(space)))
}
