#' Sampling family: the observation model
#'
#' A `sampling_family` owns the model \eqn{p_\theta(x)}{p_theta(x)} and nothing
#' else; no knowledge of null hypotheses, alternatives, or any optimisation.
#' Families supply a compiled log-likelihood, the score, a sampler, and support
#' enumeration (for finite families).
#'
#' @examples
#' # `sampling_family` is abstract; families subclass it, e.g.
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' S7::S7_inherits(fam, sampling_family)
#' @export
sampling_family <- new_class("sampling_family", abstract = TRUE)


#' Dimension of the parameter vector
#' @param family A [sampling_family].
#' @return Integer parameter dimension.
#' @examples
#' param_dim(multinomial_family(n_trials = 4L, k = 3L))
#' @export
param_dim <- new_generic("param_dim", "family", function(family) {
  S7::S7_dispatch()
})


#' Compile the log-likelihood function for a fixed set of outcomes
#'
#' The single density method a family must implement; [log_density_batch()] and
#' [log_density()] are thin wrappers over it. Returns a function of `theta_mat`
#' alone.
#'
#' Compiling lets a family precompute whatever depends on `x` alone, which pays
#' off because the outcomes stay fixed for a whole fit while `theta` moves
#' thousands of times. Closing over `x` rather than taking a cache argument also
#' means precomputed constants cannot be paired with the wrong outcomes.
#'
#' @param family A [sampling_family].
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
#' @param family A [sampling_family].
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
#'
#' @param family A [sampling_family].
#' @param theta Parameter vector of length [param_dim()].
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
#'
#' @param family A [sampling_family].
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
#'
#' Named `draw` rather than `simulate` to avoid masking [stats::simulate()].
#' @param family A [sampling_family].
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

#' Enumerated support of a finite family
#'
#' Every outcome with positive probability under some parameter. Families with
#' infinite sample spaces error: exact integration over the sample space is
#' undefined for them, not merely slow.
#' @param family A [sampling_family].
#' @return `(M, K)` matrix of outcomes.
#' @examples
#' support(multinomial_family(n_trials = 3L, k = 2L))
#' @export
support <- new_generic("support", "family", function(family) S7::S7_dispatch())


method(support, sampling_family) <- function(family) {
  stop(
    "`",
    S7_class(family)@name,
    "` has no enumerable support. Use a Monte ",
    "Carlo or quadrature engine instead of an exact one.",
    call. = FALSE
  )
}

#' Dimension of one element of the sample space
#'
#' Distinct from [param_dim()] in principle, though equal for every family here.
#' @param family A [sampling_family].
#' @return Integer.
#' @examples
#' outcome_dim(multinomial_family(n_trials = 4L, k = 3L))
#' @export
outcome_dim <- new_generic("outcome_dim", "family", function(family) {
  S7::S7_dispatch()
})


method(outcome_dim, sampling_family) <- function(family) param_dim(family)


#' Coerce and check elements of the sample space
#'
#' Accepts one element as a length-`d` vector or `n` of them as an `(n, d)`
#' matrix, and returns the `(n, d)` form. Anything else is an error: a random
#' variable is only defined on its own sample space, and silently reshaping the
#' wrong thing would give a number rather than a complaint.
#' @param family A [sampling_family].
#' @param x A length-`d` vector or `(n, d)` matrix.
#' @return `(n, d)` numeric matrix.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' as_outcomes(fam, c(2L, 1L, 1L))
#' @export
as_outcomes <- new_generic("as_outcomes", "family", function(family, x) {
  S7::S7_dispatch()
})

#' Shape checks common to every family
#'
#' A plain function rather than a method, so the per-family methods can call it
#' without dispatching to a parent. `S7::super()` would do the same thing, but
#' its behaviour varies between S7 versions, and when dispatch fails S7 reports
#' it by describing the offending object -- a path that can itself fail, hiding
#' the real error behind a formatting one.
#' @keywords internal
#' @noRd
check_outcome_shape <- function(x, d) {
  if (!is.numeric(x)) {
    stop("outcomes must be numeric.", call. = FALSE)
  }
  if (!is.matrix(x)) {
    if (length(x) != d) {
      stop(
        "one outcome must be a length-",
        d,
        " vector; got length ",
        length(x),
        ".",
        call. = FALSE
      )
    }
    x <- matrix(x, nrow = 1L)
  }
  if (ncol(x) != d) {
    stop(
      "outcomes must have ",
      d,
      " columns; got ",
      ncol(x),
      ".",
      call. = FALSE
    )
  }
  if (anyNA(x)) {
    stop("outcomes must not be missing.", call. = FALSE)
  }
  x
}


method(as_outcomes, sampling_family) <- function(family, x) {
  check_outcome_shape(x, outcome_dim(family))
}


#' Does this family have a finite, enumerable sample space?
#' @param family A [sampling_family].
#' @return `TRUE` or `FALSE`.
#' @examples
#' is_finite_support(multinomial_family(n_trials = 4L, k = 3L))
#' is_finite_support(gaussian_family(dim = 1L))
#' @export
is_finite_support <- new_generic(
  "is_finite_support",
  "family",
  function(family) {
    S7::S7_dispatch()
  }
)

method(is_finite_support, sampling_family) <- function(family) FALSE

#' A reference point for the parameter space.
#'
#' Used as a fallback for initialising the atoms for a RIPr optimisation run,
#' where the alternative does not take the form of a mixture.
#' @keywords internal
#' @noRd
reference_parameter <- new_generic("reference_parameter", "family", \(family) {
  S7::S7_dispatch()
})
