#' RIPr e-variable
#'
#' Bundles a numerator distribution `Q` and its reverse information projection
#' `P*` into a single object that evaluates the e-variable on data. The rescaled
#' e-variable `(Q / P*) / (1 + gap)` has expectation at most 1 under every null
#' distribution; `gap` is the certified duality gap the projection achieved
#' (`0` leaves the raw likelihood ratio uncorrected).
#'
#' Usually obtained by fitting with [run_ripr()] (as the `$e_variable` element of
#' its result), but can also be assembled directly from the pieces.
#'
#' @param numerator The numerator distribution `Q` (a [distribution]).
#' @param projection The reverse information projection `P*` (a [distribution],
#'   typically a [mixture_dist]).
#' @param family The sampling [family] both distributions live over.
#' @param gap The certified duality gap used to rescale. Default `0`.
#' @return An `e_variable`.
#' @export
e_variable <- new_class(
  "e_variable",
  properties = list(
    numerator = distribution,
    projection = distribution,
    family = family,
    gap = class_numeric
  ),
  constructor = function(numerator, projection, family, gap = 0) {
    new_object(
      S7_object(),
      numerator = numerator,
      projection = projection,
      family = family,
      gap = as.numeric(gap)
    )
  }
)

#' Evaluate an e-variable on outcomes
#'
#' Computes the e-value `e(x) = (Q(x) / P*(x)) / (1 + gap)` at each outcome `x`.
#'
#' @param ev An `e_variable`.
#' @param x Outcomes: an `(N, d)` matrix, or a length-d vector for a single
#'   outcome.
#' @param log Return the log e-value instead? Default `FALSE`.
#' @param corrected Apply the `1 / (1 + gap)` rescaling that makes the e-variable
#'   valid under the null? Default `TRUE`. With `FALSE` the raw likelihood ratio
#'   `Q(x) / P*(x)` (the e-variable before the duality-gap correction) is
#'   returned.
#' @return Numeric vector of (log) e-values, one per outcome.
#' @export
e_value <- new_generic(
  "e_value",
  "ev",
  function(ev, x, log = FALSE, corrected = TRUE) {
    S7::S7_dispatch()
  }
)

method(e_value, e_variable) <- function(ev, x, log = FALSE, corrected = TRUE) {
  log_e <- dist_log_density(ev@numerator, ev@family, x) -
    dist_log_density(ev@projection, ev@family, x)
  if (corrected) {
    log_e <- log_e - log1p(ev@gap)
  }
  if (log) log_e else exp(log_e)
}

method(print, e_variable) <- function(x, ...) {
  proj <- x@projection
  n_atoms <- if (S7_inherits(proj, mixture_dist)) {
    ncol(proj@components)
  } else {
    NA_integer_
  }
  cat("<e_variable>  e(x) = (Q / P*) / (1 + gap)\n")
  cat(sprintf("  numerator  Q  : %s\n", S7_class(x@numerator)@name))
  cat(sprintf(
    "  projection P* : %s%s\n",
    S7_class(proj)@name,
    if (is.na(n_atoms)) "" else sprintf(" (%d atom%s)", n_atoms, if (n_atoms == 1L) "" else "s")
  ))
  cat(sprintf("  gap           : %.4g   (correction 1 + gap = %.6g)\n", x@gap, 1 + x@gap))
  invisible(x)
}
