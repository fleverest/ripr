#' @include space.R
NULL

# --- Count vectors ------------------------------------------------------------

#' Every count vector with `k` categories summing to `n`, by stars and bars.
#'
#' Rows are in ascending lexicographic order. The order matches
#' `bernstein_lattice()`, so `certify()` can use values on the sample space
#' directly as Bernstein coefficients.
#' @return `(M, k)` integer matrix, `M = choose(n + k - 1, k - 1)`.
#' @keywords internal
#' @noRd
enumerate_counts <- function(n, k) {
  if (k == 1L) {
    return(matrix(as.integer(n), nrow = 1L))
  }
  bars <- utils::combn(n + k - 1L, k - 1L)
  t(diff(rbind(0L, bars, n + k)) - 1L)
}


#' The space of `k`-category count vectors summing to `n_trials`
#'
#' The sample space of a multinomial or multivariate hypergeometric: the
#' non-negative integer points of the scaled simplex.
#' @param n_trials Integer total count per outcome (the number of trials).
#' @param k Integer number of categories.
#' @return A `count_space`.
#' @examples
#' count_space(n_trials = 4L, k = 3L)
#' enumerate_space(count_space(n_trials = 2L, k = 3L))
#' @export
count_space <- new_class(
  "count_space",
  parent = space,
  properties = list(
    n_trials = class_numeric,
    k = class_numeric,
    outcomes = class_any
  ),
  constructor = function(n_trials, k) {
    rlang::check_number_whole(n_trials, min = 0)
    rlang::check_number_whole(k, min = 1)
    n_trials <- as.integer(n_trials)
    k <- as.integer(k)
    new_object(
      S7_object(),
      n_trials = n_trials,
      k = k,
      outcomes = enumerate_counts(n_trials, k)
    )
  }
)


method(space_dim, count_space) <- function(space) as.integer(space@k)


method(is_finite_space, count_space) <- function(space) TRUE


#' @description A count vector belongs when its entries are non-negative whole
#'   numbers summing to `n_trials`; `tol` is unused.
#' @rdname contains
#' @usage NULL
method(contains, count_space) <- function(space, theta, tol = 1e-8) {
  length(theta) == space_dim(space) &&
    all(is.finite(theta)) &&
    all(theta >= 0) &&
    all(theta == trunc(theta)) &&
    sum(theta) == space@n_trials
}


method(enumerate_space, count_space) <- function(space) space@outcomes


method(validate_outcome, count_space) <- function(space, x) {
  x <- check_outcome_shape(x, space_dim(space))
  if (any(x < 0) || any(x != trunc(x))) {
    stop("outcomes must be non-negative whole numbers.", call. = FALSE)
  }
  totals <- rowSums(x)
  if (any(totals != space@n_trials)) {
    stop(
      "outcomes must be counts summing to ",
      space@n_trials,
      "; got ",
      paste(unique(totals[totals != space@n_trials]), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  x
}

# --- Real vectors -------------------------------------------------------------
