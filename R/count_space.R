#' @include space.R
NULL

# --- Count vectors ------------------------------------------------------------

#' Enumerate every count vector with `k` categories summing to `n`
#'
#' Stars and bars: each count vector is a choice of `k - 1` bar positions among
#' `n + k - 1` slots, so enumeration costs `O(M * k)`.
#'
#' @param n Total count.
#' @param k Number of categories.
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


#' The space of `k`-category count vectors summing to `n`
#'
#' This is the sample space for multinomial or multivariate hypergeometric: the
#' non-negative integer lattice points of the scaled simplex.
#' @param n Integer total count per outcome.
#' @param k Integer number of categories.
#' @return A `count_space`.
#' @examples
#' count_space(n = 4L, k = 3L)
#' enumerate_space(count_space(n = 2L, k = 3L))
#' @export
count_space <- new_class(
  "count_space",
  parent = space,
  properties = list(n = class_numeric, k = class_numeric, outcomes = class_any),
  constructor = function(n, k) {
    n <- as.integer(n)
    k <- as.integer(k)
    stopifnot(
      "`n` must be a single non-negative integer" = length(n) == 1L &&
        !is.na(n) &&
        n >= 0L,
      "`k` must be a single integer >= 1" = length(k) == 1L &&
        !is.na(k) &&
        k >= 1L
    )
    new_object(S7_object(), n = n, k = k, outcomes = enumerate_counts(n, k))
  }
)


method(space_dim, count_space) <- function(space) as.integer(space@k)


method(is_finite_space, count_space) <- function(space) TRUE


#' @description A count vector belongs when its entries are non-negative whole
#'   numbers summing to `n`. This is the predicate form of the checks
#'   `validate_outcome()` makes; `tol` is accepted for the generic's signature
#'   and unused, since the conditions are exact.
#' @rdname contains
#' @usage NULL
method(contains, count_space) <- function(space, theta, tol = 1e-8) {
  length(theta) == space_dim(space) &&
    all(is.finite(theta)) &&
    all(theta >= 0) &&
    all(theta == trunc(theta)) &&
    sum(theta) == space@n
}


method(enumerate_space, count_space) <- function(space) space@outcomes


#' The space names the total in its message rather than calling it `n_trials`:
#' it cannot know what the family reading it calls that number.
#' @keywords internal
#' @noRd
method(validate_outcome, count_space) <- function(space, x) {
  x <- check_outcome_shape(x, space_dim(space))
  if (any(x < 0) || any(x != trunc(x))) {
    stop("outcomes must be non-negative whole numbers.", call. = FALSE)
  }
  totals <- rowSums(x)
  if (any(totals != space@n)) {
    stop(
      "outcomes must be counts summing to ",
      space@n,
      "; got ",
      paste(unique(totals[totals != space@n]), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  x
}

# --- Real vectors -------------------------------------------------------------
