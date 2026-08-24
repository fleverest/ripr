#' Sample spaces
#'
#' A `sample_space` defines the set \eqn{\mathcal{X}}{X} that observation are
#' drawn from, and nothing about how it is drawn. It answers three questions:
#' how many coordinates one element has ([space_dim()]), whether a given object
#' is an element ([validate_outcome()]) and, when the space is finite,
#' ([is_finite_space()]), what all of the possible outcomes are
#' ([enumerate_space()]).
#'
#' This is deliberately separate from [parametric_family]. The role that
#' sample spaces take in this package is effectively just defining the outcome
#' space for a distribution, and the object is in charge of validation and
#' carrying data needed upstream.
#' @examples
#' # `sample_space` is abstract; spaces subclass it, e.g.
#' S7::S7_inherits(count_space(n = 4L, k = 3L), sample_space)
#' @export
sample_space <- new_class("sample_space", abstract = TRUE)


#' The dimension of one element of the sample space
#' @param space A [sample_space].
#' @return An integer.
#' @examples
#' space_dim(count_space(n = 4L, k = 3L))
#' space_dim(real_space(2L))
#' @export
space_dim <- new_generic("space_dim", "space", function(space) {
  S7::S7_dispatch()
})


#' Coerce and check elements of a sample space
#'
#' Accepts one element as a length-`d` vector or `n` of them as an `(n, d)`
#' matrix, and returns the `(n, d)` form. Anything else is an error: a random
#' variable is only defined on its own sample space.
#' @param space A [sample_space].
#' @param x A length-`d` vector or `(n, d)` matrix.
#' @return `(n, d)` numeric matrix.
#' @examples
#' validate_outcome(count_space(n = 4L, k = 3L), c(2L, 1L, 1L))
#' @export
validate_outcome <- new_generic(
  "validate_outcome",
  "space",
  function(space, x) {
    S7::S7_dispatch()
  }
)


#' Shape checks common to every sample space
#'
#' A plain function rather than a method, so the per-space methods can call it
#' without dispatching to a parent.
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


method(validate_outcome, sample_space) <- function(space, x) {
  check_outcome_shape(x, space_dim(space))
}


#' Every element of a finite sample space
#'
#' Infinite spaces error.
#' @param space A [sample_space].
#' @return `(M, d)` matrix, one outcome per row.
#' @examples
#' enumerate_space(count_space(n = 3L, k = 2L))
#' @export
enumerate_space <- new_generic("enumerate_space", "space", function(space) {
  S7::S7_dispatch()
})


method(enumerate_space, sample_space) <- function(space) {
  stop(
    "`",
    S7_class(space)@name,
    "` cannot be enumerated. Use a Monte Carlo or ",
    "quadrature engine instead of an exact one.",
    call. = FALSE
  )
}


#' Is this sample space finite?
#'
#' Use to determine whether a sample space may be enumerated.
#' @seealso [enumerate_space()]
#' @param space A [sample_space].
#' @return `TRUE` or `FALSE`.
#' @examples
#' is_finite_space(count_space(n = 4L, k = 3L))
#' is_finite_space(real_space(1L))
#' @export
is_finite_space <- new_generic("is_finite_space", "space", function(space) {
  S7::S7_dispatch()
})


method(is_finite_space, sample_space) <- function(space) FALSE


#' A short description of a space, for printing
#'
#' Works for a [region] too: all it needs is a class name and
#' [space_dim()].
#' @keywords internal
#' @noRd
space_label <- function(space) {
  name <- attr(S7_class(space), "name")
  sprintf("%s, dimension %d", name, space_dim(space))
}


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
  parent = sample_space,
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

#' The space of `d`-dimensional real numbers: \eqn{\mathbb{R}^d}{R^d}
#'
#' @param d Integer dimension.
#' @return A `real_space`.
#' @examples
#' real_space(2L)
#' @export
real_space <- new_class(
  "real_space",
  parent = sample_space,
  properties = list(n_dim = class_numeric),
  constructor = function(d) {
    d <- as.integer(d)
    stopifnot(
      "`d` must be a single positive integer" = length(d) == 1L &&
        !is.na(d) &&
        d >= 1L
    )
    new_object(S7_object(), n_dim = d)
  }
)


method(space_dim, real_space) <- function(space) as.integer(space@n_dim)


#' All of \eqn{\mathbb{R}^d}{R^d} is admissible, so the inherited checks (for
#' being numeric, correct shape, no missing values) is basically all we need.
#' All we need to add is finiteness.
#' @keywords internal
#' @noRd
method(validate_outcome, real_space) <- function(space, x) {
  x <- check_outcome_shape(x, space_dim(space))
  if (any(!is.finite(x))) {
    stop("outcomes must be finite.", call. = FALSE)
  }
  x
}
