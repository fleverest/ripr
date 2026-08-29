#' Spaces
#'
#' A `space` is a set that a measure can live on. A `space` necessarily tracks
#' two things: its dimension ([space_dim()], the length of its elements),
#' and whether a given point is an element or not ([contains()]).
#'
#' There are really only two types of spaces in this package. First, there are
#' [count_space()]s, which are finite and can be listed ([enumerate_space()]).
#' Then there are [region]s that define subsets of the real numbers
#' \eqn{\mathbb{R}^d}{R^d} carrying convex geometry, so it can be charted,
#' projected onto, and composed via set algebra.
#'
#' There are two roles that `space`s fill; they may be the sample space for
#' a [distribution], or the parameter space of a [family]. Take, for example,
#' the [simplex_region()], which is the parameter space for the multinomial
#' family states, and also the sample space for [dirichlet] distributions.
#' @examples
#' space_dim(count_space(n = 4L, k = 3L))
#' contains(count_space(n = 4L, k = 3L), c(2L, 1L, 1L))
#' @export
space <- new_class("space", abstract = TRUE)


#' The dimension of one element of a space
#' @param space A [space].
#' @return An integer.
#' @examples
#' space_dim(count_space(n = 4L, k = 3L))
#' space_dim(simplex_region(vertices = diag(3)))
#' @export
space_dim <- new_generic("space_dim", "space", function(space) {
  S7::S7_dispatch()
})


#' Does a point belong to a space?
#'
#' Checks whether `theta` belongs to `space`.
#'
#' Contrast [validate_outcome()], which asks the same question of a batch and
#' answers by coercing and erroring rather than by returning `TRUE` or `FALSE`;
#' it is built on this.
#' @param space A [space].
#' @param theta A point of the space.
#' @param tol Tolerance.
#' @return `TRUE` or `FALSE`.
#' @examples
#' s <- simplex_region(vertices = diag(3))
#' contains(s, c(1 / 3, 1 / 3, 1 / 3))
#' contains(s, c(2, -1, 0))
#'
#' # A count space is a space too, and answers the same question.
#' contains(count_space(n = 4L, k = 3L), c(2L, 1L, 1L))
#' contains(count_space(n = 4L, k = 3L), c(2L, 1L, 0L))
#' @export
contains <- new_generic(
  "contains",
  "space",
  function(space, theta, tol = 1e-8) S7::S7_dispatch()
)


#' Coerce and check elements of a sample space
#'
#' Accepts one element as a length-`d` vector or `n` of them as an `(n, d)`
#' matrix, and returns the `(n, d)` form. Anything else is an error: a random
#' variable is only defined on its own sample space.
#' @param space A [space].
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


method(validate_outcome, space) <- function(space, x) {
  check_outcome_shape(x, space_dim(space))
}


#' Every element of a finite sample space
#'
#' Infinite spaces error.
#' @param space A [space].
#' @return `(M, d)` matrix, one outcome per row.
#' @examples
#' enumerate_space(count_space(n = 3L, k = 2L))
#' @export
enumerate_space <- new_generic("enumerate_space", "space", function(space) {
  S7::S7_dispatch()
})


method(enumerate_space, space) <- function(space) {
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
#' @param space A [space].
#' @return `TRUE` or `FALSE`.
#' @examples
#' is_finite_space(count_space(n = 4L, k = 3L))
#' is_finite_space(real_region(1L))
#' @export
is_finite_space <- new_generic("is_finite_space", "space", function(space) {
  S7::S7_dispatch()
})


method(is_finite_space, space) <- function(space) FALSE


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
