#' Spaces
#'
#' A `space` is a measurable space, carrying a dimension ([space_dim()]) along
#' with membership tests ([contains()]). It is either a finite [count_space()]
#' that can be listed ([enumerate_space()]), or a [region] of
#' \eqn{\mathbb{R}^d}{R^d}. A space may serve as the sample space of a
#' [distribution] or the parameter space of a [family]; [simplex_region()] is
#' both, for the multinomial family and for [dirichlet].
#' @examples
#' space_dim(count_space(n_trials = 4L, k = 3L))
#' contains(count_space(n_trials = 4L, k = 3L), c(2L, 1L, 1L))
#' @export
space <- new_class("space", abstract = TRUE)


#' The dimension of one element of a space
#' @param space A [space].
#' @return An integer.
#' @examples
#' space_dim(count_space(n_trials = 4L, k = 3L))
#' space_dim(simplex_region(vertices = diag(3)))
#' @export
space_dim <- new_generic("space_dim", "space", function(space) {
  S7::S7_dispatch()
})


#' Does a point belong to a space?
#'
#' Checks whether `theta` belongs to `space`.
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
#' contains(count_space(n_trials = 4L, k = 3L), c(2L, 1L, 1L))
#' contains(count_space(n_trials = 4L, k = 3L), c(2L, 1L, 0L))
#' @export
contains <- new_generic(
  "contains",
  "space",
  function(space, theta, tol = 1e-8) S7::S7_dispatch()
)


#' Coerce and check elements of a sample space
#'
#' Takes a length-`d` vector or `(n, d)` matrix and returns the `(n, d)` form;
#' anything outside the sample space is an error.
#' @keywords internal
#' @noRd
validate_outcome <- new_generic(
  "validate_outcome",
  "space",
  function(space, x) {
    S7::S7_dispatch()
  }
)


#' Shape checks common to every sample space
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
#' enumerate_space(count_space(n_trials = 3L, k = 2L))
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
#' Whether a sample space may be enumerated.
#' @seealso [enumerate_space()]
#' @param space A [space].
#' @return `TRUE` or `FALSE`.
#' @examples
#' is_finite_space(count_space(n_trials = 4L, k = 3L))
#' is_finite_space(real_region(1L))
#' @export
is_finite_space <- new_generic("is_finite_space", "space", function(space) {
  S7::S7_dispatch()
})


method(is_finite_space, space) <- function(space) FALSE


#' A short description of a space, for printing
#' @keywords internal
#' @noRd
space_label <- function(space) {
  sprintf("%s, dimension %d", class_name(space), space_dim(space))
}


#' Name a class as it should appear in a message
#' @keywords internal
#' @noRd
class_name <- function(x) attr(S7_class(x), "name")
