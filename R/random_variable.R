#' @include family.R mixture.R
NULL

# Random variables on a sample space, and arithmetic over them.
#
# A random variable maps one element of the sample space to one number, or `n`
# elements to `n` numbers. Instances are callable, and check their input before
# mapping it: a variable is only defined on its own sample space, and silently
# reshaping the wrong thing would return a number rather than a complaint.
#
# Nothing here relates to e-values, other than the fact that e-values are
# a random variable. The mixture likelihood ratio `Q / P_star` can be turned
# into an e-variable (for H_0) by dividing by sup_theta E_theta[Q/P_star],
# where the supremum is taken over the null.

#' A random variable on a sample space
#'
#' A `random_variable` `X` takes one element of the sample space (a length-`d`
#' vector for a `d`-dimensional sample space) and returns a number, or `n`
#' elements as an `(n, d)` matrix and returns `n` numbers. The input is checked
#' by [as_outcomes()] first.
#'
#' `Inf` can be a legitimate value here. A likelihood ratio can be genuinely
#' infinite when the denominator is not absolutely continuous with respect to
#' the numerator.
#'
#' @param f The mapping that defines the random variable, accepting an `(n, d)`
#'   matrix and returning `n` numbers.
#' @param family The [sampling_family] whose sample space this is defined on.
#' @param label How to name this variable when printing. Ignored when `op` is
#'   given, since the expression is then built from the operands.
#' @param op The operator that produced this variable, or `NA` for a leaf. Set
#'   by [random_variable_arithmetic]; there is rarely a reason to pass it.
#' @return A callable `random_variable`.
#' @seealso [mixture_likelihood()], [random_variable_arithmetic]
#' @examples
#' X <- random_variable(\(x) dnorm(x, 1), family = gaussian_family(dim = 1))
#' X(as.matrix(0:2))
#' @export
random_variable <- new_class(
  "random_variable",
  parent = class_function,
  properties = list(
    family = sampling_family,
    label = class_character,
    op = class_character,
    operands = class_list
  ),
  constructor = function(
    f,
    family,
    label = "<rv>",
    op = NA_character_,
    operands = list()
  ) {
    # Forced so the closure captures values, not promises: `saveRDS` on a
    # `class_function` parent serialises whatever the environment holds.
    force(f)
    force(family)
    if (!is.function(f)) {
      stop("`f` must be a function.", call. = FALSE)
    }

    new_object(
      function(x) {
        force(x)
        out <- f(as_outcomes(family, x))
        if (!is.numeric(out)) {
          stop("a random variable must return numbers.", call. = FALSE)
        }
        as.vector(out)
      },
      family = family,
      label = label,
      op = op,
      operands = operands
    )
  }
)


# --- Printing -----------------------------------------------------------------

#' Shorten a label that came from deparsing an inline argument
#'
#' `mixture_likelihood(Q)` gives a short label, but an expression written in
#' place gives back the whole expression, which swamps the printed line.
#' @keywords internal
#' @noRd
short_label <- function(label, width = 24L) {
  if (is.na(label) || nchar(label) <= width) {
    return(label)
  }
  paste0(substr(label, 1L, width - 3L), "...")
}


#' Binding strength, for deciding where brackets are needed
#' @keywords internal
#' @noRd
op_precedence <- function(op) if (op %in% c("*", "/")) 2L else 1L


#' Render a random variable as the expression that built it
#'
#' Leaves show their label; derived variables show their operands joined by the
#' operator. Brackets appear only where they change the reading: around a weaker
#' operand, and around an equally strong right operand of `-` or `/`, neither of
#' which associates.
#' @keywords internal
#' @noRd
rv_expression <- function(x) {
  if (!S7_inherits(x, random_variable)) {
    return(format(x, digits = 7L))
  }
  if (is.na(x@op)) {
    return(short_label(x@label))
  }
  here <- op_precedence(x@op)
  side <- function(operand, right) {
    text <- rv_expression(operand)
    weaker <- S7_inherits(operand, random_variable) &&
      !is.na(operand@op) &&
      (op_precedence(operand@op) < here ||
        (right && op_precedence(operand@op) == here && x@op %in% c("-", "/")))
    if (weaker) paste0("(", text, ")") else text
  }
  paste(
    side(x@operands[[1L]], right = FALSE),
    x@op,
    side(x@operands[[2L]], right = TRUE)
  )
}


#' A short description of a sample space
#' @keywords internal
#' @noRd
family_label <- function(family) {
  name <- attr(S7_class(family), "name")
  sprintf("%s, dimension %d", name, outcome_dim(family))
}


#' @param x A [random_variable].
#' @param ... Ignored.
#' @rdname random_variable
#' @export
method(print, random_variable) <- function(x, ...) {
  cat("<random_variable>", rv_expression(x), "\n")
  cat("  on", family_label(x@family), "\n")
  invisible(x)
}


# --- Mixture Likelihood RV ----------------------------------------------------

#' The mixture likelihood, as a random variable
#'
#' For a mixing measure W, \eqn{X(x) = P_W(x)}{X(x) = P_W(x)}. A likelihood
#' ratio can be defines as a quotient of two such random variables.
#' @param dist An [outcome_distribution].
#' @param label How to name it when printing. `NULL` takes how the argument was
#'   written, so `mixture_likelihood(Q)` prints as `Q`. This degenerates when
#'   the call is behind a helper or inside a loop, since it records how the
#'   variable was written rather than what it is. Name it yourself if that
#'   matters to you.
#' @return A [random_variable].
#' @seealso [random_variable_arithmetic]
#' @examples
#' fam <- gaussian_family(dim = 2)
#' Q <- mixture_likelihood(mixture(point_mixing(c(0.5, 0.5)), fam))
#' Q(c(2,2))
#' @export
mixture_likelihood <- function(dist, label = NULL) {
  if (is.null(label)) {
    label <- deparse1(substitute(dist))
  }
  if (!is.character(label) || length(label) != 1L || is.na(label)) {
    stop("`label` must be a single string, or NULL.", call. = FALSE)
  }
  force(dist)
  random_variable(
    function(x) exp(dist_log_density(dist, x)),
    family = dist@family,
    label = label
  )
}


# --- Arithmetic ---------------------------------------------------------------

#' Arithmetic on random variables and scalars
#'
#' In the following document `X`, `Y` and `Z` are random variables; `a` and `b`
#' are length-1 numerics.
#'
#' Random variables may be scaled, transformed and combined by the arithmetic
#' operations `+`, `-`, `*` and `/`. For instance, we may define a random
#' variable `Y` from `X` by the linear transformation `Y <- 2 * X + 3`.
#' Similarly `Z <- X + Y` is the variable defined by `Z(x) = X(x) + Y(x)`, and
#' likewise for `-`, `*` and `/`. Comparison and other operators are not
#' currently supported.
#'
#' If both operands are random variables, they must be defined on the same sample
#' space, checked here rather than at evaluation, so a mismatch is reported where
#' it was written.
#'
#' @param e1,e2 A [random_variable] or a single number, at least one of them a
#'   random variable.
#' @return A [random_variable].
#' @name random_variable_arithmetic
NULL


#' Both operands must live on the same sample space
#' @keywords internal
#' @noRd
shared_family <- function(e1, e2) {
  if (!S7_inherits(e1, random_variable)) {
    return(e2@family)
  }
  if (!S7_inherits(e2, random_variable)) {
    return(e1@family)
  }
  if (!identical(e1@family, e2@family)) {
    stop(
      "random variables are defined on different sample spaces, so they ",
      "cannot be combined.",
      call. = FALSE
    )
  }
  e1@family
}


#' Build the derived variable for a binary operator
#'
#' Evaluates each operand in turn, so the input is checked once per operand
#' rather than once. That is a few microseconds against a mixture density, and
#' it keeps every variable independently valid rather than trusting a caller.
#' @keywords internal
#' @noRd
combine_rv <- function(e1, e2, op, symbol) {
  family <- shared_family(e1, e2)
  force(op)
  left <- if (S7_inherits(e1, random_variable)) e1 else NULL
  right <- if (S7_inherits(e2, random_variable)) e2 else NULL
  const_left <- if (is.null(left)) as.numeric(e1) else NULL
  const_right <- if (is.null(right)) as.numeric(e2) else NULL

  if (!is.null(const_left) && length(const_left) != 1L) {
    stop(
      "only a single number may be combined with a random variable.",
      call. = FALSE
    )
  }
  if (!is.null(const_right) && length(const_right) != 1L) {
    stop(
      "only a single number may be combined with a random variable.",
      call. = FALSE
    )
  }

  random_variable(
    function(x) {
      force(x)
      a <- if (is.null(left)) const_left else left(x)
      b <- if (is.null(right)) const_right else right(x)
      op(a, b)
    },
    family = family,
    op = symbol,
    operands = list(e1, e2)
  )
}


#' @rdname random_variable_arithmetic
method(`+`, list(random_variable, random_variable)) <- function(e1, e2) {
  combine_rv(e1, e2, `+`, "+")
}
#' @rdname random_variable_arithmetic
method(`+`, list(random_variable, class_numeric)) <- function(e1, e2) {
  combine_rv(e1, e2, `+`, "+")
}
#' @rdname random_variable_arithmetic
method(`+`, list(class_numeric, random_variable)) <- function(e1, e2) {
  combine_rv(e1, e2, `+`, "+")
}

#' @rdname random_variable_arithmetic
method(`-`, list(random_variable, random_variable)) <- function(e1, e2) {
  combine_rv(e1, e2, `-`, "-")
}
#' @rdname random_variable_arithmetic
method(`-`, list(random_variable, class_numeric)) <- function(e1, e2) {
  combine_rv(e1, e2, `-`, "-")
}
#' @rdname random_variable_arithmetic
method(`-`, list(class_numeric, random_variable)) <- function(e1, e2) {
  combine_rv(e1, e2, `-`, "-")
}

#' @rdname random_variable_arithmetic
method(`*`, list(random_variable, random_variable)) <- function(e1, e2) {
  combine_rv(e1, e2, `*`, "*")
}
#' @rdname random_variable_arithmetic
method(`*`, list(random_variable, class_numeric)) <- function(e1, e2) {
  combine_rv(e1, e2, `*`, "*")
}
#' @rdname random_variable_arithmetic
method(`*`, list(class_numeric, random_variable)) <- function(e1, e2) {
  combine_rv(e1, e2, `*`, "*")
}

#' @rdname random_variable_arithmetic
method(`/`, list(random_variable, random_variable)) <- function(e1, e2) {
  combine_rv(e1, e2, `/`, "/")
}
#' @rdname random_variable_arithmetic
method(`/`, list(random_variable, class_numeric)) <- function(e1, e2) {
  combine_rv(e1, e2, `/`, "/")
}
#' @rdname random_variable_arithmetic
method(`/`, list(class_numeric, random_variable)) <- function(e1, e2) {
  combine_rv(e1, e2, `/`, "/")
}
