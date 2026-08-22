#' @include sample_space.R distribution.R
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
#' by [validate_outcome()] first.
#'
#' `Inf` can be a legitimate value here. A likelihood ratio can be genuinely
#' infinite when the denominator is not absolutely continuous with respect to
#' the numerator.
#'
#' @param f The mapping that defines the random variable, accepting an `(n, d)`
#'   matrix and returning `n` numbers.
#' @param sample_space The [sample_space] this variable is defined on.
#' @param label How to name this variable when printing. Ignored when `op` is
#'   given, since the expression is then built from the operands.
#' @param op The operator that produced this variable, or `NA` for a leaf. Set
#'   by [random_variable_arithmetic]; there is rarely a reason to pass it.
#' @param operands The operands `op` combined, or an empty list for a leaf. Set
#'   alongside `op`, and used only for printing.
#' @return A callable `random_variable`.
#' @seealso [likelihood()], [random_variable_arithmetic]
#' @examples
#' X <- random_variable(\(x) dnorm(x, 1), sample_space = real_space(1))
#' X(as.matrix(0:2))
#' @export
random_variable <- new_class(
  "random_variable",
  parent = class_function,
  properties = list(
    sample_space = sample_space,
    label = class_character,
    op = class_character,
    operands = class_list
  ),
  constructor = function(
    f,
    sample_space,
    label = "<rv>",
    op = NA_character_,
    operands = list()
  ) {
    # Forced so the closure captures values, not promises: `saveRDS` on a
    # `class_function` parent serialises whatever the environment holds.
    force(f)
    force(sample_space)
    if (!is.function(f)) {
      stop("`f` must be a function.", call. = FALSE)
    }

    new_object(
      function(x) {
        force(x)
        out <- f(validate_outcome(sample_space, x))
        if (!is.numeric(out)) {
          stop("a random variable must return numbers.", call. = FALSE)
        }
        as.vector(out)
      },
      sample_space = sample_space,
      label = label,
      op = op,
      operands = operands
    )
  }
)


# --- Printing -----------------------------------------------------------------

#' Shorten a label that came from deparsing an inline argument
#'
#' `likelihood(Q)` gives a short label, but an expression written in
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


#' @rdname random_variable
#' @usage NULL
#' @export
method(print, random_variable) <- function(x, ...) {
  cat("<random_variable>", format(x), "\n")
  cat("  on", space_label(x@sample_space), "\n")
  invisible(x)
}


#' @description `format()` gives the expression alone, without the class
#'   banner `print()` adds.
#'
#' Needed rather than inherited: the parent is `class_function`, so
#' `format.default()` reaches `deparse()`. `rv_expression()` already calls
#' `format()` on non-`random_variable` operands, which is how `X / 4.27`
#' renders its divisor, so the generic has to work on these too.
#' @rdname random_variable
#' @usage NULL
#' @export
method(format, random_variable) <- function(x, ...) rv_expression(x)


# --- Likelihood RV ------------------------------------------------------------

#' The likelihood of a distribution, as a random variable
#'
#' \eqn{X(x) = P(x)}{X(x) = P(x)}. A likelihood ratio is then a quotient of two
#' of these: `likelihood(Q) / likelihood(P_star)`.
#' @param dist A [distribution].
#' @param label How to name it when printing. `NULL` takes how the argument was
#'   written, so `likelihood(Q)` prints as `Q`. This degenerates when
#'   the call is behind a helper or inside a loop, since it records how the
#'   variable was written rather than what it is. Name it yourself if that
#'   matters to you.
#' @return A [random_variable].
#' @seealso [random_variable_arithmetic]
#' @examples
#' fam <- gaussian_family(dim = 2)
#' Q <- likelihood(induced_distribution(fam, point_mixing(c(0.5, 0.5))))
#' Q(c(2,2))
#' @export
likelihood <- function(dist, label = NULL) {
  if (is.null(label)) {
    label <- deparse1(substitute(dist))
  }
  if (!is.character(label) || length(label) != 1L || is.na(label)) {
    stop("`label` must be a single string, or NULL.", call. = FALSE)
  }
  force(dist)
  random_variable(
    function(x) exp(log_density(dist, x)),
    sample_space = dist@sample_space,
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
#' @examples
#' fam <- gaussian_family(dim = 1)
#' X <- random_variable(\(x) dnorm(x, 1), sample_space = fam@sample_space)
#' Y <- 2 * X + 3
#' Y(as.matrix(0:2))
#' Z <- X / X
#' Z(as.matrix(0:2))
#' @name random_variable_arithmetic
NULL


#' Both operands must live on the same sample space
#'
#' Compared by value, not identity: two separately built `count_space(20, 3)`
#' objects are `identical()`, so variables from unrelated families over the same
#' space combine freely.
#' @keywords internal
#' @noRd
shared_space <- function(e1, e2) {
  if (!S7_inherits(e1, random_variable)) {
    return(e2@sample_space)
  }
  if (!S7_inherits(e2, random_variable)) {
    return(e1@sample_space)
  }
  if (!identical(e1@sample_space, e2@sample_space)) {
    stop(
      "random variables are defined on different sample spaces, so they ",
      "cannot be combined.",
      call. = FALSE
    )
  }
  e1@sample_space
}


#' Build the derived variable for a binary operator
#'
#' Evaluates each operand in turn, so the input is checked once per operand
#' rather than once. That is a few microseconds against a mixture density, and
#' it keeps every variable independently valid rather than trusting a caller.
#' @keywords internal
#' @noRd
combine_rv <- function(e1, e2, symbol) {
  space <- shared_space(e1, e2)
  op <- get(symbol, envir = baseenv())
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
    sample_space = space,
    op = symbol,
    operands = list(e1, e2)
  )
}


#' @rdname random_variable_arithmetic
#' @usage NULL
method(`+`, list(random_variable, random_variable)) <- function(e1, e2) {
  combine_rv(e1, e2, "+")
}
#' @rdname random_variable_arithmetic
#' @usage NULL
method(`+`, list(random_variable, class_numeric)) <- function(e1, e2) {
  combine_rv(e1, e2, "+")
}
#' @rdname random_variable_arithmetic
#' @usage NULL
method(`+`, list(class_numeric, random_variable)) <- function(e1, e2) {
  combine_rv(e1, e2, "+")
}

#' @rdname random_variable_arithmetic
#' @usage NULL
method(`-`, list(random_variable, random_variable)) <- function(e1, e2) {
  combine_rv(e1, e2, "-")
}
#' @rdname random_variable_arithmetic
#' @usage NULL
method(`-`, list(random_variable, class_numeric)) <- function(e1, e2) {
  combine_rv(e1, e2, "-")
}
#' @rdname random_variable_arithmetic
#' @usage NULL
method(`-`, list(class_numeric, random_variable)) <- function(e1, e2) {
  combine_rv(e1, e2, "-")
}

#' @rdname random_variable_arithmetic
#' @usage NULL
method(`*`, list(random_variable, random_variable)) <- function(e1, e2) {
  combine_rv(e1, e2, "*")
}
#' @rdname random_variable_arithmetic
#' @usage NULL
method(`*`, list(random_variable, class_numeric)) <- function(e1, e2) {
  combine_rv(e1, e2, "*")
}
#' @rdname random_variable_arithmetic
#' @usage NULL
method(`*`, list(class_numeric, random_variable)) <- function(e1, e2) {
  combine_rv(e1, e2, "*")
}

#' @rdname random_variable_arithmetic
#' @usage NULL
method(`/`, list(random_variable, random_variable)) <- function(e1, e2) {
  combine_rv(e1, e2, "/")
}
#' @rdname random_variable_arithmetic
#' @usage NULL
method(`/`, list(random_variable, class_numeric)) <- function(e1, e2) {
  combine_rv(e1, e2, "/")
}
#' @rdname random_variable_arithmetic
#' @usage NULL
method(`/`, list(class_numeric, random_variable)) <- function(e1, e2) {
  combine_rv(e1, e2, "/")
}
