#' @include space.R distribution.R
NULL

# Random variables on a sample space, and arithmetic over them.

#' A random variable on a sample space
#'
#' A callable `random_variable` maps a length-`d` outcome to a number, or an
#' `(n, d)` matrix of outcomes to `n` numbers, after checking they belong to
#' its sample space. `Inf` is a legitimate value, e.g. a likelihood ratio whose
#' denominator is zero.
#'
#' @param f The mapping that defines the random variable, accepting an `(n, d)`
#'   matrix and returning `n` numbers.
#' @param sample_space The [space] this variable is defined on.
#' @param label How to name this variable when printing. Ignored when `op` is
#'   given, since the expression is then built from the operands.
#' @param op The operator that produced this variable, or `NA` for a leaf. Set
#'   by [random_variable_arithmetic]; unnecessary externally.
#' @param operands The operands `op` combined, or an empty list for a leaf.
#'   Used only for printing.
#' @param log_f Optional: the same mapping in log space, `log_f(x) ==
#'   log(f(x))`, for a non-negative variable. `NULL` otherwise.
#' @return A callable `random_variable`.
#' @seealso [likelihood()], [random_variable_arithmetic]
#' @examples
#' X <- random_variable(\(x) dnorm(x, 1), sample_space = real_region(1))
#' X(as.matrix(0:2))
#' @export
random_variable <- new_class(
  "random_variable",
  parent = class_function,
  properties = list(
    sample_space = space,
    label = class_character,
    op = class_character,
    operands = class_list,
    log_f = class_any
  ),
  constructor = function(
    f,
    sample_space,
    label = "<rv>",
    op = NA_character_,
    operands = list(),
    log_f = NULL
  ) {
    # Forced so the closure captures values, not promises: `saveRDS` on a
    # `class_function` parent serialises whatever the environment holds.
    force(f)
    force(sample_space)
    force(log_f)
    if (!is.function(f)) {
      stop("`f` must be a function.", call. = FALSE)
    }
    if (!is.null(log_f) && !is.function(log_f)) {
      stop("`log_f` must be a function, or NULL.", call. = FALSE)
    }

    new_object(
      checked_mapping(f, sample_space),
      sample_space = sample_space,
      label = label,
      op = op,
      operands = operands,
      log_f = if (!is.null(log_f)) checked_mapping(log_f, sample_space)
    )
  }
)


#' Wrap `f` to validate its input against `sample_space` and require numbers out.
#' @keywords internal
#' @noRd
checked_mapping <- function(f, sample_space) {
  force(f)
  force(sample_space)
  function(x) {
    force(x)
    out <- f(validate_outcome(sample_space, x))
    if (!is.numeric(out)) {
      stop("a random variable must return numbers.", call. = FALSE)
    }
    as.vector(out)
  }
}


#' Does this variable have a known log form?
#' @keywords internal
#' @noRd
has_log_form <- function(x) !is.null(x@log_f)


# --- Printing -----------------------------------------------------------------

#' Binding strength, for deciding where brackets are needed
#' @keywords internal
#' @noRd
op_precedence <- function(op) if (op %in% c("*", "/")) 2L else 1L


#' Render a random variable as the expression that built it. Brackets only
#' around a weaker operand, or an equally strong right operand of `-` or `/`.
#' @keywords internal
#' @noRd
rv_expression <- function(x) {
  if (!S7_inherits(x, random_variable)) {
    return(format(x, digits = 7L))
  }
  if (is.na(x@op)) {
    return(x@label)
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
method(print, random_variable) <- function(x, ...) {
  cat("<random_variable>", format(x), "\n")
  cat("  on", space_label(x@sample_space), "\n")
  invisible(x)
}


#' @description `format()` gives the expression alone, without the class
#'   banner `print()` adds.
#' @rdname random_variable
#' @usage NULL
method(format, random_variable) <- function(x, ...) rv_expression(x)


# --- Likelihood RV ------------------------------------------------------------

#' The likelihood of a distribution, as a random variable
#'
#' \eqn{X(x) = P(x)}{X(x) = P(x)}. A likelihood ratio is then a quotient of two
#' of these: `likelihood(Q) / likelihood(P_star)`.
#' @param dist A [distribution].
#' @param label How to name it when printing. `NULL` describes the
#'   distribution, e.g. `P[theta = (0.5, 0.5)]` or `P[mixed over 3 atoms]`;
#'   pass a short name such as `"Q"` when printing ratios.
#' @return A [random_variable].
#' @seealso [random_variable_arithmetic]
#' @examples
#' fam <- gaussian_family(d = 2)
#' Q <- likelihood(mixture(fam, dirac(c(0.5, 0.5))))
#' Q
#' Q(c(2, 2))
#' likelihood(fam(c(0.5, 0.5)), label = "Q")
#' @export
likelihood <- function(dist, label = NULL) {
  if (is.null(label)) {
    label <- likelihood_label(dist)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label)) {
    stop("`label` must be a single string, or NULL.", call. = FALSE)
  }
  force(dist)
  random_variable(
    function(x) exp(log_density(dist, x)),
    sample_space = dist@sample_space,
    label = label,
    log_f = function(x) log_density(dist, x)
  )
}


#' The default label for `likelihood()`. A mixture omits its family's name;
#' long parameter vectors are elided so the label fits inside an expression.
#' @keywords internal
#' @noRd
likelihood_label <- function(dist) {
  inner <- if (!S7_inherits(dist, mixture)) {
    format(dist)
  } else {
    n <- n_atoms(dist@mixing)
    if (identical(n, 1L)) {
      theta <- signif(atoms(dist@mixing)[1L, ], 3L)
      if (length(theta) > 4L) {
        theta <- c(theta[1:3], "...")
      }
      paste0("theta = (", toString(theta), ")")
    } else if (is.na(n)) {
      paste("mixed over", class_name(dist@mixing))
    } else {
      paste("mixed over", count_label(n, "atom"))
    }
  }
  paste0("P[", inner, "]")
}


# --- Arithmetic ---------------------------------------------------------------

#' Arithmetic on random variables and scalars
#'
#' Random variables combine with each other and with single numbers under `+`,
#' `-`, `*` and `/`, pointwise: `Z <- X + Y` is `Z(x) = X(x) + Y(x)`, and
#' `2 * X + 3` rescales `X`. Other operators are not supported.
#'
#' Two random variables must share a sample space; a mismatch is an error when
#' the expression is built, not when it is evaluated.
#'
#' @param e1,e2 A [random_variable] or a single number, at least one of them a
#'   random variable.
#' @return A [random_variable].
#' @examples
#' fam <- gaussian_family(d = 1)
#' X <- random_variable(\(x) dnorm(x, 1), sample_space = fam@sample_space)
#' Y <- 2 * X + 3
#' Y(as.matrix(0:2))
#' Z <- X / X
#' Z(as.matrix(0:2))
#' @name random_variable_arithmetic
NULL


#' Same class and properties. Not `identical()`, which fails after a `saveRDS()`
#' round trip.
#' @keywords internal
#' @noRd
same_space <- function(x, y) {
  identical(class(x), class(y)) && identical(S7::props(x), S7::props(y))
}


#' The operands' common sample space, compared by value so separately built
#' equal spaces combine freely.
#' @keywords internal
#' @noRd
shared_space <- function(e1, e2) {
  if (!S7_inherits(e1, random_variable)) {
    return(e2@sample_space)
  }
  if (!S7_inherits(e2, random_variable)) {
    return(e1@sample_space)
  }
  if (!same_space(e1@sample_space, e2@sample_space)) {
    stop(
      "random variables are defined on different sample spaces, so they ",
      "cannot be combined.",
      call. = FALSE
    )
  }
  e1@sample_space
}


#' An operand's log form, or `NULL` if it has none (a negative constant, or a
#' variable without `log_f`). A derived variable has one only if both do.
#' @keywords internal
#' @noRd
log_operand <- function(e) {
  if (S7_inherits(e, random_variable)) {
    return(e@log_f)
  }
  value <- as.numeric(e)
  if (length(value) != 1L || is.na(value) || value < 0) {
    return(NULL)
  }
  function(x) log(value)
}


#' How an operator acts in log space; `NULL` for `-`, whose result may be
#' negative.
#' @keywords internal
#' @noRd
log_combiner <- function(symbol) {
  switch(
    symbol,
    "*" = `+`,
    "/" = `-`,
    "+" = function(a, b) row_logsumexp(cbind(a, b)),
    NULL
  )
}

#' An operand as a function of the outcomes; a constant becomes a scalar.
#' @keywords internal
#' @noRd
value_operand <- function(e) {
  if (S7_inherits(e, random_variable)) {
    return(e)
  }
  value <- as.numeric(e)
  if (length(value) != 1L) {
    stop(
      "only a single number may be combined with a random variable.",
      call. = FALSE
    )
  }
  function(x) value
}


#' Build the derived variable for a binary operator. The log form is carried
#' along, so a likelihood ratio stays finite where `X / Y` is near `0 / 0`.
#' @keywords internal
#' @noRd
combine_rv <- function(e1, e2, symbol) {
  space <- shared_space(e1, e2)
  op <- get(symbol, envir = baseenv())
  left <- value_operand(e1)
  right <- value_operand(e2)

  log_left <- log_operand(e1)
  log_right <- log_operand(e2)
  log_op <- log_combiner(symbol)
  log_f <- if (!is.null(log_left) && !is.null(log_right) && !is.null(log_op)) {
    function(x) log_op(log_left(x), log_right(x))
  }

  random_variable(
    function(x) {
      force(x)
      op(left(x), right(x))
    },
    sample_space = space,
    op = symbol,
    operands = list(e1, e2),
    log_f = log_f
  )
}


# Define random variable arithmetic
local({
  signatures <- list(
    list(random_variable, random_variable),
    list(random_variable, class_numeric),
    list(class_numeric, random_variable)
  )
  for (symbol in c("+", "-", "*", "/")) {
    for (signature in signatures) {
      generic <- get(symbol, envir = baseenv())
      method(generic, signature) <- local({
        symbol <- symbol
        function(e1, e2) combine_rv(e1, e2, symbol)
      })
    }
  }
})
