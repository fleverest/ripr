#' @include enclosure.R bernstein.R multinomial.R point_region.R
NULL

# Our bounding methods: exact evaluation for a point null on any family
# with an enumerable sample space, and the Bernstein enclosure for multinomial
# expectations over simplices inside the standard simplex. Each extends the
# `enclosure` class, and is registered via `bounding_method()` on the
# (family, cell) pairs it covers.

# --- Choosing a method --------------------------------------------------------

#' Nothing is wrong with the cell: no method covers this family over it
#' @keywords internal
#' @noRd
no_method_message <- function(family, cell) {
  geometry <- class_name(cell)
  paste0(
    "No bounding method is implemented for ",
    class_name(family),
    " expectations over ",
    geometry,
    ".\n",
    "Certifying this requires deriving and implementing one. ",
    "In the meantime, `search_sup()` can search and approximate the least ",
    "upper bound."
  )
}


method(bounding_method, list(parametric_family, region)) <- function(
  family,
  cell
) {
  no_method_message(family, cell)
}


# Default to evaluation at a point null for enumerable families
method(bounding_method, list(parametric_family, point_region)) <- function(
  family,
  cell
) {
  if (!is_finite_space(family@sample_space)) {
    return(refusal_message(
      "Evaluation at a point",
      point_unenumerable(family),
      family,
      cell
    ))
  }
  if (contains(family@parameter_space, cell@theta)) {
    return(point_enclosure())
  }
  no_method_message(family, cell)
}

# Bernstein for multinomial
method(bounding_method, list(multinomial_family, convex_region)) <- function(
  family,
  cell
) {
  why <- bernstein_obstruction(cell)
  if (is.null(why)) {
    return(bernstein_enclosure())
  }
  refusal_message("The Bernstein enclosure", why, family, cell)
}


# --- Refusals -----------------------------------------------------------------

#' "<subject> cannot bound <family> expectations over this <cell>, because
#' <why$because>." followed by `why$remedy`
#' @keywords internal
#' @noRd
refusal_message <- function(subject, why, family, cell) {
  paste0(
    subject,
    " cannot bound ",
    class_name(family),
    " expectations over this ",
    class_name(cell),
    ", because ",
    why$because,
    ".\n",
    why$remedy
  )
}


#' Why evaluation at a point cannot certify under this family
#' @keywords internal
#'@noRd
point_unenumerable <- function(family) {
  list(
    because = paste0(
      "its expectation is an integral over a `",
      class_name(family@sample_space),
      "` rather than a sum over an enumerable one"
    ),
    remedy = paste0(
      "The expectation has to be computed exactly (i.e. not via quadrature or ",
      "Monte Carlo estimates) for certification. `search_sup()` can still ",
      "approximate a least upper bound."
    )
  )
}


#' Why the Bernstein enclosure cannot handle a region, or `NULL` if it can
#'
#' `NULL` when `reparametrise_to()` will accept the vertices, i.e. it is
#' inside the standard simplex. Returns `list(because, remedy)`, where `because`
#' completes "... cannot bound this region, because ...", `remedy` is a complete
#' sentence.
#' @keywords internal
#' @noRd
bernstein_obstruction <- function(space) {
  only_cert_affected_msg <- paste0(
    "`search_sup()` can still approximate a least upper bound."
  )
  # Before the class tests: it applies to any unbounded region.
  if (!is_bounded(space)) {
    return(list(
      because = paste0(
        "it is unbounded, and no finite set of simplices covers an unbounded ",
        "region"
      ),
      remedy = paste0(
        "The Bernstein enclosure only works for bounded polytopes. State the ",
        "null over a bounded region instead. ",
        only_cert_affected_msg
      )
    ))
  }
  not_simplex <- list(
    because = paste0("it is a `", class_name(space), "` rather than a simplex"),
    remedy = paste0(
      "This is a bug in ripr: `cells()` should have triangulated every ",
      "bounded region into simplices before getting here. Please report it at ",
      "<https://github.com/fleverest/ripr/issues>, including how the reproduce ",
      "the error. ",
      only_cert_affected_msg
    )
  )
  if (!S7_inherits(space, polytope_region)) {
    return(not_simplex)
  }
  v <- space@vertices
  departure <- simplex_departure(v)
  if (!is.null(departure)) {
    return(list(
      because = paste0("its vertices leave the standard simplex: ", departure),
      remedy = paste0(
        "A region that leaves the standard simplex cannot be bounded by ",
        "its Bernstein coefficients, as the Bernstein basis polynomials take ",
        "take negative values there."
      )
    ))
  }
  if (ncol(v) < 2L) {
    return(list(
      because = "it has a single coordinate, so its simplex is one point",
      remedy = paste0(
        "A one-category multinomial has a single outcome, and its expectation ",
        "is trivial; there is nothing to bound."
      )
    ))
  }
  if (!S7_inherits(space, simplex_region)) {
    return(not_simplex)
  }
  NULL
}


# --- Evaluation at a point ---------------------------------------------------

#' Exact evaluation at a point
#'
#' The supremum over `{theta_0}` is `E_theta0[X]`
#' @keywords internal
#' @noRd
point_enclosure <- new_class(
  "point_enclosure",
  parent = enclosure,
  properties = list(values = class_numeric, loglik = class_any),
  constructor = function(values = numeric(), loglik = NULL) {
    new_object(
      S7_object(),
      label = "point",
      values = values,
      loglik = loglik
    )
  }
)

method(prepare_enclosure, point_enclosure) <- function(
  enclosure,
  x,
  family,
  control
) {
  point_enclosure(
    values = evaluate_on_space(x, family),
    loglik = compile_loglik(family, enumerate_space(family@sample_space))
  )
}

#' `E_theta[X]` exactly, from a prepared point enclosure
#' @keywords internal
#' @noRd
point_value <- function(enc, theta) {
  log_p <- as.vector(enc@loglik(matrix(theta, nrow = 1L)))
  # `exp(-Inf) * x` is a clean 0 for finite `x`.
  sum(exp(log_p) * enc@values)
}

method(enclosure_root, point_enclosure) <- function(enclosure, cell) {
  list(
    V = matrix(cell@theta, nrow = 1L),
    theta = cell@theta,
    value = point_value(enclosure, cell@theta)
  )
}

method(enclosure_upper, point_enclosure) <- function(enclosure, box) {
  box$value
}

method(enclosure_best, point_enclosure) <- function(enclosure, box) {
  list(value = box$value, theta = box$theta)
}

# Never reached (hopefully)
method(enclosure_split, point_enclosure) <- function(enclosure, box, edge) {
  stop("a point cannot be split.", call. = FALSE)
}


# --- The Bernstein enclosure --------------------------------------------------

#' The Bernstein enclosure, for multinomial expectations over simplices
#'
#' The multinomial pmf is the Bernstein basis, so `x` on the sample space is
#' already the coefficient vector over the standard simplex; every cell
#' shares that one lattice and evaluation. A box carries its coefficients,
#' whose maximum bounds the polynomial and whose vertex entries are attained.
#' @keywords internal
#' @noRd
bernstein_enclosure <- new_class(
  "bernstein_enclosure",
  parent = enclosure,
  properties = list(lattice = class_any, values = class_numeric),
  constructor = function(lattice = NULL, values = numeric()) {
    new_object(
      S7_object(),
      label = "bernstein",
      lattice = lattice,
      values = values
    )
  }
)

method(prepare_enclosure, bernstein_enclosure) <- function(
  enclosure,
  x,
  family,
  control
) {
  # Check feasibility
  check_bernstein_size(family@n_trials, family@k, control$max_coefficients)
  # Make lattice and get realised values of x
  bernstein_enclosure(
    lattice = bernstein_lattice(family@n_trials, family@k),
    values = evaluate_on_space(x, family)
  )
}

method(enclosure_root, bernstein_enclosure) <- function(enclosure, cell) {
  lat <- enclosure@lattice
  V <- pad_vertices(cell@vertices, lat$K)
  list(V = V, coef = reparametrise_to(enclosure@values, lat, V))
}

# Convex hull property (PBP 10.2, 10.3 Remark 2).
method(enclosure_upper, bernstein_enclosure) <- function(enclosure, box) {
  box_bound(box)
}

# Vertex coefficients are exact values (PBP 10.2).
method(enclosure_best, bernstein_enclosure) <- function(enclosure, box) {
  box_best(box, enclosure@lattice)
}

method(enclosure_split, bernstein_enclosure) <- function(enclosure, box, edge) {
  bisect(box, edge[[1L]], edge[[2L]], enclosure@lattice)
}


# --- Shared -------------------------------------------------------------------

#' `x` on the whole enumerable sample space; errors unless finite everywhere
#' @keywords internal
#' @noRd
evaluate_on_space <- function(x, family) {
  values <- x(enumerate_space(family@sample_space))
  if (any(!is.finite(values))) {
    stop(
      "Cannot certify: the variable is not finite everywhere on the sample ",
      "space, so its null expectation is unbounded.",
      call. = FALSE
    )
  }
  values
}
