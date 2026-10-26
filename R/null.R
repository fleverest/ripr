#' @include region.R family.R
NULL

# --- The null hypothesis ------------------------------------------------------

#' A null hypothesis: a family together with its parameter region
#'
#' \eqn{H_0 = \{P_\theta : \theta \in \bigcup_i \Theta_{0i}\}}{H_0 = {P_theta : theta in union_i Theta_0i}}.
#' Each convex part of the region (see [parts()]) is itself a null.
#' The decomposition into [cells()] is computed once, at construction. Test
#' membership with `contains(null@region, theta)`.
#'
#' A part that reaches outside the family's parameter space is intersected
#' with it, so a [halfspace_region()] null for a multinomial becomes the
#' bounded polytope it describes within the simplex.
#'
#' @param family A [parametric_family].
#' @param region The null's geometry: any [region], or a list of
#'   [convex_region]s, which becomes their [union_region].
#' @param max_cells The number of simplices the decomposition may produce
#'   before giving up. Defaults to the `ripr.max_cells` option, or `1000L`
#'   when that is unset.
#' @return A `null_model`.
#' @section Properties:
#' \describe{
#'   \item{`family`}{The [parametric_family].}
#'   \item{`region`}{The null's [region], as declared, with any part reaching
#'   outside the family's parameter space intersected with it.}
#'   \item{`cells`}{The flat list of convex cells the region decomposes into,
#'   in part order: `cells(parts(region)[[1]])`, then those of part 2, and so
#'   on.}
#'   \item{`cell_part`}{Which part each cell came from, as an index into
#'   `parts(region)`; used to report the `part` of each atom in a fit.}
#' }
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' # The plurality null
#' null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' null <- null_model(fam, list(halfspace_region(c(1, -1, 0))))
#' contains(null@region, c(0.2, 0.5, 0.3))
#' contains(null@region, c(0.6, 0.2, 0.2))
#' @export
null_model <- new_class(
  "null_model",
  properties = list(
    family = parametric_family,
    region = region,
    cells = class_list,
    cell_part = class_integer
  ),
  constructor = function(
    family,
    region,
    max_cells = getOption("ripr.max_cells", 1000L)
  ) {
    region <- clip_to_parameters(as_region(region), family@parameter_space)
    if (S7_inherits(region, empty_region)) {
      stop(
        "the null is empty: there is no hypothesis to fit or certify.",
        call. = FALSE
      )
    }
    prts <- parts(region)
    budget <- cell_budget(max_cells)
    per_part <- lapply(seq_along(prts), function(i) {
      tryCatch(
        cells(prts[[i]], .budget = budget),
        error = function(e) {
          stop(
            "could not decompose part ",
            i,
            " of the null: ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )
    })
    new_object(
      S7_object(),
      family = family,
      region = region,
      cells = unlist(per_part, recursive = FALSE),
      cell_part = rep(seq_along(per_part), lengths(per_part))
    )
  },
  validator = function(self) {
    d <- space_dim(self@family@parameter_space)
    d_region <- space_dim(self@region)
    if (d_region != d) {
      return(paste0(
        "the null's region must have dimension ",
        d,
        ", matching the family's parameter space; got ",
        d_region
      ))
    }
    # Catches a property replaced after construction, which would file a cell
    # under the wrong part.
    if (length(self@cells) != length(self@cell_part)) {
      return("`cell_part` must have one entry per element of `cells`")
    }
    NULL
  }
)


#' Intersect each part of `region` that leaves `parameter_space` with it
#'
#' Parts already inside are kept as declared. A part with nothing inside is
#' dropped with a warning, and a region with no part inside is an error. A
#' dimension mismatch is left for the validator.
#' @keywords internal
#' @noRd
clip_to_parameters <- function(region, parameter_space) {
  if (
    S7_inherits(region, empty_region) ||
      space_dim(region) != space_dim(parameter_space)
  ) {
    return(region)
  }
  prts <- parts(region)
  inside <- vapply(prts, within_space, logical(1L), space = parameter_space)
  if (all(inside)) {
    return(region)
  }
  clipped <- prts
  clipped[!inside] <- lapply(prts[!inside], function(p) p & parameter_space)
  empty <- vapply(clipped, is_empty, logical(1L))
  if (all(empty)) {
    stop(
      "the null lies outside the family's parameter space, where the family ",
      "has no distributions: there is no hypothesis to fit or certify.",
      call. = FALSE
    )
  }
  if (any(empty)) {
    rlang::warn(paste0(
      ngettext(sum(empty), "Part ", "Parts "),
      toString(which(empty)),
      " of the null ",
      ngettext(sum(empty), "lies", "lie"),
      " outside the family's parameter space, and ",
      ngettext(sum(empty), "was", "were"),
      " dropped."
    ))
  }
  kept <- clipped[!empty]
  if (length(kept) == 1L) kept[[1L]] else union_region(kept)
}


#' Is a convex region inside `space`?
#'
#' A bounded part is inside when each vertex is, up to [contains()]'s
#' tolerance, so vertices typed in decimal (`c(0.9, 0.05, 0.05)` does not sum
#' to exactly one in binary) are not clipped. Anything with rays or lines is
#' decided exactly.
#' @keywords internal
#' @noRd
within_space <- function(part, space) {
  if (!S7_inherits(part, polyhedron_region)) {
    return(region_subset(part, space))
  }
  g <- part@generators
  if (nrow(g$r) == 0L && nrow(g$l) == 0L) {
    return(all(apply(g$v, 1L, function(v) contains(space, v))))
  }
  region_subset(part, space)
}


#' @rdname null_model
#' @usage NULL
#' @export
method(print, null_model) <- function(x, ...) {
  prts <- parts(x@region)
  n <- length(prts)
  cat("<", class_name(x), ">\n", sep = "")
  cat("  ", format(x@family), "\n", sep = "")
  cat(
    "  ",
    parts_label(n),
    ", ",
    length(x@cells),
    ngettext(length(x@cells), " cell", " cells"),
    "\n",
    sep = ""
  )
  cat_parts(prts)
  invisible(x)
}


#' @description `format()` gives the family and part count on one line.
#' @rdname null_model
#' @usage NULL
#' @export
method(format, null_model) <- function(x, ...) {
  sprintf(
    "%s: %s over %s",
    class_name(x),
    class_name(x@family),
    parts_label(n_parts(x@region))
  )
}
