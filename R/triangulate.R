#' @include region.R
#' @include polyhedron_region.R
#' @include simplex_region.R
#' @include point_region.R
NULL

# --- Triangulation ------------------------------------------------------------
#
# Vertex-fan triangulation: for a vertex `v` of polytope `P`,
#   P = union over facets F with v not in F of conv({v} u F),
# recursing into each `F` until the vertex count is one more than the
# dimension.

#' Break a bounded region into simplices
#'
#' Returns interior-disjoint [simplex_region()] cells covering `space`. Errors
#' past `max_cells` (the fan can be combinatorial in the vertex count); pass
#' `budget` to share one cap across a union's parts.
#' @keywords internal
#' @noRd
triangulate <- function(
  space,
  max_cells = getOption("ripr.max_cells", 1000L),
  budget = NULL
) {
  if (!is_bounded(space)) {
    stop(
      "only a bounded region can be triangulated; this `",
      class_name(space),
      "` has rays or lineality, and no finite set of simplices covers an ",
      "unbounded region.",
      call. = FALSE
    )
  }
  if (is.null(budget)) {
    budget <- cell_budget(max_cells)
  }
  lapply(fan_cells(q_nonredundant(q_vrep(space)), budget), simplex_from_qv)
}


#' One level of the fan, on rational V-representations
#'
#' `qv` must hold only extreme points. Each simplex is a row subset of `qv`, so
#' no coordinate is recomputed. Coning preserves the count, so decrementing
#' `budget$left` at the base case counts final cells exactly.
#' @keywords internal
#' @noRd
fan_cells <- function(qv, budget) {
  facets <- q_facets(qv)
  # `scdd()` states the affine hull as equality rows.
  hull_dim <- (ncol(qv) - 2L) - sum(facets$h[, 1L] == "1")
  if (nrow(qv) == hull_dim + 1L) {
    budget$left <- budget$left - 1L
    if (budget$left < 0L) {
      stop(
        "triangulation gave up after `max_cells = ",
        budget$max_cells,
        "` simplices.",
        call. = FALSE
      )
    }
    return(list(qv))
  }
  apex <- q_subrows(qv, 1L)
  cells <- list()
  for (r in which(facets$h[, 1L] == "0")) {
    on <- facets$on[[r]]
    # A facet through the apex cones to nothing.
    if (length(on) == 0L || 1L %in% on) {
      next
    }
    cells <- c(
      cells,
      lapply(
        fan_cells(q_subrows(qv, on), budget),
        \(cell) q_rbind(apex, cell)
      )
    )
  }
  cells
}


#' Build one simplex from an exact rational V-representation, keeping it exact
#' @keywords internal
#' @noRd
simplex_from_qv <- function(qv) {
  simplex_region(.hv = hv_fill(qv = qv))
}


# --- cells() ------------------------------------------------------------------

#' @description A bounded polyhedron's cells are its triangulation; an
#'   unbounded one, or a [simplex_region()], is its own only cell.
#' @rdname cells
#' @usage NULL
method(cells, polyhedron_region) <- function(
  space,
  max_cells = getOption("ripr.max_cells", 1000L),
  .budget = NULL
) {
  if (is_bounded(space)) triangulate(space, max_cells, .budget) else list(space)
}


#' @rdname cells
#' @usage NULL
method(cells, simplex_region) <- function(
  space,
  max_cells = getOption("ripr.max_cells", 1000L),
  .budget = NULL
) {
  list(space)
}

