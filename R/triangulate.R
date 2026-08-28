#' @include region.R
NULL

# --- Triangulation ------------------------------------------------------------
#
#  In this file we implement vertex-fan triangulation. Pick a vertex `v` of a
# polytope `P`. Every facet `F` of `P` that does not contain `v` cones back to
# it, and those cones tile `P`:
#
#   P = union over facets F with v not in F of conv({v} u F).
#
# `F` is itself a polytope, one dimension down, so the same construction inside
# its affine hull triangulates it, and coning a simplex of `F` to `v` gives a
# simplex of `P`.
# The recursion bottoms out where the vertex count is one more than the
# dimension, which is a simplex already.

#' Break a bounded region into simplices
#'
#' The vertex fan described above, returned as [simplex_region()] cells whose
#' union is `space` and whose interiors are disjoint. A simplex triangulates to
#' itself.
#'
#' @param space A bounded [convex_region].
#' @param max_cells Give up rather than produce more simplices than this. The
#'   fan is combinatorial in the vertex count at worst, and a decomposition
#'   past this size is probably not something we can afford to use downstream.
#' @return A list of [simplex_region] objects.
#' @keywords internal
#' @noRd
triangulate <- function(space, max_cells = 1000L) {
  if (!is_bounded(space)) {
    stop(
      "only a bounded region can be triangulated; this `",
      class_name(space),
      "` has rays or lineality, and no finite set of simplices covers an ",
      "unbounded region.",
      call. = FALSE
    )
  }
  budget <- new.env(parent = emptyenv())
  budget$left <- max_cells
  budget$max_cells <- max_cells
  lapply(fan_cells(q_nonredundant(q_vrep(space)), budget), simplex_from_qv)
}


#' One level of the fan, on rational V-representations
#'
#' Takes the vertices of a polytope and returns the vertex sets of its
#' simplices, each a row subset of the input, so no coordinate is ever
#' recomputed.
#' @param qv An rcdd V-representation of a bounded polyhedron, all rows
#'   extreme points.
#' @param budget An environment with `left`, the simplices still allowed, and
#'   `max_cells`, the original cap for the error message. Every base-case
#'   emission maps one-to-one onto a final cell -- coning preserves the count
#'   -- so decrementing there counts the finished triangulation exactly.
#' @return A list of rcdd V-representations, one per simplex.
#' @keywords internal
#' @noRd
fan_cells <- function(qv, budget) {
  facets <- q_facets(qv)
  # `scdd()` states the affine hull as equality rows, the ambient dimension
  # left over is the dimension of the hull itself.
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
  # At some point in the future we could add a rule for picking an apex
  # but for now the first one is fine.
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


#' Build one simplex from an exact rational V-representation
#'
#' @keywords internal
#' @noRd
simplex_from_qv <- function(qv) {
  # As in `region_from_qh()`: the cell keeps the exact representations, so
  # algebra or a further decomposition on it starts where this one left off
  # rather than from the rounding.
  simplex_region(.hv = hv_from_qv(qv))
}


# --- cells() ------------------------------------------------------------------

#' @description A polyhedron's cells are its triangulation: the fan of
#'   simplices described in [polytope_region()]. Boundedness is what decides
#'   it, not the class -- an unbounded polyhedron has no simplicial
#'   decomposition and is its own only cell, as is one that is a simplex
#'   already: [simplex_region()], or [point_region()], the degenerate one.
#' @rdname cells
#' @usage NULL
method(cells, polyhedron_region) <- function(space, ..., max_cells = 1000L) {
  if (is_bounded(space)) triangulate(space, max_cells) else list(space)
}


#' @rdname cells
#' @usage NULL
method(cells, simplex_region) <- function(space, ...) list(space)


#' @rdname cells
#' @usage NULL
method(cells, point_region) <- function(space, ...) list(space)
