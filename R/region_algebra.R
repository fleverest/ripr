#' @include region.R
NULL

# --- Set algebra on regions ---------------------------------------------------

#' Build one region from an exact rational H-representation
#'
#' The only bridge from rational land back to double country. The generators
#' come from one exact double-description step, so a coordinate that is 1/3 in
#' rationals lands on the nearest double once rather than accumulating error
#' over the chain.
#'
#' @keywords internal
#' @noRd
region_from_qh <- function(qh) {
  qh <- q_nonredundant(qh)
  qv <- q_scdd(qh)
  g <- with_origin_vertex(from_vmatrix(qv))
  facets <- from_hmatrix(qh)
  n_v <- ncol(g$v)
  d <- nrow(g$v)

  # Check if the region is a polytope or simplex
  bounded <- ncol(g$r) == 0L && ncol(g$l) == 0L
  simplex <- n_v <= d + 1L && d - sum(facets$eq) == n_v - 1L
  out <- if (bounded) {
    if (simplex) {
      simplex_region(vertices = g$v, facets = facets)
    } else {
      polytope_region(vertices = g$v, facets = facets)
    }
  } else {
    polyhedron_region(
      vertices = g$v,
      rays = g$r,
      lines = g$l,
      facets = facets
    )
  }
  # The cell remembers its exact representations, so that algebra composed on
  # this cell uses the known rational representation rather than rounding.
  out@q_cache <- list(h = qh, v = qv)
  out
}


# --- Union region -------------------------------------------------------------

#' A finite union of convex regions
#'
#' The union \eqn{\bigcup_i \Theta_{0i}}{union_i Theta_0i} of finitely many
#' [convex_region]s, which generally is not convex. A null hypothesis may be one
#' such union, and so may be the support of a truncated prior.
#'
#' A `union_region` is a [region] but deliberately **not** a [convex_region].
#' [chart()], [project()], `maximise_over()` assume convexity, and a union of
#' convex sets is not convex. What this class does implement is [space_dim()],
#' [contains()], [parts()] and [cells()], so that optimisation procedures that
#' require certain properties may operate on the individual components that
#' comply.
#'
#' For instance, `maximise_over` requires a single continuous coordinate system
#' for the entire space, so it runs on a loop over the regions [parts()].
#' `certify` for a [multinomial_family] runs only over [simplex_region]s, so
#' we may compute triangulation accessible via [cells()].
#'
#' Given exactly one convex region, `union_region()` returns it unchanged.
#'
#' @param ... [convex_region] objects, other `union_region` objects, and lists
#'   of either, in any combination and any nesting. A `union_region` argument
#'   flattens rather than nests.
#' @return A `union_region`, or the lone [convex_region] it was given.
#' @section Properties:
#' \describe{
#'   \item{`parts`}{The flat list of convex cells, as declared.}
#' }
#' @examples
#' # The K = 3 plurality null: two overlapping sub-simplices.
#' union_region(
#'   simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'   simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#' )
#'
#' # Nesting is flattened, so these agree:
#' s <- simplex_region(vertices = diag(3))
#' h <- halfspace_region(normal = c(1, -1, 0))
#' n_parts(union_region(s, h))
#' n_parts(union_region(list(s, h)))
#' n_parts(union_region(union_region(s), list(h)))
#'
#' # One cell is already a region, so it is handed back as it came:
#' identical(union_region(s), s)
#' @seealso [region_algebra] for the set-operation verbs: `union()` dispatches
#'   to this constructor, and `intersect()` computes new regions from old.
#' @export
union_region <- new_class(
  "union_region",
  parent = region,
  properties = list(parts = class_list),
  constructor = function(...) {
    flat <- flatten_parts(list(...))
    if (length(flat) == 1L && S7_inherits(flat[[1L]], convex_region)) {
      return(flat[[1L]])
    }
    new_object(S7_object(), parts = flat)
  },
  validator = function(self) {
    if (length(self@parts) == 0L) {
      return("`parts` must be a non-empty list")
    }
    ok <- vapply(
      self@parts,
      \(p) S7_inherits(p, convex_region),
      logical(1)
    )
    if (!all(ok)) {
      return("every element of `parts` must be a `convex_region`")
    }
    # Ambient dimension only; shape and codimension are free. Cells of
    # differing ambient dimension have no common space to union in, and
    # comparing one against a parameter would silently recycle rather than
    # complain. Both dimensions are named, since neither is more wrong.
    dims <- vapply(self@parts, space_dim, integer(1))
    if (length(unique(dims)) > 1L) {
      return(paste0(
        "every element of `parts` must have the same ambient dimension; got ",
        paste(unique(dims), collapse = ", ")
      ))
    }
    NULL
  }
)


#' Flatten union-ish input into a list of convex parts
#'
#' Descends bare lists, unwraps unions into their own parts, and leaves anything
#' else alone as a leaf for the validator to name.
#' @keywords internal
#' @noRd
flatten_parts <- function(x) {
  if (S7_inherits(x, union_region)) {
    return(x@parts)
  }
  if (S7_inherits(x, convex_region)) {
    return(list(x))
  }
  if (is.list(x) && !S7_inherits(x)) {
    return(c(list(), unlist(lapply(x, flatten_parts), recursive = FALSE)))
  }
  list(x)
}


#' Coerce region-ish input to a [region]
#'
#' A [region] passes through untouched; a list becomes a [union_region] of its
#' elements, which for a one-element list is that element itself.
#' @param x A [region], or a list of them.
#' @return A [region].
#' @keywords internal
#' @noRd
as_region <- function(x) {
  if (S7_inherits(x, region)) x else union_region(x)
}


method(space_dim, union_region) <- function(space) {
  # The validator has already established that there is at least one part and
  # that they agree, so the first one speaks for all of them.
  space_dim(space@parts[[1L]])
}


method(contains, union_region) <- function(space, theta, tol = 1e-8) {
  any(vapply(space@parts, \(p) contains(p, theta, tol), logical(1)))
}

#' The refusal both representations owe a union
#' @keywords internal
#' @noRd
refuse_union <- function(what) {
  stop(
    "`",
    what,
    "()` is not defined for a `union_region`: a union is not an ",
    "intersection of half-spaces and has no single generator set. Take the ",
    "representation of each of `parts()` or `cells()` instead.",
    call. = FALSE
  )
}

method(h_rep, union_region) <- function(space) refuse_union("h_rep")
method(v_rep, union_region) <- function(space) refuse_union("v_rep")
method(q_hrep, union_region) <- function(space) refuse_union("q_hrep")
method(q_vrep, union_region) <- function(space) refuse_union("q_vrep")


method(is_empty, union_region) <- function(space) {
  all(vapply(space@parts, is_empty, logical(1)))
}


method(is_bounded, union_region) <- function(space) {
  all(vapply(space@parts, is_bounded, logical(1)))
}


method(parts, union_region) <- function(space) space@parts


#' @description A union's cells are its parts' cells, flattened: the parts are
#'   what was declared, the cells are what the algorithms run on.
#' @rdname cells
#' @usage NULL
method(cells, union_region) <- function(space) {
  unlist(lapply(space@parts, cells), recursive = FALSE)
}


#' The count of cells, as it should read in a message
#' @keywords internal
#' @noRd
parts_label <- function(n) sprintf("%d part%s", n, if (n == 1L) "" else "s")


#' @rdname union_region
#' @usage NULL
#' @export
method(print, union_region) <- function(x, ...) {
  n <- length(x@parts)
  cat("<", attr(S7_class(x), "name"), ">\n", sep = "")
  cat("  ", parts_label(n), ", dimension ", space_dim(x), "\n", sep = "")
  # The cells share a dimension, so the header has already said it. Each
  # part's own format() line while that is readable, a tally beyond it: a
  # triangulated null can hold hundreds of cells, and listing them tells the
  # reader nothing the tally does not.
  if (n <= 6L) {
    for (p in x@parts) {
      cat("    ", format(p), "\n", sep = "")
    }
  } else {
    named <- vapply(x@parts, \(p) attr(S7_class(p), "name"), character(1))
    tally <- table(named)
    for (nm in names(tally)) {
      cat("    ", tally[[nm]], " x ", nm, "\n", sep = "")
    }
  }
  invisible(x)
}


#' @description `format()` gives the same summary on one line, without the class
#'   banner and the per-cell listing that `print()` adds.
#' @rdname union_region
#' @usage NULL
#' @export
method(format, union_region) <- function(x, ...) {
  sprintf(
    "%s: %s, dimension %d",
    attr(S7_class(x), "name"),
    parts_label(length(x@parts)),
    space_dim(x)
  )
}


#' Set algebra on regions
#'
#' Performs set union, intersection, asymmetric difference on two [region]s.
#'
#' Union is structural: a [union_region()] *is* its parts, nothing is
#' computed, and `union()` forwards to its constructor.
#'
#' Intersection is computed precisely. It distributes over union, so the
#' result is the union over every intersection of one part from each argument,
#' where empty intersections are pruned via exact rational feasibility. Parts
#' of a [union_region] may overlap, and so may the parts of the result; nothing
#' here checks to make sure they are disjoint.
#'
#' `setdiff(x, y)` is the closure of `x` minus `y`. Equivalently, the complement
#' of `y` within `x`.  A part of `y` that never meets `x` subtracts nothing. A
#' part `y` that meets `x` in a lower-dimensional slice (e.g. `x` is R^2 and
#' `y` a line segment) also subtracts nothing, becuase we compute only a closed
#' difference and closure reverses the subtraction. Doing so will raise a
#' warning. The decomposition is a product across the parts of `y`, guarded by
#' `max_cells` (default `1000L`, passed through `...`).
#'
#' `setequal(x, y)` decides whether two regions are the same set, exactly.
#' Containment each way is what it tests, and where the containing side is
#' convex that is one exact linear program per facet of it, with nothing
#' decomposed.
#'
#' An empty difference is not what equality means here. The difference is
#' closed, so a region covered exactly by a decomposition of itself still leaves
#' behind the boundaries its pieces share: `setdiff()` of a square and its two
#' triangles is the diagonal between them, where `setequal()` of the same two is
#' `TRUE`.
#'
#' @param x,y [region]s, or anything the base R namesake accepts.
#' @param ... Further regions.
#' @return A [region], except from `setequal()`, which returns `TRUE` or
#'   `FALSE`. `union()` returns what [union_region()] would.
#'   `intersect()` returns a [union_region] of the surviving cells, the lone
#'   cell itself when a single one survives, or `NULL` when the intersection
#'   is empty.
#' @examples
#' # The K = 3 plurality null, by verb rather than constructor:
#' union(
#'   simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'   simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#' )
#'
#' # Its two cells meet in the region where candidate 1 trails both others:
#' intersect(
#'   simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'   simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#' )
#'
#' # Disjoint regions intersect in nothing:
#' intersect(
#'   point_region(theta = c(1, 0, 0)),
#'   point_region(theta = c(0, 1, 0))
#' )
#'
#' # The complement of the K = 3 plurality null within the simplex is the
#' # region where candidate 1 wins -- the alternative, as a region:
#' plurality <- union(
#'   simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'   simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#' )
#' setdiff(simplex_region(vertices = diag(3)), plurality)
#'
#' # A square, and the same square cut into two triangles: different objects,
#' # the same set. One of those triangles alone is not.
#' square <- polytope_region(vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)))
#' setequal(square, union_region(cells(square)))
#' setequal(square, cells(square)[[1]])
#'
#' # Not regions, so base R behaviour, untouched:
#' union(c(1, 2), c(2, 3))
#' intersect(c(1, 2), c(2, 3))
#' setdiff(c(1, 2), c(2, 3))
#' setequal(c(1, 2), c(2, 1))
#' @name region_algebra
NULL


#' @rdname region_algebra
#' @export
union <- function(x, y, ...) UseMethod("union")


#' @rdname region_algebra
#' @export
union.default <- function(x, y, ...) base::union(x, y)


method(union, region) <- function(x, y, ...) union_region(x, y, ...)


#' @rdname region_algebra
#' @export
intersect <- function(x, y, ...) UseMethod("intersect")


#' @rdname region_algebra
#' @export
intersect.default <- function(x, y, ...) base::intersect(x, y)


method(intersect, region) <- function(x, y, ...) {
  regions <- lapply(c(list(x, y), list(...)), as_region)
  dims <- vapply(regions, space_dim, integer(1))
  if (length(unique(dims)) > 1L) {
    stop(
      "every region must have the same ambient dimension; got ",
      paste(unique(dims), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  acc <- lapply(parts(regions[[1L]]), q_hrep)
  for (region_i in regions[-1L]) {
    hs <- lapply(parts(region_i), q_hrep)
    acc <- unlist(
      lapply(acc, \(x) lapply(hs, \(y) q_rbind(x, y))),
      recursive = FALSE
    )
    # Prune before building anything: `prod(n_parts)` combinations will be
    # mostly empty for practical nulls, and every one kept costs a
    # V-representation computation.
    acc <- Filter(Negate(q_is_empty), acc)
    if (length(acc) == 0L) {
      return(NULL)
    }
  }
  union_region(lapply(acc, region_from_qh))
}


#' @rdname region_algebra
#' @export
setdiff <- function(x, y, ...) UseMethod("setdiff")


#' @rdname region_algebra
#' @export
setdiff.default <- function(x, y, ...) base::setdiff(x, y)


method(setdiff, region) <- function(x, y, ..., max_cells = 1000L) {
  y <- as_region(y)
  if (space_dim(x) != space_dim(y)) {
    stop(
      "every region must have the same ambient dimension; got ",
      space_dim(x),
      " and ",
      space_dim(y),
      ".",
      call. = FALSE
    )
  }
  # Difference distributes over the minuend's parts:
  # (A1 u A2) \ y = (A1 \ y) u (A2 \ y).
  results <- lapply(
    parts(x),
    \(ambient) part_difference(q_hrep(ambient), parts(y), max_cells)
  )
  n_sliced <- sum(vapply(results, \(r) r$n_sliced, integer(1)))
  if (n_sliced > 0L) {
    warning(slice_warning(paste0(
      "in ",
      count_label(n_sliced, "case"),
      ", a part of `y` met a part of `x` only in a lower-dimensional ",
      "slice; nothing was subtracted there, since a closed difference ",
      "removes nothing from a slice."
    )))
  }
  cells <- unlist(lapply(results, \(r) r$cells), recursive = FALSE)
  if (length(cells) == 0L) {
    return(NULL)
  }
  union_region(lapply(cells, region_from_qh))
}


#' The warning `setdiff()` raises when a part subtracts nothing
#'
#' Classed, so that an internal caller who is subtracting only to answer a
#' question (`setequal()` asking whether a difference is empty, `disjoin()`
#' peeling parts apart for a measure) can silence just this warning without
#' suppressing other warnings upstream.
#' @keywords internal
#' @noRd
slice_warning <- function(message) {
  structure(
    class = c("ripr_slice_warning", "warning", "condition"),
    list(message = message, call = NULL)
  )
}


#' Run an expression with the slice warning suppressed
#' @keywords internal
#' @noRd
without_slice_warning <- function(expr) {
  withCallingHandlers(
    expr,
    ripr_slice_warning = function(w) invokeRestart("suppressWarning")
  )
}


#' @rdname region_algebra
#' @export
setequal <- function(x, y, ...) UseMethod("setequal")


#' @rdname region_algebra
#' @export
setequal.default <- function(x, y, ...) base::setequal(x, y)


method(setequal, region) <- function(x, y, ..., max_cells = 1000L) {
  y <- as_region(y)
  if (space_dim(x) != space_dim(y)) {
    # Not an error, unlike `intersect()` and `setdiff()`. Those have no answer
    # to give for regions of different ambient dimensions; this one does, and
    # it is that two sets living in different spaces are not the same set.
    return(FALSE)
  }
  region_subset(x, y, max_cells) && region_subset(y, x, max_cells)
}


#' Is every point of one region in another?
#'
#' Checks if `inner` is a subset of `whole`, where both are [region]s.
#' @keywords internal
#' @noRd
region_subset <- function(inner, whole, max_cells = 1000L) {
  # If whole is just one convex_region, we can check with a single
  # call for each part of inner.
  if (S7_inherits(whole, convex_region)) {
    qh <- q_hrep(whole)
    return(all(vapply(
      parts(inner),
      \(p) q_subset(q_hrep(p), qh),
      logical(1)
    )))
  }
  # Otherwise we subtract each convex part of whole from each part of
  # inner, then check the dimension of the remainders.
  subtract <- parts(whole)
  for (p in parts(inner)) {
    qh <- q_hrep(p)
    dim_p <- q_dim(qh)
    leftover <- part_difference(qh, subtract, max_cells)$cells
    if (any(vapply(leftover, \(cell) q_dim(cell) == dim_p, logical(1)))) {
      return(FALSE)
    }
  }
  TRUE
}


# --- Disjoining ---------------------------------------------------------------

#' Transform a region's parts into a disjoint cover of the union
#'
#' Sequential differences: leave the first part as is, then subtract the first
#' from the second, subtract both from the third, and so on. The result covers
#' the same set and its parts meet only on shared boundaries, so a measure can
#' be summed over them where the declared parts would double-count their
#' overlaps.
#'
#' Internally, this is only used for evaluating measures, not during
#' optimisation or certification. Both require just a supremum, and a supremum
#' over a union is equivalently a maxmium of the suprema of its parts. Often the
#' declared cover is a simpler one to search over anyway; its parts are the ones
#' the caller stated. Disjoining produces cells that are smaller, more numerous
#' and cut along facets that may not be interesting in the problem setting.
#'
#' Parts meeting in a lower-dimensional slice are left overlapping, since
#' `setdiff()` computes closed differences and a slice has no measure to
#' double-count.
#'
#' @param x A [region].
#' @param ... Passed to `setdiff()`, e.g. `max_cells`.
#' @return A [region] covering the same set, whose parts have disjoint
#'   interiors, or `NULL` if `x` is empty.
#' @examples
#' # The two cells of the K = 3 plurality null overlap where candidate 1 trails
#' # both others. Peeling them apart leaves that region in one of the two.
#' plurality <- union(
#'   simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'   simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#' )
#' peeled <- disjoin(plurality)
#' n_parts(peeled)
#' setequal(peeled, plurality)
#' @seealso [region_algebra]
#' @export
disjoin <- new_generic("disjoin", "x", function(x, ...) S7::S7_dispatch())


#' @rdname disjoin
#' @usage NULL
method(disjoin, convex_region) <- function(x, ...) x


#' @rdname disjoin
#' @usage NULL
method(disjoin, union_region) <- function(x, ...) {
  kept <- list()
  for (part in x@parts) {
    remainder <- if (length(kept)) {
      without_slice_warning(setdiff(part, union_region(kept), ...))
    } else {
      part
    }
    # A part wholly covered by the ones before it contributes nothing, and an
    # empty one was never going to. Dropping them is the point: what comes back
    # is a cover with no redundant piece in it.
    if (!is.null(remainder) && !is_empty(remainder)) {
      kept <- c(kept, parts(remainder))
    }
  }
  if (length(kept) == 0L) {
    return(NULL)
  }
  union_region(kept)
}


#' One convex ambient part minus a list of parts, as exact H-matrices
#'
#' First-violated-facet decomposition: a point is outside `B` exactly when
#' some facet of `B` is violated, and taking the *first* violated facet makes
#' the pieces interior-disjoint:
#'
#'   B^c = union over facets f of
#'         { s_1 leq x_1, ..., s_(f-1) leq x_(f-1), s_f geq x_f }   (closures)
#'
#' Only facets that can actually be violated inside the ambient take part: a
#' a constraint that is implied by the ambient is dropped (by an exact LP),
#' which keeps the K-candidate plurality complement at one cell, for example.
#' The ambient's own rows are stacked into every piece, so cells never leave
#' it.
#'
#' Subtracting several parts multiplies: a cell of the difference picks one
#' piece per subtracted part (the Cartesian product), pruned via feasibility
#' `max_cells` bounds the product before it is expanded.
#' @keywords internal
#' @noRd
part_difference <- function(h_ambient, subtract, max_cells) {
  # A part that never meets the ambient subtracts nothing.
  hs <- Filter(
    \(h) !q_is_empty(q_rbind(h_ambient, h)),
    lapply(subtract, q_hrep)
  )
  pieces <- lapply(hs, \(h) complement_pieces(h_ambient, h))

  # NULL marks a part meeting the ambient only in a lower-dimensional slice,
  # which removes nothing there, so the part is skipped and the caller warns.
  sliced <- vapply(pieces, is.null, logical(1))
  pieces <- pieces[!sliced]
  n_sliced <- sum(sliced)

  # Every part missed the ambient (or only sliced it): the difference is the
  # ambient itself.
  if (length(pieces) == 0L) {
    return(list(cells = list(h_ambient), n_sliced = n_sliced))
  }
  # A part with no violable facet covers the ambient: nothing is left.
  n_cells <- prod(lengths(pieces))
  if (n_cells == 0L) {
    return(list(cells = list(), n_sliced = n_sliced))
  }
  if (n_cells > max_cells) {
    stop(
      "the difference would decompose into ",
      n_cells,
      " cells before pruning, above `max_cells = ",
      max_cells,
      "`. Subtract within a tighter region, or raise `max_cells`.",
      call. = FALSE
    )
  }

  combos <- expand.grid(lapply(pieces, seq_along))
  cells <- lapply(seq_len(nrow(combos)), function(i) {
    cell <- h_ambient
    for (j in seq_along(pieces)) {
      cell <- q_rbind(cell, pieces[[j]][[combos[i, j]]])
    }
    cell
  })
  list(cells = Filter(Negate(q_is_empty), cells), n_sliced = n_sliced)
}


#' The interior-disjoint pieces of one part's complement within an ambient
#' @keywords internal
#' @noRd
complement_pieces <- function(h_ambient, h_part) {
  implied <- function(row) {
    # The row stores `(l, b, -a)`; the test is `max { a . x : ambient } <= b`.
    peak <- q_maximum(h_ambient, q_neg(row[-(1:2)]))
    !is.null(peak) && q_leq(peak, row[2L])
  }
  eq <- h_part[, 1L] == "1"
  for (r in which(eq)) {
    row <- h_part[r, ]
    flipped <- q_reverse_ineq(q_subrows(h_part, r), 1L)
    if (!implied(row) || !implied(flipped[1L, ])) {
      # The part meets the ambient only in a lower-dimensional slice, whose
      # closed complement is the whole ambient: subtracting it removes
      # nothing. NULL tells the caller to skip the part and say so.
      return(NULL)
    }
  }

  ineq <- which(!eq)
  surviving <- ineq[!vapply(ineq, \(r) implied(h_part[r, ]), logical(1))]
  lapply(seq_along(surviving), function(j) {
    block <- q_subrows(h_part, surviving[seq_len(j)])
    q_reverse_ineq(block, j)
  })
}
