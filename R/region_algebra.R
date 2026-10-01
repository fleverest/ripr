#' @include region.R
NULL

# --- Set algebra on regions ---------------------------------------------------

#' Build one region from an exact rational H-representation
#'
#' The one place algebra results return to doubles, rounding once at the end of
#' the chain. Picks the most specific class the cell qualifies for.
#' @keywords internal
#' @noRd
region_from_qh <- function(qh) {
  hv <- hv_fill(qh = qh)
  if (nrow(hv@v$r) > 0L || nrow(hv@v$l) > 0L) {
    polyhedron_region(.hv = hv)
  } else if (is.null(simplex_defect(hv))) {
    simplex_region(.hv = hv)
  } else {
    polytope_region(.hv = hv)
  }
}


# --- Empty region -------------------------------------------------------------

#' The region with nothing in it
#'
#' What the set algebra returns when nothing remains, e.g. `x & y` for disjoint
#' regions or `x - y` when `y` covers `x`. It has no [parts()] or [cells()],
#' contains no point, and is bounded. [null_model()] refuses one.
#'
#' @return An `empty_region`.
#' @examples
#' # Disjoint regions intersect in nothing:
#' nothing <- point_region(theta = c(1, 0, 0)) &
#'   point_region(theta = c(0, 1, 0))
#' is_empty(nothing)
#' length(parts(nothing))
#'
#' # And the algebra keeps going from there:
#' identical(
#'   nothing | simplex_region(vertices = diag(3)),
#'   simplex_region(vertices = diag(3))
#' )
#' @seealso [region_algebra]
#' @export
empty_region <- new_class(
  "empty_region",
  parent = region,
  constructor = function() {
    new_object(S7_object())
  }
)


method(space_dim, empty_region) <- function(space) {
  stop(
    "`space_dim()` is not defined for an `empty_region`: the empty set is ",
    "the same set in every ambient dimension.",
    call. = FALSE
  )
}


method(parts, empty_region) <- function(space) list()


#' @description An empty region has no cells.
#' @rdname cells
#' @usage NULL
method(cells, empty_region) <- function(
  space,
  max_cells = getOption("ripr.max_cells", 1000L),
  .budget = NULL
) {
  list()
}


method(contains, empty_region) <- function(space, theta, tol = 1e-8) FALSE


method(is_empty, empty_region) <- function(space) TRUE


method(is_bounded, empty_region) <- function(space) TRUE


method(region_phrase, empty_region) <- function(space) {
  "the empty region"
}


#' @rdname empty_region
#' @usage NULL
method(format, empty_region) <- function(x, ...) {
  sprintf("%s: %s", class_name(x), region_phrase(x))
}


#' @rdname empty_region
#' @usage NULL
method(print, empty_region) <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}


# --- Union region -------------------------------------------------------------

#' A finite union of convex regions
#'
#' The union \eqn{\bigcup_i \Theta_{0i}}{union_i Theta_0i} of finitely many
#' [convex_region]s, e.g. a null hypothesis or a truncated prior's support. It
#' is not a [convex_region], so has no [chart()] or [project()]; algorithms
#' run over its [parts()] or [cells()] instead. Given exactly one convex
#' region, `union_region()` returns it unchanged.
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
#'   simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'   simplex_region(vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#' )
#'
#' # Nesting is flattened, so these agree:
#' s <- simplex_region(vertices = diag(3))
#' h <- halfspace_region(normal = c(1, -1, 0))
#' length(parts(union_region(s, h)))
#' length(parts(union_region(list(s, h))))
#' length(parts(union_region(union_region(s), list(h))))
#'
#' # One cell is already a region, so it is handed back as it came:
#' identical(union_region(s), s)
#' @seealso [region_algebra] for the set-operation operators: `x | y`
#'   forwards to this constructor, while `x & y` and `x - y` compute new
#'   regions from old.
#' @export
union_region <- new_class(
  "union_region",
  parent = region,
  properties = list(parts = class_list),
  constructor = function(...) {
    flat <- flatten_parts(list(...))
    empties <- vapply(flat, \(p) S7_inherits(p, empty_region), logical(1))
    if (any(empties)) {
      flat <- flat[!empties]
      if (length(flat) == 0L) {
        return(empty_region())
      }
    }
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
    # Ambient dimension only; codimension is free.
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
#' Non-regions are left as leaves for the validator to name.
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


#' Coerce a [region], or a list of them, to a [region]
#' @keywords internal
#' @noRd
as_region <- function(x) {
  if (S7_inherits(x, region)) x else union_region(x)
}


method(space_dim, union_region) <- function(space) {
  space_dim(space@parts[[1L]])
}


method(contains, union_region) <- function(space, theta, tol = 1e-8) {
  any(vapply(space@parts, \(p) contains(p, theta, tol), logical(1)))
}


method(is_empty, union_region) <- function(space) {
  all(vapply(space@parts, is_empty, logical(1)))
}


method(is_bounded, union_region) <- function(space) {
  all(vapply(space@parts, is_bounded, logical(1)))
}


method(parts, union_region) <- function(space) space@parts


#' @description A union's cells are its parts' cells, flattened, without
#'   exceeding the shared `max_cells` budget.
#' @rdname cells
#' @usage NULL
method(cells, union_region) <- function(
  space,
  max_cells = getOption("ripr.max_cells", 1000L),
  .budget = NULL
) {
  budget <- if (is.null(.budget)) cell_budget(max_cells) else .budget
  unlist(
    lapply(space@parts, \(p) cells(p, .budget = budget)),
    recursive = FALSE
  )
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
  cat("<", class_name(x), ">\n", sep = "")
  cat("  ", parts_label(n), ", dimension ", space_dim(x), "\n", sep = "")
  cat_parts(x@parts)
  invisible(x)
}


#' List parts under a print banner
#'
#' One line per part, or a tally by class beyond six.
#' @keywords internal
#' @noRd
cat_parts <- function(parts) {
  if (length(parts) <= 6L) {
    for (p in parts) {
      cat("    ", format(p), "\n", sep = "")
    }
    return(invisible())
  }
  tally <- table(vapply(parts, \(p) class_name(p), character(1)))
  for (nm in names(tally)) {
    cat("    ", tally[[nm]], " x ", nm, "\n", sep = "")
  }
}


#' @description `format()` gives the same summary on one line.
#' @rdname union_region
#' @usage NULL
#' @export
method(format, union_region) <- function(x, ...) {
  sprintf(
    "%s: %s, dimension %d",
    class_name(x),
    parts_label(length(x@parts)),
    space_dim(x)
  )
}


#' Set algebra on regions
#'
#' Regions combine with `|` (union), `&` (intersection), `-` (difference) and
#' `==` (set equality), decided in exact rational arithmetic. Base R's
#' [union()], [intersect()], [setdiff()] and [setequal()] do not accept regions.
#'
#' `x | y` is structural: it forwards to [union_region()]. `x & y` intersects
#' every part of `x` with every part of `y` and drops the empty ones; parts of
#' the result may overlap.
#'
#' `x - y` is the closed difference. A part of `y` that meets `x` only in a
#' lower-dimensional slice subtracts nothing and raises a warning. To remove
#' several regions, subtract their union: `x - (y1 | y2)`. `-` errors past
#' `getOption("ripr.max_cells", 1000L)` cells; the same option sets the default
#' `max_cells` for [cells()] and [null_model()].
#'
#' `x == y` returns a single `TRUE` or `FALSE`, testing containment both ways
#' (exactly, from generators, when the containing side is convex; by a
#' difference otherwise). It compares sets, not objects (use [identical()] for
#' that); regions of different dimensions are unequal.
#'
#' @param e1,e2 [region]s.
#' @return `==` and `!=` return `TRUE` or `FALSE`. The other operators return a
#'   [region]: `x | y` returns the same as [union_region()]; `x & y` returns a
#'   [union_region] of the surviving cells, a single cell alone when it is the
#'   only survivor, or an [empty_region()] if the intersection is empty; and
#'   `x - y` likewise returns an [empty_region()] when nothing remains.
#' @examples
#' # The two cells of the K = 3 plurality null: the union of the regions where
#' # candidate 2 beats candidate 1, and where candidate 3 beats candidate 1.
#' loses_2 <- simplex_region(
#'   vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))
#' )
#' loses_3 <- simplex_region(
#'   vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1))
#' )
#'
#' # The null itself, by operator rather than constructor:
#' plurality <- loses_2 | loses_3
#' plurality
#'
#' # Its two cells meet where candidate 1 trails both others:
#' loses_2 & loses_3
#'
#' # Disjoint regions intersect in nothing:
#' point_region(theta = c(1, 0, 0)) & point_region(theta = c(0, 1, 0))
#'
#' # The complement of the null within the simplex is the region where
#' # candidate 1 wins -- the alternative, as a region:
#' simplex <- simplex_region(vertices = diag(3))
#' simplex - plurality
#'
#' # Several regions are subtracted as their union, so this is the same set:
#' (simplex - (trails_2 | trails_3)) == (simplex - plurality)
#'
#' # A square, and the same square cut into two triangles: different objects,
#' # the same set. One of those triangles alone is not.
#' square <- polytope_region(vertices = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)))
#' square == union_region(cells(square))
#' square == cells(square)[[1]]
#' square != cells(square)[[1]]
#'
#' # The cap on a difference's decomposition is an option:
#' old <- options(ripr.max_cells = 5000L)
#' simplex - plurality
#' options(old)
#' @seealso [union_region()], [empty_region()], and [disjoin()], which turns a
#'   union into one whose parts do not overlap.
#' @name region_algebra
NULL


#' The default cap on a decomposition's cell count, for the operators, which
#' have no argument to put it in
#' @keywords internal
#' @noRd
max_cells_option <- function() getOption("ripr.max_cells", 1000L)


#' @rdname region_algebra
#' @usage NULL
method(`|`, list(region, region)) <- function(e1, e2) union_region(e1, e2)


#' @rdname region_algebra
#' @usage NULL
method(`&`, list(region, region)) <- function(e1, e2) {
  region_intersect(e1, e2)
}


#' @rdname region_algebra
#' @usage NULL
method(`-`, list(region, region)) <- function(e1, e2) {
  region_difference(e1, e2, max_cells_option())
}


#' @rdname region_algebra
#' @usage NULL
method(`==`, list(region, region)) <- function(e1, e2) {
  region_equal(e1, e2, max_cells_option())
}


#' @rdname region_algebra
#' @usage NULL
method(`!=`, list(region, region)) <- function(e1, e2) {
  !region_equal(e1, e2, max_cells_option())
}


#' Refuse two regions of different ambient dimensions, naming both
#' @keywords internal
#' @noRd
check_same_dim <- function(x, y) {
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
  invisible(NULL)
}


#' The intersection of two regions, behind `&`
#' @keywords internal
#' @noRd
region_intersect <- function(x, y) {
  if (S7_inherits(x, empty_region) || S7_inherits(y, empty_region)) {
    return(empty_region())
  }
  check_same_dim(x, y)
  hs <- lapply(parts(y), q_hrep)
  acc <- unlist(
    lapply(parts(x), \(p) {
      hp <- q_hrep(p)
      lapply(hs, \(h) q_rbind(hp, h))
    }),
    recursive = FALSE
  )
  acc <- Filter(Negate(q_is_empty), acc)
  if (length(acc) == 0L) {
    return(empty_region())
  }
  union_region(lapply(acc, region_from_qh))
}


#' The closed difference of two regions, behind `-`
#' @keywords internal
#' @noRd
region_difference <- function(
  x,
  y,
  max_cells = getOption("ripr.max_cells", 1000L)
) {
  if (S7_inherits(x, empty_region) || S7_inherits(y, empty_region)) {
    return(x)
  }
  check_same_dim(x, y)
  # (A1 u A2) \ y = (A1 \ y) u (A2 \ y).
  results <- lapply(
    parts(x),
    \(ambient) part_difference(q_hrep(ambient), parts(y), max_cells)
  )
  n_sliced <- sum(vapply(results, \(r) r$n_sliced, integer(1)))
  if (n_sliced > 0L) {
    warning(slice_warning(paste0(
      "in `x - y`, in ",
      count_label(n_sliced, "case"),
      ", a part of `y` met a part of `x` only in a lower-dimensional slice; ",
      "nothing was subtracted there, since a closed difference removes ",
      "nothing from a slice."
    )))
  }
  cells <- unlist(lapply(results, \(r) r$cells), recursive = FALSE)
  if (length(cells) == 0L) {
    return(empty_region())
  }
  union_region(lapply(cells, region_from_qh))
}


#' The warning `-` raises when a part subtracts nothing
#'
#' Classed warning so that `disjoin()` can silence it.
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
    ripr_slice_warning = function(w) invokeRestart("muffleWarning")
  )
}


#' Whether two regions are the same set, behind `==`
#' @keywords internal
#' @noRd
region_equal <- function(
  x,
  y,
  max_cells = getOption("ripr.max_cells", 1000L)
) {
  if (S7_inherits(x, empty_region) || S7_inherits(y, empty_region)) {
    return(is_empty(x) && is_empty(y))
  }
  if (space_dim(x) != space_dim(y)) {
    # Unlike `&` and `-`, not an error: the sets are simply different.
    return(FALSE)
  }
  region_subset(x, y, max_cells) && region_subset(y, x, max_cells)
}


#' Is every point of `inner` in `whole`?
#' @keywords internal
#' @noRd
region_subset <- function(
  inner,
  whole,
  max_cells = getOption("ripr.max_cells", 1000L)
) {
  if (S7_inherits(whole, convex_region)) {
    qh <- q_hrep(whole)
    return(all(vapply(
      parts(inner),
      \(p) q_holds(q_scdd(q_hrep(p)), qh),
      logical(1)
    )))
  }
  subtract <- parts(whole)
  for (p in parts(inner)) {
    if (length(part_difference(q_hrep(p), subtract, max_cells)$cells) > 0L) {
      return(FALSE)
    }
  }
  TRUE
}


# --- Disjoining ---------------------------------------------------------------

#' Transform a region's parts into a disjoint cover of the union
#'
#' Subtracts from each part all the parts before it, so the result covers the
#' same set with parts meeting only on boundaries, and a measure can be summed
#' over them without double-counting. Parts meeting only in a lower-dimensional
#' slice are left overlapping (see [region_algebra]).
#'
#' @param x A [region].
#' @param ... For a [union_region], `max_cells`: the cap on the number of
#'   cells each difference may decompose into before giving up. Defaults to
#'   the `ripr.max_cells` option, or `1000L` when that is unset.
#' @return A [region] covering the same set, whose parts have disjoint
#'   interiors, or an [empty_region()] if `x` is empty.
#' @examples
#' # The two cells of the K = 3 plurality null overlap where candidate 1 trails
#' # both others. Peeling them apart leaves that region in one of the two.
#' plurality <-
#'   simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))) |
#'   simplex_region(vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#' peeled <- disjoin(plurality)
#' length(parts(peeled))
#' peeled == plurality
#' @seealso [region_algebra]
#' @export
disjoin <- new_generic("disjoin", "x", function(x, ...) S7::S7_dispatch())


#' @rdname disjoin
#' @usage NULL
method(disjoin, convex_region) <- function(x, ...) x


#' @rdname disjoin
#' @usage NULL
method(disjoin, empty_region) <- function(x, ...) x


#' @rdname disjoin
#' @usage NULL
method(disjoin, union_region) <- function(
  x,
  ...,
  max_cells = getOption("ripr.max_cells", 1000L)
) {
  kept <- list()
  for (part in x@parts) {
    remainder <- if (length(kept)) {
      without_slice_warning(
        region_difference(part, union_region(kept), max_cells)
      )
    } else {
      part
    }
    if (!is_empty(remainder)) {
      kept <- c(kept, parts(remainder))
    }
  }
  if (length(kept) == 0L) {
    return(empty_region())
  }
  union_region(kept)
}


#' One convex ambient part minus a list of parts, as exact H-matrices
#'
#' First-violated-facet decomposition, interior-disjoint:
#'
#'   B^c = union over facets f of
#'         { s_1 leq x_1, ..., s_(f-1) leq x_(f-1), s_f geq x_f }   (closures)
#'
#' Facets implied by the ambient are dropped. Parts are subtracted one at a
#' time from every cell so far, keeping only full-dimensional cells, and at
#' most `max_cells` of them. Returns `list(cells, n_sliced)`.
#' @keywords internal
#' @noRd
part_difference <- function(h_ambient, subtract, max_cells) {
  hs <- Filter(
    \(h) !q_is_empty(q_rbind(h_ambient, h)),
    lapply(subtract, q_hrep)
  )
  # A part meeting the ambient only in a lower-dimensional slice removes
  # nothing there, so it is skipped and the caller warns.
  sliced <- vapply(
    hs,
    \(h) is.null(complement_pieces(h_ambient, h)),
    logical(1)
  )
  dim_ambient <- q_dim(h_ambient)
  full <- \(m) !q_is_empty(m) && q_dim(m) == dim_ambient

  cells <- list(h_ambient)
  for (h in hs[!sliced]) {
    cells <- unlist(
      lapply(cells, function(cell) {
        if (!full(q_rbind(cell, h))) {
          return(list(cell))
        }
        pieces <- lapply(complement_pieces(cell, h), \(p) q_rbind(cell, p))
        Filter(full, pieces)
      }),
      recursive = FALSE
    )
    if (length(cells) > max_cells) {
      stop(
        "the difference decomposes into more than `max_cells = ",
        max_cells,
        "` cells. Subtract within a tighter region, or raise `max_cells`.",
        call. = FALSE
      )
    }
  }
  list(cells = cells, n_sliced = sum(sliced))
}


#' The interior-disjoint pieces of one part's complement within an ambient
#' @keywords internal
#' @noRd
complement_pieces <- function(h_ambient, h_part) {
  qv <- q_scdd(h_ambient)
  implied <- function(r) q_holds(qv, h_part[r, , drop = FALSE])
  eq <- h_part[, 1L] == "1"
  for (r in which(eq)) {
    if (!implied(r)) {
      # Meets the ambient only in a slice: removes nothing (see caller).
      return(NULL)
    }
  }

  ineq <- which(!eq)
  surviving <- ineq[!vapply(ineq, implied, logical(1))]
  lapply(seq_along(surviving), function(j) {
    block <- q_subrows(h_part, surviving[seq_len(j)])
    q_reverse_ineq(block, j)
  })
}
