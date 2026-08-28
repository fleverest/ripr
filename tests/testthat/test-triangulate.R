# Tests for the vertex fan in R/triangulate.R.

# The `k`-volume of a simplex given by `k + 1` vertices, by Gram determinant
simplex_volume <- function(space) {
  edges <- space@vertices[, -1L, drop = FALSE] - space@vertices[, 1L]
  sqrt(max(det(crossprod(edges)), 0)) / factorial(ncol(edges))
}

# A regular `n`-gon on the unit circle, whose area is known in closed form.
polygon_region <- function(n) {
  angle <- 2 * pi * seq_len(n) / n
  polytope_region(vertices = rbind(cos(angle), sin(angle)))
}

polygon_area <- function(n) n * sin(2 * pi / n) / 2


# --- The fan ------------------------------------------------------------------

test_that("a pentagon fans into three triangles that tile it", {
  cells <- cells(polygon_region(5L))

  # Five edges, two of them through the apex: three cones, one per remaining
  # edge. Every one is a triangle, not merely a polytope.
  expect_length(cells, 3L)
  for (cell in cells) {
    expect_true(S7_inherits(cell, simplex_region))
    expect_identical(ncol(cell@vertices), 3L)
  }
  # Areas summing to the pentagon's is the tiling itself: an overlap would
  # overshoot and a gap would undershoot.
  expect_equal(
    sum(vapply(cells, simplex_volume, numeric(1))),
    polygon_area(5L)
  )
})


test_that("the fan tiles in three dimensions too", {
  cube <- polytope_region(
    vertices = t(as.matrix(expand.grid(c(0, 1), c(0, 1), c(0, 1))))
  )
  cells <- cells(cube)

  # The pulling triangulation of the cube is the minimal one: six tetrahedra.
  expect_length(cells, 6L)
  expect_equal(sum(vapply(cells, simplex_volume, numeric(1))), 1)
})


test_that("every cell of the fan is a cell of the region it came from", {
  hexagon <- polygon_region(6L)
  for (cell in cells(hexagon)) {
    for (i in seq_len(ncol(cell@vertices))) {
      expect_true(contains(hexagon, cell@vertices[, i]))
    }
    # The centroid is interior, so this catches a cell reflected or built from
    # vertices of the wrong facet, which the vertex test alone would not.
    expect_true(contains(hexagon, rowMeans(cell@vertices)))
  }
})


# --- Regions that are already simplices ---------------------------------------

test_that("a simplex is its own only cell, handed back untouched", {
  s <- simplex_region(vertices = diag(3))
  expect_identical(cells(s), list(s))

  # A point is the degenerate simplex and takes the same path, which matters
  # because `point_region` is a `polytope_region` and would otherwise be
  # triangulated at every call for no gain.
  p <- point_region(theta = c(0.5, 0.3, 0.2))
  expect_identical(cells(p), list(p))
})


test_that("a polytope that happens to be a simplex still fans to one cell", {
  # Declared as a polytope, so it takes the fan rather than the identity, and
  # the fan's recursion floor is what stops it at one cell.
  triangle <- polytope_region(vertices = cbind(c(0, 0), c(1, 0), c(0, 1)))
  cells <- cells(triangle)
  expect_length(cells, 1L)
  expect_true(S7_inherits(cells[[1L]], simplex_region))
  expect_equal(simplex_volume(cells[[1L]]), 0.5)
})


# --- Degenerate inputs --------------------------------------------------------

test_that("a lower-dimensional polytope is triangulated within its own hull", {
  # A square lying in the plane `theta_3 = 0`. Its cells are 2-simplices in
  # `R^3`, so the recursion has to read the hull's dimension off the equality
  # row rather than assume the ambient one.
  square <- polytope_region(
    vertices = cbind(c(0, 0, 0), c(1, 0, 0), c(1, 1, 0), c(0, 1, 0))
  )
  cells <- cells(square)
  expect_length(cells, 2L)
  for (cell in cells) {
    expect_true(S7_inherits(cell, simplex_region))
    expect_identical(ncol(cell@vertices), 3L)
  }
})


test_that("a segment in R^3 is one cell, and a point is one cell", {
  segment <- polytope_region(vertices = cbind(c(0, 0, 0), c(1, 1, 1)))
  cells <- cells(segment)
  expect_length(cells, 1L)
  expect_true(S7_inherits(cells[[1L]], simplex_region))
  expect_identical(ncol(cells[[1L]]@vertices), 2L)

  single <- polytope_region(vertices = matrix(c(2, 3, 4), ncol = 1L))
  expect_length(cells(single), 1L)
})


test_that("collinear vertices drop to the two that are extreme", {
  # The midpoint is not a vertex of the hull, and a cell carrying it would be
  # three affinely dependent points, which `simplex_region()` refuses. The fan
  # drops it before it can reach a cell.
  collinear <- polytope_region(
    vertices = cbind(c(0, 0), c(1, 1), c(2, 2))
  )
  cells <- cells(collinear)
  expect_length(cells, 1L)
  expect_equal(
    cells[[1L]]@vertices[, order(cells[[1L]]@vertices[1L, ])],
    cbind(c(0, 0), c(2, 2))
  )
})


test_that("an unbounded region is refused, naming what it has", {
  expect_error(
    triangulate(halfspace_region(normal = c(1, -1, 0))),
    "rays or lineality"
  )
  expect_error(
    triangulate(unconstrained_region(2L)),
    "only a bounded region can be triangulated"
  )
})


# --- Exactness ----------------------------------------------------------------

test_that("triangulating an algebra cell starts from its exact form", {
  # The complement of the K = 3 plurality null within the simplex is a
  # quadrilateral with the barycentre as one vertex. Its exact coordinate is
  # 1/3, which is not a double, so a fan that worked from the region's rounded
  # vertices would carry the rounding into the cells and into any algebra done
  # on them afterwards.
  plurality_cell <- function(j) {
    v <- diag(3)
    v[, 1L] <- replace(numeric(3), c(1L, j), 0.5)
    simplex_region(vertices = v)
  }
  quad <- setdiff(
    simplex_region(vertices = diag(3)),
    union(plurality_cell(2L), plurality_cell(3L))
  )
  cells <- cells(quad)
  expect_length(cells, 2L)

  # The barycentre survives as the rational 1/3 rather than as `d2q(1/3)`,
  # which is the 53-bit fraction the double rounds to.
  bary <- vapply(
    cells,
    function(cell) {
      any(apply(cell@hv@qv[, -(1:2), drop = FALSE], 1L, \(r) {
        all(r == "1/3")
      }))
    },
    logical(1)
  )
  expect_true(any(bary))

  # And the cells still tile the quadrilateral, whose area is the simplex's
  # less the two plurality cells' shares of it.
  expect_equal(
    sum(vapply(cells, simplex_volume, numeric(1))),
    simplex_volume(simplex_region(vertices = diag(3))) / 3
  )
})


test_that("cells keep facets exact rather than re-deriving them from doubles", {
  cell <- cells(polygon_region(5L))[[1L]]
  expect_false(is.null(cell@hv@qh))
  expect_false(is.null(cell@facets))
  # The record is the exact pair for the same set the doubles approximate.
  expect_equal(from_hmatrix(cell@hv@qh)$a, cell@facets$a)
  expect_equal(from_hmatrix(cell@hv@qh)$b, cell@facets$b)
})


test_that("triangulation refuses to fan past max_cells", {
  # A pentagon fans into three triangles; a cap of two must refuse rather
  # than return a truncated tiling.
  pent <- polygon_region(5L)
  expect_length(cells(pent), 3L)
  expect_error(cells(pent, max_cells = 2L), "gave up after `max_cells = 2`")
})


test_that("max_cells caps a union's triangulation in total, not per part", {
  pent <- polygon_region(5L) # three triangles each
  u <- union_region(pent, pent)
  expect_length(cells(u, max_cells = 6L), 6L)
  expect_error(cells(u, max_cells = 5L), "gave up after `max_cells = 5`")
})


test_that("a null names the part that failed to decompose", {
  local_mocked_bindings(triangulate = function(space, ...) {
    stop("boom", call. = FALSE)
  })
  fam <- multinomial_family(n_trials = 2L, k = 3L)
  square <- polytope_region(
    vertices = cbind(
      c(0.5, 0.5, 0),
      c(0, 0.5, 0.5),
      c(0, 0, 1),
      c(0.5, 0, 0.5)
    )
  )
  expect_error(
    null_model(fam, list(simplex_region(vertices = diag(3)), square)),
    "part 2 of the null: boom"
  )
})
