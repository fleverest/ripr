# Tests for the set algebra in R/region_algebra.R.
#
# The algebra composes on the rational layer and converts once, which is what
# the exactness assertions here pin: derived coordinates equal their rational
# values to the last bit, not to a tolerance.

# One cell of the K-candidate plurality null: candidate 1 trails candidate j.
plurality_cell <- function(k, j) {
  vertices <- diag(k)
  vertices[, 1L] <- replace(numeric(k), c(1L, j), 0.5)
  simplex_region(vertices = vertices)
}

# Columns sorted lexicographically, so vertex sets compare independently of
# the order cddlib emits them in.
sorted_columns <- function(m) {
  m[, do.call(order, asplit(m, 1L)), drop = FALSE]
}

# --- union / intersect dispatch -----------------------------------------------

test_that("union() and intersect() fall through to base R off regions", {
  # Just a sanity check for me because I live in fear of breaking stuff
  expect_identical(union(c(1, 2), c(2, 3)), c(1, 2, 3))
  expect_identical(intersect(c(1, 2), c(2, 3)), 2)
  expect_identical(union(letters[1:2], letters[2:3]), letters[1:3])
})


test_that("union() on regions is the union_region constructor", {
  s <- plurality_cell(3L, 2L)
  h <- halfspace_region(normal = c(1, -1, 0))
  expect_identical(union(s, h), union_region(s, h))
  # A union argument flattens rather than nests, exactly as the constructor.
  expect_identical(n_parts(union(union(s, h), plurality_cell(3L, 3L))), 3L)
})


# --- intersect ----------------------------------------------------------------

test_that("two overlapping plurality cells intersect in the exact triangle", {
  ab <- intersect(plurality_cell(3L, 2L), plurality_cell(3L, 3L))

  # {theta1 <= theta2} meets {theta1 <= theta3} where candidate 1 trails both:
  # the triangle spanned by the two loser vertices and the barycentre. Its
  # coordinates are exactly 0, 1 and double(1/3) -- the acceptance test for
  # composing in rationals and converting once.
  expect_true(S7_inherits(ab, polytope_region))
  expect_identical(
    sorted_columns(ab@vertices),
    sorted_columns(cbind(c(0, 1, 0), c(0, 0, 1), c(1, 1, 1) / 3))
  )
})


test_that("intersect() distributes over the parts of unions", {
  u <- union_region(plurality_cell(3L, 2L), plurality_cell(3L, 3L))
  h <- halfspace_region(normal = c(0, 1, -1)) # theta2 <= theta3

  r <- intersect(u, h)
  expect_true(S7_inherits(r, region))

  # Pointwise agreement with the definition of intersection.
  set.seed(11)
  for (i in seq_len(200L)) {
    theta <- as.numeric(stats::rexp(3L))
    theta <- theta / sum(theta)
    expect_identical(
      contains(r, theta),
      contains(u, theta) && contains(h, theta)
    )
  }
})


test_that("a disjoint intersection is empty, not a degenerate cell", {
  expect_null(intersect(
    point_region(theta = c(1, 0, 0)),
    point_region(theta = c(0, 1, 0))
  ))

  # Cells meeting only in a shared face are not empty: the closed cells of a
  # cover genuinely intersect in that face.
  edge <- intersect(
    simplex_region(vertices = cbind(c(0, 0), c(1, 0), c(0, 1))),
    simplex_region(vertices = cbind(c(1, 1), c(1, 0), c(0, 1)))
  )
  expect_false(is.null(edge))
  expect_true(contains(edge, c(0.5, 0.5)))
  expect_false(contains(edge, c(0.25, 0.25)))
})


test_that("an unbounded intersection returns a polyhedron_region", {
  quadrant <- intersect(
    halfspace_region(normal = c(1, 0)),
    halfspace_region(normal = c(0, 1))
  )
  expect_true(S7_inherits(quadrant, polyhedron_region))
  expect_false(is_bounded(quadrant))
  expect_true(contains(quadrant, c(-3, -5)))
  expect_false(contains(quadrant, c(1, 1)))
})


test_that("intersect() takes more than two regions and refuses mismatched dimensions", {
  # Three halfspaces cutting the plane down to a bounded triangle.
  tri <- intersect(
    halfspace_region(normal = c(-1, 0)), # x >= 0
    halfspace_region(normal = c(0, -1)), # y >= 0
    halfspace_region(normal = c(1, 1), offset = 1) # x + y <= 1
  )
  expect_true(is_bounded(tri))
  expect_identical(
    sorted_columns(v_rep(tri)$v),
    sorted_columns(cbind(c(0, 0), c(1, 0), c(0, 1)))
  )

  expect_error(
    intersect(unconstrained_region(2L), unconstrained_region(3L)),
    "ambient dimension"
  )
})


# --- setdiff ------------------------------------------------------------------

test_that("setdiff() falls through to base R off regions", {
  expect_identical(setdiff(c(1, 2), c(2, 3)), 1)
})


test_that("subtracting one plurality cell from the simplex gives one part", {
  # Every facet of the cell except its own boundary is a wall of the simplex,
  # dropped by the exact LP; what remains is the single cell where candidate 1
  # beats candidate j.
  for (k in c(3L, 4L, 5L)) {
    ambient <- simplex_region(vertices = diag(k))
    left <- setdiff(ambient, plurality_cell(k, 2L))
    expect_true(S7_inherits(left, convex_region))
    expect_identical(n_parts(left), 1L)
  }
})


test_that("the complement of the plurality null is the candidate-1-wins region", {
  ambient <- simplex_region(vertices = diag(3))
  null_region <- union_region(plurality_cell(3L, 2L), plurality_cell(3L, 3L))
  wins <- setdiff(ambient, null_region)

  expect_true(contains(wins, c(0.5, 0.3, 0.2)))
  expect_false(contains(wins, c(0.2, 0.5, 0.3), tol = 1e-12))

  # Every cell of a difference inside the simplex still lives on the simplex:
  # exactly one equality row, and every vertex sums to one.
  for (cell in parts(wins)) {
    h <- h_rep(cell)
    expect_identical(sum(h$eq), 1L)
    expect_true(max(abs(colSums(v_rep(cell)$v) - 1)) <= 1e-12)
  }
})


test_that("K = 5 complement of the whole plurality null is fast and small", {
  ambient <- simplex_region(vertices = diag(5L))
  null_region <- union_region(lapply(2:5, \(j) plurality_cell(5L, j)))
  elapsed <- system.time(wins <- setdiff(ambient, null_region))[["elapsed"]]
  expect_lt(elapsed, 1)
  expect_identical(n_parts(wins), 1L)
  expect_true(contains(wins, c(0.6, 0.1, 0.1, 0.1, 0.1)))
})


test_that("the double difference agrees with the original", {
  ambient <- simplex_region(vertices = diag(3))
  null_region <- union_region(plurality_cell(3L, 2L), plurality_cell(3L, 3L))
  back <- setdiff(ambient, setdiff(ambient, null_region))

  set.seed(13)
  agree <- vapply(
    seq_len(10000L),
    function(i) {
      theta <- as.numeric(stats::rexp(3L))
      theta <- theta / sum(theta)
      contains(back, theta) == contains(null_region, theta)
    },
    logical(1)
  )
  expect_gte(mean(agree), 0.999)
})


test_that("subtracting a lower-dimensional slice warns and removes nothing", {
  expect_warning(
    back <- setdiff(
      unconstrained_region(3L),
      simplex_region(vertices = diag(3))
    ),
    "lower-dimensional"
  )
  expect_identical(n_parts(back), 1L)
  expect_true(contains(back, c(5, -3, 2)))

  # A shared affine hull is fine: both live on the simplex.
  expect_no_warning(
    setdiff(simplex_region(vertices = diag(3)), plurality_cell(3L, 2L))
  )
})


test_that("a coarser ambient gives strictly more cells", {
  # The triangle shares two edges with the small square but none with the big
  # one, so fewer of its facets are dropped as ambient-implied.
  triangle <- simplex_region(vertices = cbind(c(0, 0), c(1, 0), c(0, 1)))
  small <- polytope_region(vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)))
  big <- polytope_region(
    vertices = 4 * cbind(c(-1, -1), c(1, -1), c(1, 1), c(-1, 1))
  )

  n_small <- n_parts(setdiff(small, triangle))
  n_big <- n_parts(setdiff(big, triangle))
  expect_gt(n_big, n_small)
})


test_that("x inside y leaves nothing, and max_cells errors rather than hangs", {
  expect_null(setdiff(
    plurality_cell(3L, 2L),
    simplex_region(vertices = diag(3))
  ))

  big <- polytope_region(
    vertices = 4 * cbind(c(-1, -1), c(1, -1), c(1, 1), c(-1, 1))
  )
  squares <- union_region(lapply(0:2, function(k) {
    polytope_region(
      vertices = cbind(c(k, 0), c(k + 0.5, 0), c(k + 0.5, 0.5), c(k, 0.5))
    )
  }))
  expect_error(setdiff(big, squares, max_cells = 10L), "max_cells")
})


test_that("difference cells are interior-disjoint by construction", {
  big <- polytope_region(
    vertices = 2 * cbind(c(-1, -1), c(1, -1), c(1, 1), c(-1, 1))
  )
  triangle <- simplex_region(vertices = cbind(c(0, 0), c(1, 0), c(0, 1)))
  left <- setdiff(big, triangle)

  # Interior points of one cell belong to no other cell.
  set.seed(17)
  for (cell in parts(left)) {
    ch <- chart(cell)
    for (i in seq_len(20L)) {
      theta <- ch$to_theta(stats::rnorm(ch$n_par))
      others <- Filter(\(p) !identical(p, cell), parts(left))
      strictly_inside <- all(vapply(
        seq_len(nrow(h_rep(cell)$a)),
        \(r) {
          h <- h_rep(cell)
          h$eq[r] || sum(h$a[r, ] * theta) < h$b[r] - 1e-9
        },
        logical(1)
      ))
      if (strictly_inside) {
        expect_false(any(vapply(
          others,
          \(p) contains(p, theta, tol = 1e-9),
          logical(1)
        )))
      }
    }
  }
})


test_that("algebra cells keep their exact rational representation", {
  # `q_hrep()` of a produced cell must serve the exact rows, not re-derive
  # from the rounded vertices: a quadrilateral on the simplex with a 1/3
  # vertex rounds to four points that exact arithmetic sees as an ulp-thin
  # tetrahedron -- full-dimensional, no equality row, sliver facets.
  ambient <- simplex_region(vertices = diag(3))
  null_region <- union_region(plurality_cell(3L, 2L), plurality_cell(3L, 3L))
  wins <- setdiff(ambient, null_region)

  h <- from_hmatrix(q_hrep(wins))
  expect_identical(sum(h$eq), 1L)
  expect_false(is.null(wins@q_cache$h))

  # A declared region carries the exact rationals of its own doubles, kept
  # from the derivation its constructor ran anyway.
  expect_false(is.null(ambient@q_cache$h))
  expect_identical(sum(from_hmatrix(ambient@q_cache$h)$eq), 1L)
})


test_that("algebra cells that are simplices come back as simplex_region", {
  # The K = 3 single-cell complement has three vertices on the simplex: a
  # certifiable simplex, and classified as one.
  left <- setdiff(simplex_region(vertices = diag(3)), plurality_cell(3L, 2L))
  expect_true(S7_inherits(left, simplex_region))

  # The intersection triangle of the two plurality cells likewise.
  ab <- intersect(plurality_cell(3L, 2L), plurality_cell(3L, 3L))
  expect_true(S7_inherits(ab, simplex_region))

  # The candidate-1-wins region is a quadrilateral: a polytope, not a simplex.
  wins <- setdiff(
    simplex_region(vertices = diag(3)),
    union_region(plurality_cell(3L, 2L), plurality_cell(3L, 3L))
  )
  expect_false(S7_inherits(wins, simplex_region))
  expect_true(S7_inherits(wins, polytope_region))
})


test_that("intersection reduces dimension; only setdiff refuses slices", {
  # The positive cone meets the sum-one hyperplane in the probability
  # simplex, exactly.
  cone <- polyhedron_region(rays = diag(3))
  plane <- h_region(a = matrix(1, 1L, 3L), b = 1, eq = TRUE)
  simplex <- intersect(cone, plane)
  expect_true(S7_inherits(simplex, simplex_region))
  expect_identical(
    sorted_columns(simplex@vertices),
    sorted_columns(diag(3) + 0)
  )
})


test_that("a subtracted part that never meets x subtracts nothing", {
  square <- polytope_region(
    vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  # A segment on a hyperplane elsewhere: lower-dimensional, but disjoint from
  # the square, so the difference is the square itself rather than a refusal.
  far_segment <- simplex_region(vertices = cbind(c(5, 0), c(5, 1)))
  left <- setdiff(square, far_segment)
  expect_identical(n_parts(left), 1L)
  expect_true(contains(left, c(0.5, 0.5)))

  # A slice actually through the square warns and leaves the square whole.
  through <- simplex_region(vertices = cbind(c(0.5, 0), c(0.5, 1)))
  expect_warning(whole <- setdiff(square, through), "lower-dimensional")
  expect_identical(n_parts(whole), 1L)
  expect_true(contains(whole, c(0.5, 0.5)))

  # And a disjoint full-dimensional part no longer multiplies cells.
  far_square <- polytope_region(
    vertices = cbind(c(5, 5), c(6, 5), c(6, 6), c(5, 6))
  )
  triangle <- simplex_region(vertices = cbind(c(0, 0), c(1, 0), c(0, 1)))
  expect_identical(
    n_parts(setdiff(square, union_region(triangle, far_square))),
    n_parts(setdiff(square, triangle))
  )
})


# --- setequal -----------------------------------------------------------------

test_that("setequal() falls through to base R off regions", {
  expect_true(setequal(c(1, 2), c(2, 1)))
  expect_false(setequal(c(1, 2), c(2, 3)))
  expect_true(setequal(letters[1:3], rev(letters[1:3])))
})


test_that("setequal() decides convex regions from their facets alone", {
  # The same halfspace written two ways. Scaling a normal changes nothing about
  # the set, and nothing here is decomposed to find that out.
  expect_true(setequal(
    halfspace_region(normal = c(1, -1, 0)),
    halfspace_region(normal = c(2, -2, 0), offset = 0)
  ))
  expect_false(setequal(
    halfspace_region(normal = c(1, -1, 0)),
    halfspace_region(normal = c(1, -1, 0), offset = 1)
  ))

  # Reflexivity across every geometry, including the ones whose *generators*
  # are a rounded frame rather than an exact description of them.
  for (region in list(
    simplex_region(vertices = diag(3)),
    polytope_region(vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))),
    halfspace_region(normal = c(1, -1, 0)),
    point_region(theta = c(0.5, 0.3, 0.2)),
    unconstrained_region(2L)
  )) {
    expect_true(setequal(region, region))
  }

  expect_false(setequal(
    simplex_region(vertices = diag(3)),
    unconstrained_region(3L)
  ))
  # A different ambient dimension is a `FALSE`, not an error: two sets in
  # different spaces are answerably not the same set.
  expect_false(setequal(
    simplex_region(vertices = diag(3)),
    simplex_region(vertices = diag(4))
  ))
})


test_that("setequal() sees through a decomposition into cells", {
  # The case an emptiness test gets wrong. `setdiff(square, its triangles)` is
  # the diagonal, not nothing, because the difference is closed; the square is
  # nonetheless exactly the union of its triangles.
  square <- polytope_region(
    vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  expect_length(cells(square), 2L)
  expect_false(is.null(setdiff(square, union_region(cells(square)))))
  expect_true(setequal(square, union_region(cells(square))))

  # And one triangle alone is not the square, which is what stops the
  # dimension test from calling everything equal.
  expect_false(setequal(square, cells(square)[[1L]]))
})


test_that("setequal() puts a region back together from its complement", {
  ambient <- simplex_region(vertices = diag(3))
  null_region <- union_region(plurality_cell(3L, 2L), plurality_cell(3L, 3L))
  wins <- setdiff(ambient, null_region)

  expect_true(setequal(union(null_region, wins), ambient))
  expect_false(setequal(null_region, ambient))
  expect_false(setequal(wins, ambient))

  # The overlapping declared cover and the peeled one are the same set.
  expect_true(setequal(
    union(plurality_cell(3L, 2L), setdiff(ambient, plurality_cell(3L, 2L))),
    ambient
  ))
})


test_that("setequal() does not inherit setdiff's slice warning", {
  # A part meeting the ambient in a slice makes `setdiff()` warn that it
  # subtracted nothing. For an equality test that is not news: under-
  # subtracting leaves the difference larger, and a slice was never going to
  # cover a full-dimensional piece anyway.
  square <- polytope_region(
    vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  diagonal <- polytope_region(vertices = cbind(c(0, 0), c(1, 1)))
  padded <- union_region(square, diagonal)

  expect_warning(setdiff(square, padded), "lower-dimensional")
  expect_silent(expect_true(setequal(square, padded)))
})


# --- disjoin ------------------------------------------------------------------

test_that("disjoin() covers the same set with parts that do not overlap", {
  plurality <- union_region(plurality_cell(3L, 2L), plurality_cell(3L, 3L))
  peeled <- disjoin(plurality)

  expect_true(setequal(peeled, plurality))
  # The declared parts genuinely overlap; the peeled ones meet in nothing of
  # full dimension, which is what lets a measure be summed over them.
  expect_false(is.null(intersect(
    plurality_cell(3L, 2L),
    plurality_cell(3L, 3L)
  )))
  for (i in seq_len(n_parts(peeled) - 1L)) {
    for (j in (i + 1L):n_parts(peeled)) {
      shared <- intersect(peeled[[i]], peeled[[j]])
      if (is.null(shared)) {
        next
      }
      # They may meet, but only in a face: the plurality cells are
      # 2-dimensional in `R^3`, so an overlap of dimension 2 would be the
      # double-counting the peeling exists to remove.
      full <- min(q_dim(q_hrep(peeled[[i]])), q_dim(q_hrep(peeled[[j]])))
      for (cell in parts(shared)) {
        expect_lt(q_dim(q_hrep(cell)), full)
      }
    }
  }
})


test_that("disjoin() of a convex region is that region", {
  s <- plurality_cell(3L, 2L)
  expect_identical(disjoin(s), s)
})


test_that("disjoin() drops a part its predecessors already cover", {
  ambient <- simplex_region(vertices = diag(3))
  # The second part is inside the first, so it survives as nothing at all.
  peeled <- disjoin(union_region(ambient, plurality_cell(3L, 2L)))
  expect_identical(peeled, ambient)
  expect_true(setequal(peeled, ambient))
})
