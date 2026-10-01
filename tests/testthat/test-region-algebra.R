# Tests for the set algebra in R/region_algebra.R.
#
# The algebra composes on the rational layer and converts once, which is what
# the exactness assertions here pin: derived coordinates equal their rational
# values to the last bit, not to a tolerance.

# One cell of the K-candidate plurality null: candidate 1 trails candidate j.
plurality_cell <- function(k, j) {
  vertices <- diag(k)
  vertices[1L, ] <- replace(numeric(k), c(1L, j), 0.5)
  simplex_region(vertices = vertices)
}

# Rows sorted lexicographically, so vertex sets compare independently of
# the order cddlib emits them in.
sorted_rows <- function(m) {
  m[do.call(order, asplit(m, 2L)), , drop = FALSE]
}

# The dimension of a convex region's affine hull, read off its generators: the
# rank of the vertex differences together with the rays and lines.
affine_dim <- function(cell) {
  g <- cell@generators
  spread <- rbind(
    sweep(g$v, 2L, g$v[1L, ])[-1L, , drop = FALSE],
    g$r,
    g$l
  )
  if (nrow(spread) == 0L) {
    return(0L)
  }
  qr(spread, tol = 1e-9)$rank
}

# Evaluate `code` with the `ripr.max_cells` option set to `n`.
with_max_cells <- function(n, code) {
  old <- options(ripr.max_cells = n)
  on.exit(options(old))
  code
}

# --- Operators and base R ----------------------------------------------------

test_that("ripr leaves base R's set functions alone", {
  # The algebra is the operators on `region`; the base verbs are not masked.
  expect_false(any(
    c("union", "intersect", "setdiff", "setequal") %in%
      getNamespaceExports("ripr")
  ))
  expect_identical(union(c(1, 2), c(2, 3)), c(1, 2, 3))
  expect_identical(intersect(c(1, 2), c(2, 3)), 2)
  expect_identical(setdiff(c(1, 2), c(2, 3)), 1)
  expect_true(setequal(c(1, 2), c(2, 1)))
})


test_that("the operators still mean what base R says off regions", {
  expect_identical(c(TRUE, FALSE) | c(FALSE, FALSE), c(TRUE, FALSE))
  expect_identical(c(TRUE, TRUE) & c(FALSE, TRUE), c(FALSE, TRUE))
  expect_identical(3 - 1, 2)
  expect_identical(c(1, 2) == c(1, 3), c(TRUE, FALSE))
})


test_that("a region does not combine with a non-region", {
  s <- simplex_region(vertices = diag(3))
  expect_error(s - 1)
  expect_error(s | TRUE)
  expect_error(s == 1)
})


# --- Intersection -------------------------------------------------------------

test_that("two overlapping plurality cells intersect in the exact triangle", {
  ab <- plurality_cell(3L, 2L) & plurality_cell(3L, 3L)

  # {theta1 <= theta2} meets {theta1 <= theta3} where candidate 1 trails both:
  # the triangle spanned by the two loser vertices and the barycentre. Its
  # coordinates are exactly 0, 1 and double(1/3) -- the acceptance test for
  # composing in rationals and converting once.
  expect_true(S7_inherits(ab, polytope_region))
  expect_identical(
    sorted_rows(ab@vertices),
    sorted_rows(rbind(c(0, 1, 0), c(0, 0, 1), c(1, 1, 1) / 3))
  )
})


test_that("`&` distributes over the parts of unions", {
  u <- union_region(plurality_cell(3L, 2L), plurality_cell(3L, 3L))
  h <- halfspace_region(normal = c(0, 1, -1)) # theta2 <= theta3

  r <- u & h
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
  nothing <- point_region(theta = c(1, 0, 0)) &
    point_region(theta = c(0, 1, 0))
  expect_true(S7_inherits(nothing, empty_region))

  # Cells meeting only in a shared face are not empty: the closed cells of a
  # cover genuinely intersect in that face.
  edge <- simplex_region(vertices = rbind(c(0, 0), c(1, 0), c(0, 1))) &
    simplex_region(vertices = rbind(c(1, 1), c(1, 0), c(0, 1)))
  expect_false(is_empty(edge))
  expect_true(contains(edge, c(0.5, 0.5)))
  expect_false(contains(edge, c(0.25, 0.25)))
})


test_that("an unbounded intersection returns a polyhedron_region", {
  quadrant <- halfspace_region(normal = c(1, 0)) &
    halfspace_region(normal = c(0, 1))
  expect_true(S7_inherits(quadrant, polyhedron_region))
  expect_false(is_bounded(quadrant))
  expect_true(contains(quadrant, c(-3, -5)))
  expect_false(contains(quadrant, c(1, 1)))
})


test_that("`&` chains over several regions and refuses mismatched dimensions", {
  # Three halfspaces cutting the plane down to a bounded triangle.
  tri <- halfspace_region(normal = c(-1, 0)) & # x >= 0
    halfspace_region(normal = c(0, -1)) & # y >= 0
    halfspace_region(normal = c(1, 1), offset = 1) # x + y <= 1
  expect_true(is_bounded(tri))
  expect_identical(
    sorted_rows(tri@generators$v),
    sorted_rows(rbind(c(0, 0), c(1, 0), c(0, 1)))
  )

  expect_error(
    real_region(2L) & real_region(3L),
    "ambient dimension"
  )
})


# --- Difference ---------------------------------------------------------------

test_that("subtracting a union subtracts each of its parts", {
  ambient <- simplex_region(vertices = diag(3))
  at_once <- ambient -
    union_region(plurality_cell(3L, 2L), plurality_cell(3L, 3L))
  in_turn <- (ambient - plurality_cell(3L, 2L)) - plurality_cell(3L, 3L)
  expect_true(at_once == in_turn)
  by_operator <- ambient - (plurality_cell(3L, 2L) | plurality_cell(3L, 3L))
  expect_true(at_once == by_operator)
})


test_that("subtracting one plurality cell from the simplex gives one part", {
  # Every facet of the cell except its own boundary is a wall of the simplex,
  # dropped by the exact LP; what remains is the single cell where candidate 1
  # beats candidate j.
  for (k in c(3L, 4L, 5L)) {
    ambient <- simplex_region(vertices = diag(k))
    left <- ambient - plurality_cell(k, 2L)
    expect_true(S7_inherits(left, convex_region))
    expect_identical(length(parts(left)), 1L)
  }
})


test_that("the complement of the plurality null is the candidate-1-wins region", {
  ambient <- simplex_region(vertices = diag(3))
  null_region <- union_region(plurality_cell(3L, 2L), plurality_cell(3L, 3L))
  wins <- ambient - null_region

  expect_true(contains(wins, c(0.5, 0.3, 0.2)))
  expect_false(contains(wins, c(0.2, 0.5, 0.3), tol = 1e-12))

  # Every cell of a difference inside the simplex still lives on the simplex:
  # exactly one equality row, and every vertex sums to one.
  for (cell in parts(wins)) {
    h <- cell@facets
    expect_identical(sum(h$eq), 1L)
    expect_true(max(abs(rowSums(cell@generators$v) - 1)) <= rounding_tol(1))
  }
})


test_that("K = 5 complement of the whole plurality null is fast and small", {
  ambient <- simplex_region(vertices = diag(5L))
  null_region <- union_region(lapply(2:5, \(j) plurality_cell(5L, j)))
  elapsed <- system.time(wins <- ambient - null_region)[["elapsed"]]
  expect_lt(elapsed, 1)
  expect_identical(length(parts(wins)), 1L)
  expect_true(contains(wins, c(0.6, 0.1, 0.1, 0.1, 0.1)))
})


test_that("the double difference agrees with the original", {
  ambient <- simplex_region(vertices = diag(3))
  null_region <- union_region(plurality_cell(3L, 2L), plurality_cell(3L, 3L))
  back <- ambient - (ambient - null_region)

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
    back <- real_region(3L) - simplex_region(vertices = diag(3)),
    "lower-dimensional"
  )
  expect_identical(length(parts(back)), 1L)
  expect_true(contains(back, c(5, -3, 2)))

  # A shared affine hull is fine: both live on the simplex.
  expect_no_warning(
    simplex_region(vertices = diag(3)) - plurality_cell(3L, 2L)
  )
})


test_that("a coarser ambient gives strictly more cells", {
  # The triangle shares two edges with the small square but none with the big
  # one, so fewer of its facets are dropped as ambient-implied.
  triangle <- simplex_region(vertices = rbind(c(0, 0), c(1, 0), c(0, 1)))
  small <- polytope_region(vertices = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)))
  big <- polytope_region(
    vertices = 4 * rbind(c(-1, -1), c(1, -1), c(1, 1), c(-1, 1))
  )

  n_small <- length(parts(small - triangle))
  n_big <- length(parts(big - triangle))
  expect_gt(n_big, n_small)
})


test_that("x inside y leaves nothing, and max_cells errors rather than hangs", {
  nothing <- plurality_cell(3L, 2L) - simplex_region(vertices = diag(3))
  expect_true(S7_inherits(nothing, empty_region))

  big <- polytope_region(
    vertices = 4 * rbind(c(-1, -1), c(1, -1), c(1, 1), c(-1, 1))
  )
  squares <- union_region(lapply(0:2, function(k) {
    polytope_region(
      vertices = rbind(c(k, 0), c(k + 0.5, 0), c(k + 0.5, 0.5), c(k, 0.5))
    )
  }))
  # Six cells in all.
  with_max_cells(5L, expect_error(big - squares, "max_cells = 5"))
})


test_that("difference cells are interior-disjoint by construction", {
  big <- polytope_region(
    vertices = 2 * rbind(c(-1, -1), c(1, -1), c(1, 1), c(-1, 1))
  )
  triangle <- simplex_region(vertices = rbind(c(0, 0), c(1, 0), c(0, 1)))
  left <- big - triangle

  # Interior points of one cell belong to no other cell.
  set.seed(17)
  for (cell in parts(left)) {
    ch <- chart(cell)
    for (i in seq_len(20L)) {
      theta <- ch$to_theta(stats::rnorm(ch$n_par))
      others <- Filter(\(p) !identical(p, cell), parts(left))
      strictly_inside <- all(vapply(
        seq_len(nrow(cell@facets$a)),
        \(r) {
          h <- cell@facets
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


test_that("algebra cells keep their exact equality rows", {
  # A produced cell's facets must be the exact rows, not re-derived from the
  # rounded vertices: a quadrilateral on the simplex with a 1/3 vertex rounds
  # to four points that exact arithmetic sees as an ulp-thin tetrahedron --
  # full-dimensional, no equality row, sliver facets.
  ambient <- simplex_region(vertices = diag(3))
  null_region <- union_region(plurality_cell(3L, 2L), plurality_cell(3L, 3L))
  wins <- ambient - null_region

  expect_identical(sum(wins@facets$eq), 1L)
  expect_identical(affine_dim(wins), 2L)

  # A declared region on the simplex carries its equality row likewise.
  expect_identical(sum(ambient@facets$eq), 1L)
})


test_that("algebra cells that are simplices come back as simplex_region", {
  # The K = 3 single-cell complement has three vertices on the simplex: a
  # certifiable simplex, and classified as one.
  left <- simplex_region(vertices = diag(3)) - plurality_cell(3L, 2L)
  expect_true(S7_inherits(left, simplex_region))

  # The intersection triangle of the two plurality cells likewise.
  ab <- plurality_cell(3L, 2L) & plurality_cell(3L, 3L)
  expect_true(S7_inherits(ab, simplex_region))

  # The candidate-1-wins region is a quadrilateral: a polytope, not a simplex.
  wins <- simplex_region(vertices = diag(3)) -
    union_region(plurality_cell(3L, 2L), plurality_cell(3L, 3L))
  expect_false(S7_inherits(wins, simplex_region))
  expect_true(S7_inherits(wins, polytope_region))
})


test_that("intersection reduces dimension; only a difference refuses slices", {
  # The positive cone meets the sum-one hyperplane in the probability
  # simplex, exactly.
  cone <- polyhedron_region(rays = diag(3))
  plane <- h_region(a = matrix(1, 1L, 3L), b = 1, eq = TRUE)
  simplex <- cone & plane
  expect_true(S7_inherits(simplex, simplex_region))
  expect_identical(
    sorted_rows(simplex@vertices),
    sorted_rows(diag(3) + 0)
  )
})


test_that("a subtracted part that never meets x subtracts nothing", {
  square <- polytope_region(
    vertices = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  # A segment on a hyperplane elsewhere: lower-dimensional, but disjoint from
  # the square, so the difference is the square itself rather than a refusal.
  far_segment <- simplex_region(vertices = rbind(c(5, 0), c(5, 1)))
  left <- square - far_segment
  expect_identical(length(parts(left)), 1L)
  expect_true(contains(left, c(0.5, 0.5)))

  # A slice actually through the square warns and leaves the square whole.
  through <- simplex_region(vertices = rbind(c(0.5, 0), c(0.5, 1)))
  expect_warning(whole <- square - through, "lower-dimensional")
  expect_identical(length(parts(whole)), 1L)
  expect_true(contains(whole, c(0.5, 0.5)))

  # And a disjoint full-dimensional part adds no cells.
  far_square <- polytope_region(
    vertices = rbind(c(5, 5), c(6, 5), c(6, 6), c(5, 6))
  )
  triangle <- simplex_region(vertices = rbind(c(0, 0), c(1, 0), c(0, 1)))
  expect_identical(
    length(parts(square - union_region(triangle, far_square))),
    length(parts(square - triangle))
  )
})


# --- Equality -----------------------------------------------------------------

test_that("`==` decides convex regions from their facets alone", {
  # The same halfspace written two ways. Scaling a normal changes nothing about
  # the set, and nothing here is decomposed to find that out.
  h <- halfspace_region(normal = c(1, -1, 0))
  expect_true(h == halfspace_region(normal = c(2, -2, 0), offset = 0))
  expect_false(h == halfspace_region(normal = c(1, -1, 0), offset = 1))
  expect_true(h != halfspace_region(normal = c(1, -1, 0), offset = 1))

  # Reflexivity across every geometry, including the ones whose *generators*
  # are a rounded frame rather than an exact description of them.
  for (region in list(
    simplex_region(vertices = diag(3)),
    polytope_region(vertices = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))),
    halfspace_region(normal = c(1, -1, 0)),
    point_region(theta = c(0.5, 0.3, 0.2)),
    real_region(2L)
  )) {
    expect_true(region == region)
  }

  expect_false(simplex_region(vertices = diag(3)) == real_region(3L))
  # A different ambient dimension is a `FALSE`, not an error: two sets in
  # different spaces are answerably not the same set.
  expect_false(
    simplex_region(vertices = diag(3)) == simplex_region(vertices = diag(4))
  )
})


test_that("`==` distinguishes a region from its boundary face", {
  # The face's equality row is what rejects the body, and it has to be tested
  # in both directions: the square satisfies `y <= 0` everywhere, so only the
  # reverse `y >= 0` can turn it away.
  square <- polytope_region(
    vertices = rbind(c(0, -1), c(1, -1), c(1, 0), c(0, 0))
  )
  edge <- polytope_region(vertices = rbind(c(0, 0), c(1, 0)))
  expect_false(square == edge)
  expect_false(edge == square)

  # And from the other side of the hyperplane, so whichever way cddlib
  # orients the equality row, both directions of the test get exercised.
  above <- polytope_region(
    vertices = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  expect_false(above == edge)
})


test_that("`-` refuses mismatched ambient dimensions", {
  expect_error(
    polytope_region(vertices = rbind(c(0, 0), c(1, 0), c(0, 1))) -
      simplex_region(vertices = diag(3)),
    "ambient dimension"
  )
})


test_that("`==` sees through a decomposition into cells", {
  square <- polytope_region(
    vertices = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  expect_length(cells(square), 2L)
  expect_true(is_empty(square - union_region(cells(square))))
  expect_true(square == union_region(cells(square)))

  # And one triangle alone is not the square.
  expect_false(square == cells(square)[[1L]])
})


test_that("`==` puts a region back together from its complement", {
  ambient <- simplex_region(vertices = diag(3))
  null_region <- union_region(plurality_cell(3L, 2L), plurality_cell(3L, 3L))
  wins <- ambient - null_region

  expect_true((null_region | wins) == ambient)
  expect_false(null_region == ambient)
  expect_false(wins == ambient)

  # The overlapping declared cover and the peeled one are the same set.
  expect_true(
    (plurality_cell(3L, 2L) | (ambient - plurality_cell(3L, 2L))) == ambient
  )
})


test_that("`==` does not inherit the slice warning of `-`", {
  # A part meeting the ambient in a slice makes `-` warn that it
  # subtracted nothing. For an equality test that is not news: under-
  # subtracting leaves the difference larger, and a slice was never going to
  # cover a full-dimensional piece anyway.
  square <- polytope_region(
    vertices = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  diagonal <- polytope_region(vertices = rbind(c(0, 0), c(1, 1)))
  padded <- union_region(square, diagonal)

  expect_warning(square - padded, "lower-dimensional")
  expect_silent(expect_true(square == padded))
})


# --- disjoin ------------------------------------------------------------------

test_that("disjoin() covers the same set with parts that do not overlap", {
  plurality <- union_region(plurality_cell(3L, 2L), plurality_cell(3L, 3L))
  peeled <- disjoin(plurality)

  expect_true(peeled == plurality)
  # The declared parts genuinely overlap; the peeled ones meet in nothing of
  # full dimension, which is what lets a measure be summed over them.
  expect_false(is_empty(plurality_cell(3L, 2L) & plurality_cell(3L, 3L)))
  for (i in seq_len(length(parts(peeled)) - 1L)) {
    for (j in (i + 1L):length(parts(peeled))) {
      shared <- parts(peeled)[[i]] & parts(peeled)[[j]]
      if (is_empty(shared)) {
        next
      }
      # They may meet, but only in a face: the plurality cells are
      # 2-dimensional in `R^3`, so an overlap of dimension 2 would be the
      # double-counting the peeling exists to remove.
      full <- min(affine_dim(parts(peeled)[[i]]), affine_dim(parts(peeled)[[j]]))
      for (cell in parts(shared)) {
        expect_lt(affine_dim(cell), full)
      }
    }
  }
})


test_that("disjoin() drops a part its predecessors already cover", {
  ambient <- simplex_region(vertices = diag(3))
  # The second part is inside the first, so it survives as nothing at all.
  peeled <- disjoin(union_region(ambient, plurality_cell(3L, 2L)))
  expect_identical(peeled, ambient)
  expect_true(peeled == ambient)
})


test_that("disjoin() is silent about a part that only slices another", {
  # Peeling the square away from the diagonal subtracts a slice, which `-`
  # would warn about; disjoin() is subtracting only to build a cover, so it
  # says nothing and keeps both.
  square <- polytope_region(
    vertices = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  diagonal <- polytope_region(vertices = rbind(c(0, 0), c(1, 1)))
  expect_silent(peeled <- disjoin(union_region(diagonal, square)))
  expect_true(peeled == square)
})


test_that("the difference cap is the ripr.max_cells option", {
  big <- polytope_region(
    vertices = 4 * rbind(c(-1, -1), c(1, -1), c(1, 1), c(-1, 1))
  )
  squares <- union_region(lapply(0:2, function(k) {
    polytope_region(
      vertices = rbind(c(k, 0), c(k + 0.5, 0), c(k + 0.5, 0.5), c(k, 0.5))
    )
  }))
  with_max_cells(5L, {
    expect_error(big - squares, "max_cells = 5")
    expect_error(disjoin(union_region(squares, big)), "max_cells = 5")
    # An explicit argument outranks the option.
    expect_no_error(disjoin(union_region(squares, big), max_cells = 1000L))
  })
  with_max_cells(1000L, expect_no_error(big - squares))
})
