# Properties of the `union_region` half of R/region.R: the container for a
# finite union of convex cells, and the part/cell split every region answers.
#
# The invariant under test throughout is that a union records what was
# *declared*. It flattens nesting, because nesting carries no information, and
# it changes nothing else -- in particular it never merges or drops an overlap.

# The K = 3 plurality null: {theta_1 <= theta_2} and {theta_1 <= theta_3}
# as sub-simplices of the standard 2-simplex. They overlap along the tie.
plurality_cell <- function(k, j) {
  basis <- lapply(setdiff(seq_len(k), 1L), function(i) {
    v <- numeric(k)
    v[i] <- 1
    v
  })
  tie <- numeric(k)
  tie[c(1L, j)] <- 0.5
  simplex_region(vertices = do.call(cbind, c(basis, list(tie))))
}

# Every geometry the package currently has, for the sweeps below.
every_geometry <- function() {
  list(
    polytope = polytope_region(
      vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
    ),
    simplex = simplex_region(vertices = diag(3)),
    halfspace = halfspace_region(normal = c(1, -1, 0), offset = 0),
    point = point_region(theta = c(0.5, 0.3, 0.2)),
    unconstrained = unconstrained_region(2L),
    union = union_region(plurality_cell(3, 2), plurality_cell(3, 3))
  )
}

# --- The hierarchy ------------------------------------------------------------

test_that("every convex geometry is both a region and a convex_region", {
  for (g in every_geometry()[-6L]) {
    expect_true(S7_inherits(g, region))
    expect_true(S7_inherits(g, convex_region))
  }
})

test_that("a union is a region but not a convex one", {
  # The whole point of the split: `chart()`, `project()` and `maximise_over()`
  # dispatch on or type as `convex_region`, so a union cannot reach them.
  u <- union_region(plurality_cell(3, 2), plurality_cell(3, 3))
  expect_true(S7_inherits(u, region))
  expect_false(S7_inherits(u, convex_region))
})

# --- Construction -------------------------------------------------------------

test_that("cells, lists, unions and any nesting of them all flatten alike", {
  a <- simplex_region(vertices = diag(3))
  b <- halfspace_region(normal = c(1, -1, 0), offset = 0)

  built <- list(
    union_region(a, b),
    union_region(list(a, b)),
    union_region(union_region(a), b),
    union_region(list(union_region(a), list(b))),
    union_region(union_region(a, b)),
    union_region(list(list(a), list(list(b))))
  )
  for (r in built) {
    expect_equal(n_parts(r), 2L)
    expect_identical(parts(r), list(a, b))
  }
})

test_that("one convex region is handed back unwrapped", {
  # A single part is already a region, so a one-element union would be a
  # distinction without a difference for every caller.
  s <- simplex_region(vertices = diag(3))
  expect_identical(union_region(s), s)
  expect_identical(union_region(list(s)), s)
  expect_identical(union_region(union_region(s)), s)
  expect_false(S7_inherits(union_region(s), union_region))
})

test_that("the decomposition caches start empty", {
  r <- union_region(plurality_cell(3, 2), plurality_cell(3, 3))
  expect_null(r@disjoint)
  expect_null(r@triangulation)
})

# --- Validation ---------------------------------------------------------------

test_that("cells of different ambient dimensions are refused, naming both", {
  err <- expect_error(
    union_region(simplex_region(vertices = diag(3)), unconstrained_region(2L)),
    "same ambient dimension"
  )
  expect_match(conditionMessage(err), "3")
  expect_match(conditionMessage(err), "2")
})

test_that("cells of differing shape and codimension are welcome", {
  # Only the ambient dimension is constrained. A triangle (intrinsic 2), a
  # point (0) and a halfspace (3) all sit in R^3 and union without complaint.
  u <- union_region(
    simplex_region(vertices = diag(3)),
    point_region(theta = c(1, 0, 0)),
    halfspace_region(normal = c(1, -1, 0), offset = 0)
  )
  expect_equal(n_parts(u), 3L)
  expect_equal(space_dim(u), 3L)
  expect_true(contains(u, rep(1 / 3, 3)))
  expect_true(contains(u, c(1, 0, 0)))
})

test_that("an empty union is refused", {
  expect_error(union_region(), "non-empty")
  expect_error(union_region(list()), "non-empty")
})

test_that("a non-convex_region element is refused", {
  expect_error(union_region("not a region"), "must be a `convex_region`")
  expect_error(
    union_region(simplex_region(vertices = diag(3)), 42),
    "must be a `convex_region`"
  )
  # An S7 object that is not a cell is a leaf, not something to descend into.
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  expect_error(union_region(fam), "must be a `convex_region`")
})

# --- parts() and cells() ------------------------------------------------------

test_that("a convex region is its own only part and its own only cell", {
  s <- simplex_region(vertices = diag(3))
  expect_identical(parts(s), list(s))
  expect_identical(cells(s), list(s))
  expect_equal(n_parts(s), 1L)
  expect_equal(n_cells(s), 1L)
})

test_that("parts() and cells() agree on every geometry in the package", {
  # A marker for the phase that adds triangulation. When we add triangulation
  # we will get more cells than parts, this test will fail for the polytope
  # row, and that failure is a reminder to check every caller picked the
  # accessor it meant. Until then, they should be the same.
  for (nm in names(every_geometry())) {
    g <- every_geometry()[[nm]]
    expect_identical(parts(g), cells(g), info = nm)
    expect_equal(n_parts(g), n_cells(g), info = nm)
  }
})

test_that("a union's parts are its members and its cells are theirs, flattened", {
  a <- plurality_cell(3, 2)
  b <- plurality_cell(3, 3)
  r <- union_region(a, b)
  expect_identical(parts(r), list(a, b))
  expect_identical(cells(r), list(a, b))
})

test_that("a union reports the dimension its cells share", {
  expect_equal(
    space_dim(union_region(plurality_cell(4, 2), plurality_cell(4, 3))),
    4L
  )
})

# --- Membership ---------------------------------------------------------------

test_that("contains() on a union is membership of any one cell", {
  r <- union_region(plurality_cell(3, 2), plurality_cell(3, 3))

  # theta_1 <= theta_2 but not theta_1 <= theta_3: the first cell only.
  expect_true(contains(r, c(0.3, 0.5, 0.2)))
  expect_true(contains(plurality_cell(3, 2), c(0.3, 0.5, 0.2)))
  expect_false(contains(plurality_cell(3, 3), c(0.3, 0.5, 0.2)))

  # The barycentre is in the overlap.
  expect_true(contains(r, rep(1 / 3, 3)))

  # theta_1 the strict plurality winner: outside both.
  expect_false(contains(r, c(0.6, 0.2, 0.2)))
})

# --- Overlaps -----------------------------------------------------------------

test_that("overlapping cells survive construction untouched", {
  # This test is just a reminder that we don't necessarily want the parts to be
  # disjoint.
  s <- simplex_region(vertices = diag(3))
  expect_equal(n_parts(union_region(s, s)), 2L)
  expect_identical(parts(union_region(s, s)), list(s, s))

  # Genuinely overlapping, distinct cells are equally untouched.
  r <- union_region(plurality_cell(3, 2), plurality_cell(3, 3))
  expect_equal(n_parts(r), 2L)
  expect_true(contains(parts(r)[[1L]], rep(1 / 3, 3)))
  expect_true(contains(parts(r)[[2L]], rep(1 / 3, 3)))
})

# --- Printing -----------------------------------------------------------------

test_that("print() and format() describe the union", {
  r <- union_region(plurality_cell(3, 2), plurality_cell(3, 3))

  expect_output(print(r), "union_region")
  expect_output(print(r), "2 parts")
  expect_output(print(r), "simplex_region")

  expect_type(format(r), "character")
  expect_length(format(r), 1L)
  expect_match(format(r), "2 parts")
  expect_match(format(r), "dimension 3")

  # Beyond a handful the cells are tallied rather than listed one by one: a
  # triangulated null can hold hundreds, and the list would say nothing extra.
  many <- union_region(lapply(2:8, \(j) plurality_cell(8, j)))
  expect_equal(n_parts(many), 7L)
  expect_output(print(many), "7 parts")
  expect_output(print(many), "7 x simplex_region")
})
