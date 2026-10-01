# Properties of R/region.R and R/null.R.
#
# The geometry is checked by round-trips and idempotence rather than against
# stored coordinates, since a chart is only required to be *a* parametrisation,
# not a particular one.

# The plurality part {theta : theta_1 <= theta_j} on the K-simplex, in both
# representations, so the two can be checked against each other.
plurality_simplex <- function(k, j) {
  basis <- lapply(setdiff(seq_len(k), 1L), function(i) {
    v <- numeric(k)
    v[i] <- 1
    v
  })
  tie <- numeric(k)
  tie[c(1L, j)] <- 0.5
  simplex_region(vertices = do.call(rbind, c(basis, list(tie))))
}

plurality_halfspace <- function(k, j) {
  a <- numeric(k)
  a[1L] <- 1
  a[j] <- -1
  halfspace_region(normal = a, offset = 0)
}

# --- Charts -------------------------------------------------------------------

test_that("a chart round-trips points in the part", {
  for (s in list(plurality_simplex(4, 2), plurality_halfspace(4, 2))) {
    ch <- chart(s)
    set.seed(1)
    for (i in 1:5) {
      u <- ch$seed(1L)[1L, ]
      theta <- ch$to_theta(u)
      expect_true(contains(s, theta))
      expect_equal(
        ch$to_theta(ch$from_theta(theta)),
        theta,
        tolerance = rounding_tol(1)
      )
    }
  }
})

test_that("the chart Jacobian matches a finite difference", {
  for (s in list(plurality_simplex(4, 2), plurality_halfspace(3, 2))) {
    ch <- chart(s)
    set.seed(3)
    u <- ch$seed(1L)[1L, ]
    jac <- ch$jacobian(u)
    eps <- 1e-6
    for (j in seq_len(ch$n_par)) {
      e <- numeric(ch$n_par)
      e[j] <- eps
      fd <- (ch$to_theta(u + e) - ch$to_theta(u - e)) / (2 * eps)
      expect_equal(jac[, j], fd, tolerance = 1e-5)
    }
  }
})

# --- Membership and projection ------------------------------------------------

test_that("projection is idempotent and lands in the part", {
  set.seed(5)
  for (s in list(plurality_simplex(4, 2), plurality_halfspace(4, 2))) {
    for (i in 1:5) {
      theta <- stats::runif(4)
      theta <- theta / sum(theta)
      p <- project(s, theta)
      expect_true(contains(s, p))
      expect_equal(project(s, p), p, tolerance = rounding_tol(1))
    }
  }
})

test_that("projection leaves points already inside untouched", {
  s <- plurality_halfspace(3, 2)
  theta <- c(0.2, 0.5, 0.3) # theta_1 <= theta_2
  expect_true(contains(s, theta))
  expect_equal(project(s, theta), theta)
})

test_that("the two representations agree on membership within the simplex", {
  # The vertex hull and the halfspace describe the same set once intersected
  # with the probability simplex, so they must agree on simplex points.
  set.seed(6)
  hull <- plurality_simplex(4, 2)
  half <- plurality_halfspace(4, 2)
  for (i in 1:40) {
    theta <- stats::rgamma(4, 1)
    theta <- theta / sum(theta)
    expect_equal(
      contains(hull, theta, tol = 1e-6),
      contains(half, theta, tol = 1e-6)
    )
  }
})

test_that("a simplex part contains its own vertices", {
  s <- plurality_simplex(4, 3)
  for (i in seq_len(nrow(s@vertices))) {
    expect_true(contains(s, s@vertices[i, ], tol = 1e-6))
  }
})

# --- The search ---------------------------------------------------------------
#
# Every oracle in the package is the same multistart search over a part's
# chart; `sup_lb()` is the one that takes an arbitrary objective, as the mean
# of a random variable. Under the multinomial `E[x_i] = n theta_i` and
# `E[x_i (x_i - 1)] = n (n - 1) theta_i^2`, so the variable below has mean
# exactly `1 - ||theta - target||^2`: a concave peak placed wherever a test
# needs it, with maximum 1.

peak_at <- function(fam, target) {
  n <- fam@n_trials
  random_variable(
    function(x) {
      x <- matrix(x, ncol = length(target))
      1 - rowSums(x * (x - 1)) / (n * (n - 1)) +
        2 * as.vector(x %*% target) / n - sum(target^2)
    },
    fam@sample_space
  )
}

test_that("the search finds a maximum interior to the chart", {
  # The peak is interior in *vertex-weight* coordinates, not merely inside the
  # set.
  fam <- multinomial_family(n_trials = 6L, k = 3L)
  s <- plurality_simplex(3, 2)
  target <- as.vector(rep(1 / 3, 3) %*% s@vertices)
  set.seed(7)
  res <- sup_lb(
    peak_at(fam, target),
    null_model(fam, s),
    n_seeds = 50L,
    n_restarts = 5L
  )
  expect_equal(res@theta, target, tolerance = 1e-5)
  expect_equal(res@sup_lb, 1, tolerance = 1e-9)
  expect_true(contains(s, res@theta, tol = 1e-6))
})

test_that("a maximum at a vertex is attained exactly", {
  # SLSQP's active constraints pin a vertex maximum exactly. The objective here
  # is concave, so the multistart is guaranteed the right basin and the test
  # is deterministic; for a non-convex oracle the exactness holds only once
  # the search finds the right face, and the search remains a lower bound on
  # the supremum.
  fam <- multinomial_family(n_trials = 6L, k = 3L)
  s <- plurality_simplex(3, 2)
  target <- s@vertices[1L, ] # vertex weight (1, 0, 0)
  set.seed(12)
  res <- sup_lb(
    peak_at(fam, target),
    null_model(fam, s),
    n_seeds = 50L,
    n_restarts = 5L
  )
  expect_equal(res@sup_lb, 1, tolerance = rounding_tol(1))
  expect_equal(res@theta, target, tolerance = 1e-9)
  expect_true(contains(s, res@theta, tol = 1e-6))
})

test_that("the chart round-trip is lossless in the interior and at a vertex", {
  s <- plurality_simplex(4, 2)
  ch <- chart(s)

  interior <- as.vector(rep(0.25, 4) %*% s@vertices)
  expect_equal(
    ch$to_theta(ch$from_theta(interior)),
    interior,
    tolerance = rounding_tol(1)
  )

  # Exact up to the floating point of the least-squares recovery: bit-identical
  # on some vertex matrices, an ulp or two off on others, and BLAS-dependent
  # either way -- so tested at 1e-12, not identical().
  vertex <- s@vertices[1L, ]
  expect_equal(
    ch$to_theta(ch$from_theta(vertex)),
    vertex,
    tolerance = rounding_tol(1)
  )
})

test_that("the search accepts seeds lying outside the part", {
  # The fit seeds every part's search with every atom, and atoms live on other
  # parts, so seeds are projected before use. Here each part's atom lies
  # outside the other part, so both searches start from a foreign seed.
  fam <- multinomial_family(n_trials = 6L, k = 3L)
  s2 <- plurality_simplex(3, 2)
  s3 <- plurality_simplex(3, 3)
  a2 <- c(0.2, 0.7, 0.1)
  a3 <- c(0.2, 0.1, 0.7)
  expect_false(contains(s3, a2))
  expect_false(contains(s2, a3))
  set.seed(9)
  st <- ripr_init(
    mixture(fam, dirac(c(0.2, 0.6, 0.2))),
    null_model(fam, list(s2, s3)),
    atoms = list(a2, a3),
    record_gap = TRUE,
    control = ripr_control(n_seeds = 20L, n_restarts = 3L)
  )
  part <- list(s2, s3)[[st@oracle$part]]
  expect_true(contains(part, st@oracle$theta, tol = 1e-6))
  expect_true(is.finite(st@oracle$value))
})

test_that("the search over a singleton evaluates the point", {
  fam <- multinomial_family(n_trials = 4L, k = 2L)
  first <- random_variable(\(x) matrix(x, ncol = 2L)[, 1L], fam@sample_space)
  res <- sup_lb(first, null_model(fam, point_region(theta = c(0.5, 0.5))))
  expect_equal(res@theta, c(0.5, 0.5))
  expect_equal(res@sup_lb, 2)
})

# --- null_model ---------------------------------------------------------------

test_that("a null takes its geometry as a part, a list or a union alike", {
  # `null_model()` coerces, so none of the three needs the caller to know
  # which form the others take.
  fam <- multinomial_family(n_trials = 10, k = 3)
  s1 <- plurality_simplex(3, 2)
  s2 <- plurality_simplex(3, 3)
  theta <- c(0.3, 0.5, 0.2)

  from_list <- null_model(fam, list(s1, s2))
  from_union <- null_model(fam, union_region(s1, s2))
  expect_length(parts(from_list@region), 2L)
  expect_identical(from_list@region, from_union@region)
  expect_equal(
    contains(from_list@region, theta),
    contains(from_union@region, theta)
  )

  # A lone convex region is already a region, so it is stored as it came --
  # not wrapped in a one-element union.
  bare <- null_model(fam, s1)
  wrapped <- null_model(fam, list(s1))
  expect_identical(bare@region, s1)
  expect_identical(wrapped@region, s1)
  expect_length(parts(bare@region), 1L)
  expect_equal(contains(bare@region, theta), contains(wrapped@region, theta))
})

test_that("a null takes its decomposition once, at construction", {
  fam <- multinomial_family(n_trials = 10, k = 3)
  # A part that is already a simplex is its own cell, and the same object:
  # nothing has been rebuilt behind the caller's back.
  simplices <- null_model(fam, lapply(2:3, \(j) plurality_simplex(3, j)))
  expect_identical(simplices@cells, parts(simplices@region))
  expect_identical(simplices@cell_part, 1:2)

  # A part that is a convex hull is several cells, all filed under it.
  square <- polytope_region(
    vertices = rbind(
      c(0.5, 0.5, 0),
      c(0, 0.5, 0.5),
      c(0, 0, 1),
      c(0.5, 0, 0.5)
    )
  )
  mixed <- null_model(fam, list(plurality_simplex(3, 2), square))
  expect_length(mixed@cells, 3L)
  expect_identical(mixed@cell_part, c(1L, 2L, 2L))

  # The cells are in part order and are that part's own `cells()`.
  expect_identical(mixed@cells[[1L]], parts(mixed@region)[[1L]])
  expect_identical(mixed@cells[2:3], cells(square))
})

test_that("null_model()'s max_cells caps the decomposition across parts", {
  fam <- multinomial_family(n_trials = 10, k = 3)
  square <- polytope_region(
    vertices = rbind(
      c(0.5, 0.5, 0),
      c(0, 0.5, 0.5),
      c(0, 0, 1),
      c(0.5, 0, 0.5)
    )
  )
  # Two triangles per square, so two squares need four in total.
  null <- null_model(fam, list(square, square), max_cells = 4L)
  expect_length(null@cells, 4L)
  # One budget spans the parts, and the refusal still names the one that
  # exhausted it.
  expect_error(
    null_model(fam, list(square, square), max_cells = 3L),
    "could not decompose part 2.*max_cells = 3"
  )
})

test_that("simplex_region rejects a malformed vertex matrix", {
  expect_error(simplex_region(vertices = c(0.5, 0.5)), "must be a matrix")
  expect_error(
    simplex_region(vertices = matrix(numeric(0), nrow = 0L, ncol = 2L)),
    "must be a matrix"
  )
})

test_that("halfspace_region rejects a zero normal", {
  expect_error(halfspace_region(normal = c(0, 0)), "non-zero")
})

# --- Dimension ----------------------------------------------------------------

test_that("a region of the wrong dimension is refused at construction", {
  # Without this the comparison inside `contains()` recycles instead of
  # complaining, and `contains()` returns TRUE for a parameter it never checked:
  #   null_model(multinomial_family(4, 3), list(halfspace_region(c(1, -1))))
  #   contains(null@region, c(0.2, 0.5, 0.3))  # TRUE, silently wrong
  fam <- multinomial_family(n_trials = 4, k = 3)
  expect_error(
    null_model(fam, list(halfspace_region(normal = c(1, -1)))),
    "dimension 3"
  )
  expect_error(
    null_model(fam, list(point_region(theta = c(0.5, 0.5)))),
    "dimension 3"
  )
  # Parts that disagree with each other fail earlier and for a better reason:
  # `region`'s validator refuses them before the family is consulted at all,
  # since a union of a plane and a line has no ambient space to live in.
  expect_error(
    null_model(
      fam,
      list(simplex_region(vertices = diag(3)), real_region(2L))
    ),
    "same ambient dimension"
  )
  expect_silent(null_model(fam, list(real_region(3L))))
})

test_that("a family's parameter space has the family's own dimension", {
  # Unlike a per-family method, `space_dim()` on the parameter space cannot
  # disagree with the geometry the null is checked against.
  fam <- multinomial_family(n_trials = 7, k = 5)
  expect_equal(space_dim(fam@parameter_space), 5L)
  expect_true(contains(fam@parameter_space, rep(1 / 5, 5)))
  expect_false(contains(fam@parameter_space, c(0.5, 0.5, 0.5, 0.5, 0.5)))
})

# --- The unconstrained region -------------------------------------------------

test_that("a real region contains every finite point and moves none", {
  r <- real_region(2L)
  expect_true(contains(r, c(1e6, -3)))
  expect_false(contains(r, c(Inf, 0)))
  expect_equal(project(r, c(3, -1)), c(3, -1))
})

test_that("the real region's chart is the identity", {
  r <- real_region(3L)
  ch <- chart(r)
  expect_equal(ch$n_par, 3L)
  theta <- c(0.4, -2, 7)
  expect_equal(ch$to_theta(ch$from_theta(theta)), theta)
  expect_equal(ch$jacobian(theta), diag(3))
  expect_equal(dim(ch$seed(5L)), c(5L, 3L))
  u <- rbind(c(1, 2, 3), c(-1, 0, 1))
  expect_equal(ch$to_theta_batch(u), u)
})

test_that("the search finds an interior optimum on a real region", {
  # An interior optimum on an unbounded region: no constraint is active, so
  # this exercises the plain quasi-Newton behaviour of the refinement. Under
  # the unit-covariance Gaussian `E[||x - t||^2] = ||theta - t||^2 + 2`, and a
  # three-point Gauss-Hermite rule integrates the quadratic exactly.
  fam <- gaussian_family(d = 2L)
  target <- c(1.5, -0.5)
  peak <- random_variable(
    function(x) 3 - rowSums(sweep(matrix(x, ncol = 2L), 2L, target)^2),
    fam@sample_space
  )
  set.seed(1)
  found <- sup_lb(
    peak,
    null_model(fam, real_region(2L)),
    engine = gh_engine(3L),
    n_seeds = 50L,
    n_restarts = 5L
  )
  expect_equal(found@theta, target, tolerance = 1e-6)
  expect_equal(found@sup_lb, 1, tolerance = 1e-10)
})

# --- Printing -----------------------------------------------------------------

test_that("a null_model prints a summary, not a property dump", {
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  null <- null_model(
    fam,
    list(plurality_simplex(3, 2), plurality_simplex(3, 3))
  )
  out <- paste(capture.output(print(null)), collapse = "\n")
  expect_match(out, "<null_model>", fixed = TRUE)
  expect_match(out, "multinomial_family", fixed = TRUE)
  expect_match(out, "2 parts, 2 cells", fixed = TRUE)
  expect_match(out, "simplex_region: 3 vertices", fixed = TRUE)
  expect_no_match(out, "@")
  expect_match(
    format(null),
    "null_model: multinomial_family over 2 parts",
    fixed = TRUE
  )
})



test_that("ripr_init() refuses an atom outside its own part", {
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  part <- simplex_region(
    vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))
  )
  null <- null_model(fam, part)
  expect_error(
    ripr_init(fam(c(0.4, 0.35, 0.25)), null, atoms = list(c(0.6, 0.2, 0.2))),
    "lie in its own part"
  )
  expect_no_error(
    ripr_init(fam(c(0.4, 0.35, 0.25)), null, atoms = list(c(0.2, 0.5, 0.3)))
  )
})
