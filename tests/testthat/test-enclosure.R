# The enclosure interface, exercised by a method that is not Bernstein: a
# standalone Lipschitz bound on a function `g`.

lipschitz_enclosure <- S7::new_class(
  "lipschitz_enclosure",
  parent = enclosure,
  properties = list(g = S7::class_function, L = S7::class_numeric),
  constructor = function(g, L) {
    S7::new_object(
      S7::S7_object(),
      label = "lipschitz",
      g = g,
      L = L
    )
  }
)

g_at_vertices <- function(enc, box) apply(box$V, 1L, enc@g)

S7::method(enclosure_root, lipschitz_enclosure) <- function(enclosure, cell) {
  list(V = cell@vertices)
}

# Every point of a simplex is within its diameter of a vertex.
S7::method(enclosure_upper, lipschitz_enclosure) <- function(enclosure, box) {
  diameter <- max(stats::dist(box$V))
  max(g_at_vertices(enclosure, box)) + enclosure@L * diameter
}

S7::method(enclosure_best, lipschitz_enclosure) <- function(enclosure, box) {
  values <- g_at_vertices(enclosure, box)
  j <- which.max(values)
  list(value = values[[j]], theta = box$V[j, ])
}

S7::method(enclosure_split, lipschitz_enclosure) <- function(
  enclosure,
  box,
  edge
) {
  mid <- colMeans(box$V[edge, , drop = FALSE])
  lapply(edge, function(i) {
    V <- box$V
    V[i, ] <- mid
    list(V = V)
  })
}


test_that("bnb brackets the supremum with the new method", {
  # g(theta) = `-|theta - c|^2` peaks (at zero) in the interior of the
  # simplex. The gradient is bounded by `2 sqrt(2)` on the simplex.
  centre <- c(0.2, 0.3, 0.5)
  enc <- lipschitz_enclosure(
    g = function(theta) -sum((theta - centre)^2),
    L = 2 * sqrt(2)
  )
  root <- enclosure_root(enc, simplex_region(vertices = diag(3)))
  res <- branch_and_bound(enc, root, tol = 1e-2, max_iter = 20000L)

  expect_true(res$converged)
  expect_gte(res$bound, 0)
  expect_lte(res$bound - res$incumbent, 1e-2)
  expect_lte(res$incumbent, 0)
  expect_gt(res$iterations, 0L)

  # The leaves tile the original simplex:
  area <- function(V) abs(det(cbind(V[, 1:2], 1))) / 2
  leaves <- Filter(function(b) b$fate != "split", res$history)
  expect_equal(
    sum(vapply(leaves, function(b) area(b$V), numeric(1L))),
    area(diag(3)),
    tolerance = 1e-12
  )
})

test_that("the default edge is the longest", {
  enc <- lipschitz_enclosure(g = function(theta) 0, L = 0)
  V <- diag(3)
  V[1L, ] <- c(0.5, 0.5, 0)
  expect_identical(enclosure_edge(enc, list(V = V)), c(2L, 3L))
})

test_that("bounding_method() errors if the pair is not implemented", {
  fam <- gaussian_family(d = 2L)
  cell <- simplex_region(vertices = rbind(c(0, 0), c(1, 0), c(0, 1)))
  why <- bounding_method(fam, cell)
  expect_type(why, "character")
  expect_match(why, "No bounding method is implemented for gaussian_family")
})

test_that("bounding_method() selects bernstein for multinomial", {
  # Points included: a point is a one-vertex simplex, bounded exactly.
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  expect_true(S7::S7_inherits(
    bounding_method(fam, simplex_region(vertices = diag(3))),
    bernstein_enclosure
  ))
  expect_true(S7::S7_inherits(
    bounding_method(fam, point_region(theta = c(0.5, 0.3, 0.2))),
    bernstein_enclosure
  ))
})

test_that("a family with no enclosure of its own gets point evaluation", {
  # Enumerable sample space and no method registered: only points can be
  # bounded, by exact evaluation.
  fam <- parametric_family(
    sample_space = count_space(n_trials = 4L, k = 3L),
    parameter_space = simplex_region(vertices = diag(3))
  )
  expect_true(S7::S7_inherits(
    bounding_method(fam, point_region(theta = c(0.5, 0.3, 0.2))),
    point_enclosure
  ))
  expect_match(
    bounding_method(fam, simplex_region(vertices = diag(3))),
    "No bounding method is implemented"
  )
})

test_that("a point enclosure bounds a point by its exact value", {
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  outcomes <- enumerate_space(fam@sample_space)
  x <- likelihood(fam(c(0.5, 0.3, 0.2))) / likelihood(fam(c(0.2, 0.5, 0.3)))
  theta <- c(0.4, 0.35, 0.25)
  enc <- prepare_enclosure(point_enclosure(), x, fam, list())
  root <- enclosure_root(enc, point_region(theta = theta))
  direct <- sum(exp(log_density(fam(theta), outcomes)) * x(outcomes))
  expect_equal(enclosure_upper(enc, root), direct)
  expect_equal(enclosure_best(enc, root)$value, direct)
  expect_error(enclosure_split(enc, root, c(1L, 1L)), "cannot be split")
})
