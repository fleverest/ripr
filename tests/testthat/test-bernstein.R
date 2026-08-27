# Properties of R/bernstein.R.
#
# Almost everything here is an indexing contract. The de Casteljau routines are
# gathers against precomputed index maps, so a mis-indexed coefficient produces
# a plausible number rather than an error, and a plausible number that is too
# small is an invalid certificate. Hence the emphasis: the arithmetic is a
# handful of convex combinations and is not what goes wrong.
#
# The reference throughout is direct evaluation of the Bernstein form,
#
#   G(lambda) = sum_alpha c_alpha * choose(n; alpha) * prod_j lambda_j^alpha_j,
#
# written out below rather than reused from the package, so that the tests
# check `dc_*()` against the definition instead of against themselves.

# `choose(n; alpha) prod lambda^alpha` for every row of `tally`.
bern_basis <- function(tally, n, lambda) {
  coef <- exp(lgamma(n + 1) - rowSums(lgamma(tally + 1)))
  coef * apply(tally, 1L, function(a) prod(lambda^a))
}

bern_eval <- function(coef, lat, lambda) {
  sum(coef * bern_basis(lat$tally, lat$n, lambda))
}

# Uniformly random barycentric weights on the K-simplex.
rand_lambda <- function(K) {
  g <- stats::rgamma(K, shape = 1)
  g / sum(g)
}

# --- Coefficient ordering -----------------------------------------------------

test_that("compositions() agrees row-for-row with enumerate_counts()", {
  # Load-bearing. `certify()` passes `x(enumerate_space(space))` straight in
  # as Bernstein coefficients, so the count-space order and the lattice
  # row order must be the same enumeration. They are computed by two unrelated
  # routines -- recursive prefixing here, stars and bars there -- and nothing
  # else forces them to agree.
  for (k in 2:5) {
    for (n in 1:6) {
      expect_identical(compositions(n, k), enumerate_counts(n, k))
    }
  }
})

test_that("compositions() is lexicographic, integer, and the right size", {
  tally <- compositions(4L, 3L)
  expect_type(tally, "integer")
  expect_null(dimnames(tally))
  expect_equal(nrow(tally), bernstein_size(4L, 3L))
  expect_true(all(rowSums(tally) == 4L))
  expect_false(anyDuplicated(apply(tally, 1L, paste, collapse = ",")) > 0L)

  keys <- apply(tally, 1L, paste, collapse = ",")
  expect_identical(keys, sort(keys, method = "radix"))
})

test_that("bernstein_size() counts the multinomial support", {
  for (k in 2:5) {
    for (n in 1:6) {
      # `choose()` returns a double, deliberately: the count overflows integer
      # long before the lattice becomes buildable.
      expect_equal(bernstein_size(n, k), nrow(enumerate_counts(n, k)))
    }
  }
})

# --- Size guard ---------------------------------------------------------------

test_that("check_bernstein_size() admits and refuses on the stated boundary", {
  size <- bernstein_size(20L, 5L)
  expect_silent(check_bernstein_size(20L, 5L, size))
  expect_error(
    check_bernstein_size(20L, 5L, size - 1L),
    "above `max_coefficients`"
  )
})

test_that("largest_batch() is the exact threshold, not an estimate", {
  # The refusal quotes this number back to the caller, so an off-by-one here is
  # advice to try something that will also be refused.
  for (k in 2:5) {
    for (budget in c(10L, 100L, 5000L)) {
      n <- largest_batch(k, budget)
      expect_lte(bernstein_size(n, k), budget)
      expect_gt(bernstein_size(n + 1L, k), budget)
    }
  }
})

test_that("largest_batch() returns 0 when even one trial will not fit", {
  expect_identical(largest_batch(10L, 2L), 0L)
})

# --- Lattice ------------------------------------------------------------------

test_that("bernstein_lattice() locates the vertex coefficients", {
  # PBP 10.2: b(a_j) = b_{n e_j}. `vertex_values()` reads these as exact values
  # of the polynomial and `box_best()` uses them as the incumbent, so a wrong
  # position here inflates the reported lower bound.
  for (K in 2:5) {
    lat <- bernstein_lattice(6L, K)
    expect_length(lat$vertex, K)
    for (j in seq_len(K)) {
      expect_identical(lat$tally[lat$vertex[j], ], replace(integer(K), j, 6L))
    }
  }
})

test_that("bernstein_lattice() enumerates every edge once", {
  lat <- bernstein_lattice(3L, 4L)
  expect_identical(dim(lat$edges), c(2L, 6L))
  expect_true(all(lat$edges[1L, ] < lat$edges[2L, ]))
  expect_identical(anyDuplicated(t(lat$edges)), 0L)
})

test_that("the degree ladder `up` steps by one basis vector", {
  lat <- bernstein_lattice(4L, 3L)
  for (m in seq_len(lat$n)) {
    # `up[[m + 1]][beta, i]` is the position of `beta + e_i` among the
    # degree-`m` multi-indices, indexed by the degree-`m - 1` row `beta`.
    lo <- compositions(m - 1L, lat$K)
    hi <- compositions(m, lat$K)
    for (i in seq_len(lat$K)) {
      shift <- rep(replace(integer(lat$K), i, 1L), each = nrow(lo))
      expect_identical(hi[lat$up[[m + 1L]][, i], , drop = FALSE], lo + shift)
    }
  }
})

# --- Hand-worked examples -----------------------------------------------------

# Some small hand-worked examples I worked through manually.

# The two helpers translate between triangles and the corresponding coefficient
# vectors so that no assertion here depends on the row order. Row orders are
# asserted separately above.
#
# Layout throughout, for degree 2 on K = 3: vertex 1 at the apex, vertex 2 at
# the lower left, vertex 3 at the lower right.
#
#              b_200
#          b_110   b_101
#      b_020   b_011   b_002
# The coefficient vector would be (b_002, b_011, b_020, b_101, b_110, b_200),
# i.e. the order of `compositions(2, 3)`.

# Build a coefficient vector from multi-indices written as strings, e.g.
# `coef_by_index(lat, "200" = 0, "110" = 6, ...)`.
coef_by_index <- function(lat, ...) {
  entries <- c(...)
  keys <- apply(lat$tally, 1L, paste, collapse = "")
  pos <- match(names(entries), keys)
  stopifnot(!anyNA(pos), length(entries) == lat$n_coef)
  out <- numeric(lat$n_coef)
  out[pos] <- entries
  out
}

# Read one entry of a pyramid by its multi-index. Level `l` of the pyramid
# holds degree `n - l`, so the degree of the index picks the level out.
pyr_at <- function(pyr, lat, index) {
  degree <- sum(as.integer(strsplit(index, "")[[1L]]))
  keys <- apply(compositions(degree, lat$K), 1L, paste, collapse = "")
  pyr[[lat$n - degree + 1L]][[match(index, keys)]]
}

test_that("the univariate pyramid matches the worked example", {
  #   b_20   b_11   b_02          2      3      1
  #      b_10   b_01        ->       2.5    2
  #         b_00                        2.25
  #
  # Cubic example: control points (2, 3, 1) at t = 1/2.
  # `compositions(2, 2)` runs (0,2), (1,1), (2,0), and row (a1, a2) carries
  # C(n; a) l1^a1 l2^a2, so with l = (1 - t, t) the standard index i is a2.
  lat <- bernstein_lattice(2L, 2L)
  coef <- coef_by_index(lat, "20" = 2, "11" = 3, "02" = 1)
  expect_identical(coef, c(1, 3, 2))

  pyr <- dc_pyramid(coef, lat, c(0.5, 0.5))
  expect_equal(pyr_at(pyr, lat, "10"), 2.5)
  expect_equal(pyr_at(pyr, lat, "01"), 2)
  # b(1/2) = 2(1/4) + 3(2)(1/4) + 1(1/4) = 2.25.
  expect_equal(pyr_at(pyr, lat, "00"), 2.25)
})

test_that("the univariate children are the two halves of the worked example", {
  # The two edges of the pyramid: (2, 2.5, 2.25) for t in [0, 1/2] and
  # (2.25, 2, 1) for t in [1/2, 1], sharing the split value 2.25.
  lat <- bernstein_lattice(2L, 2L)
  coef <- coef_by_index(lat, "20" = 2, "11" = 3, "02" = 1)
  pyr <- dc_pyramid(coef, lat, c(0.5, 0.5))

  expect_equal(
    dc_child(pyr, lat, 1L),
    coef_by_index(lat, "20" = 2.25, "11" = 2, "02" = 1)
  )
  expect_equal(
    dc_child(pyr, lat, 2L),
    coef_by_index(lat, "20" = 2, "11" = 2.5, "02" = 2.25)
  )

  # Through `bisect()`, which fixes which half comes back first: vertex 1 is
  # t = 0, so replacing it leaves the RIGHT half, t in [1/2, 1].
  kids <- bisect(list(V = diag(2L), coef = coef), 1L, 2L, lat)
  expect_equal(kids[[1L]]$V, matrix(c(0.5, 0.5, 0, 1), nrow = 2L))
  expect_equal(
    kids[[1L]]$coef,
    coef_by_index(lat, "20" = 2.25, "11" = 2, "02" = 1)
  )
  expect_equal(kids[[2L]]$V, matrix(c(1, 0, 0.5, 0.5), nrow = 2L))
  expect_equal(
    kids[[2L]]$coef,
    coef_by_index(lat, "20" = 2, "11" = 2.5, "02" = 2.25)
  )
})

test_that("the two-dimensional pyramid is a stack of shrinking triangles", {
  # Degree 2 on K = 3 at lambda = (1/2, 1/4, 1/4).
  #
  #        level 0              level 1        level 2
  #
  #           0
  #        6     2                 2
  #     8     0     4           5     2          2.75
  #
  # Each entry of the next triangle sits inside one upward-pointing triangle
  # of the current one, and is the lambda-weighted average of its three
  # corners: lambda_1 on the top corner, lambda_2 on the lower left,
  # lambda_3 on the lower right. That is the whole of de Casteljau.
  #
  #   b_100 = (1/2)(0) + (1/4)(6) + (1/4)(2) = 2
  #   b_010 = (1/2)(6) + (1/4)(8) + (1/4)(0) = 5
  #   b_001 = (1/2)(2) + (1/4)(0) + (1/4)(4) = 2
  #   b_000 = (1/2)(2) + (1/4)(5) + (1/4)(2) = 2.75
  lat <- bernstein_lattice(2L, 3L)
  coef <- coef_by_index(
    lat,
    "200" = 0,
    "110" = 6,
    "101" = 2,
    "020" = 8,
    "011" = 0,
    "002" = 4
  )
  pyr <- dc_pyramid(coef, lat, c(0.5, 0.25, 0.25))

  expect_equal(pyr_at(pyr, lat, "100"), 2)
  expect_equal(pyr_at(pyr, lat, "010"), 5)
  expect_equal(pyr_at(pyr, lat, "001"), 2)
  expect_equal(pyr_at(pyr, lat, "000"), 2.75)

  # Against the closed form. For this triangle the l2 l3 and l1^2 terms are
  # zero, leaving 4 l3^2 + 8 l2^2 + 4 l1 l3 + 12 l1 l2.
  closed <- function(l) {
    4 * l[3]^2 + 8 * l[2]^2 + 4 * l[1] * l[3] + 12 * l[1] * l[2]
  }
  expect_equal(closed(c(0.5, 0.25, 0.25)), 0.25 + 0.5 + 0.5 + 1.5)
  expect_equal(pyr_at(pyr, lat, "000"), closed(c(0.5, 0.25, 0.25)))
})

test_that("each two-dimensional child is a triangle read off one pyramid edge", {
  # V^[i] replaces vertex i by the split point. Its triangle keeps the edge
  # opposite vertex i intact, puts the apex value at vertex i, and fills the
  # rows between from the intermediate levels: the row at distance k from the
  # opposite edge comes from level k.
  #
  #    parent         level 1        level 0
  #
  #      0
  #   6     2            2
  # 8    0    4       5     2          2.75
  #
  #
  #      V^[1]           V^[2]           V^[3]
  #
  #      2.75              0               0
  #    5     2          2     2         6     2
  #  8    0    4    2.75   2     4    8     5   2.75
  #
  # Reading the child triangles: V^[1] has apex 2.75 and its base row 8 0 4 is
  # the parent's, since replacing vertex 1 does not move the edge opposite it.
  # Its middle row 5 2 is level 1. V^[2] and V^[3] are the same construction
  # rotated to the other two corners.
  lat <- bernstein_lattice(2L, 3L)
  coef <- coef_by_index(
    lat,
    "200" = 0,
    "110" = 6,
    "101" = 2,
    "020" = 8,
    "011" = 0,
    "002" = 4
  )
  pyr <- dc_pyramid(coef, lat, c(0.5, 0.25, 0.25))

  expect_equal(
    dc_child(pyr, lat, 1L),
    coef_by_index(
      lat,
      "200" = 2.75,
      "110" = 5,
      "101" = 2,
      "020" = 8,
      "011" = 0,
      "002" = 4
    )
  )
  expect_equal(
    dc_child(pyr, lat, 2L),
    coef_by_index(
      lat,
      "200" = 0,
      "110" = 2,
      "101" = 2,
      "020" = 2.75,
      "011" = 2,
      "002" = 4
    )
  )
  expect_equal(
    dc_child(pyr, lat, 3L),
    coef_by_index(
      lat,
      "200" = 0,
      "110" = 6,
      "101" = 2,
      "020" = 8,
      "011" = 5,
      "002" = 2.75
    )
  )

  # Stated as properties rather than as transcribed numbers, so a change to the
  # example does not quietly weaken the test.
  for (i in seq_len(3L)) {
    child <- dc_child(pyr, lat, i)
    # The apex value lands at the vertex that was replaced.
    expect_equal(child[[lat$vertex[i]]], 2.75)
    # The opposite face -- the entries with alpha_i = 0 -- is untouched.
    opposite <- lat$tally[, i] == 0L
    expect_equal(child[opposite], coef[opposite])
  }
})

test_that("the Bernstein form reproduces linear functions exactly", {
  # Coefficients equal to the distance from the edge opposite vertex 1:
  #
  #           2
  #        1     1
  #     0     0     0
  #
  # give sum_alpha alpha_1 C(n; alpha) l^alpha = n l_1, with no approximation
  # at any lambda. It is a fact about the basis, not about any family: there is
  # no probability anywhere in this block.
  #
  # Read alongside the pmf/basis test in `test-certify.R`, this is the half of
  # "E_theta[x_1] = n theta_1" that lives here. The other half -- that the
  # coefficients ARE X on the lattice -- is a statement about
  # `enumerate_space()` and cannot be checked without a family.
  n <- 2L
  lat <- bernstein_lattice(n, 3L)
  coef <- coef_by_index(
    lat,
    "200" = 2,
    "110" = 1,
    "101" = 1,
    "020" = 0,
    "011" = 0,
    "002" = 0
  )

  for (lambda in list(
    c(0.5, 0.25, 0.25),
    c(1, 0, 0),
    rep(1 / 3, 3),
    c(0.5, 0.5, 0)
  )) {
    expect_equal(
      pyr_at(dc_pyramid(coef, lat, lambda), lat, "000"),
      n * lambda[1L]
    )
  }

  # It is not special to this degree or this many vertices.
  for (m in 2:5) {
    for (K in 2:4) {
      other <- bernstein_lattice(as.integer(m), as.integer(K))
      lambda <- c(0.4, 0.3, 0.2, 0.1)[seq_len(K)]
      lambda <- lambda / sum(lambda)
      expect_equal(
        dc_pyramid(as.numeric(other$tally[, 1L]), other, lambda)[[m + 1L]],
        m * lambda[1L]
      )
    }
  }

  # And it survives subdivision: a child of a linear function is that same
  # linear function, read in the child's coordinates.
  kids <- bisect(list(V = diag(3L), coef = coef), 1L, 2L, lat)
  for (kid in kids) {
    for (mu in list(c(0.5, 0.25, 0.25), rep(1 / 3, 3), c(0, 1, 0))) {
      expect_equal(
        pyr_at(dc_pyramid(kid$coef, lat, mu), lat, "000"),
        n * as.vector(kid$V %*% mu)[1L]
      )
    }
  }
})


# --- de Casteljau -------------------------------------------------------------

test_that("the pyramid apex is the polynomial's value", {
  # PBP 10.4. This is the one place the arithmetic itself is checked; if the
  # apex is right then `dc_step()` and the `up` map are both right.
  set.seed(11)
  for (K in 2:4) {
    lat <- bernstein_lattice(5L, K)
    coef <- stats::rnorm(lat$n_coef)
    for (rep in 1:5) {
      lambda <- rand_lambda(K)
      apex <- dc_pyramid(coef, lat, lambda)[[lat$n + 1L]]
      expect_length(apex, 1L)
      expect_equal(apex, bern_eval(coef, lat, lambda))
    }
  }
})

test_that("the pyramid at a vertex returns that vertex's coefficient", {
  set.seed(12)
  lat <- bernstein_lattice(5L, 3L)
  coef <- stats::rnorm(lat$n_coef)
  for (j in seq_len(lat$K)) {
    lambda <- replace(numeric(lat$K), j, 1)
    expect_equal(
      dc_pyramid(coef, lat, lambda)[[lat$n + 1L]],
      coef[lat$vertex[j]]
    )
  }
})

test_that("levels of the pyramid have the degrees they claim to", {
  lat <- bernstein_lattice(4L, 3L)
  pyr <- dc_pyramid(stats::rnorm(lat$n_coef), lat, rand_lambda(3L))
  expect_length(pyr, lat$n + 1L)
  for (l in 0:lat$n) {
    expect_length(pyr[[l + 1L]], bernstein_size(lat$n - l, lat$K))
  }
})

# --- Subdivision --------------------------------------------------------------

test_that("a child's coefficients describe the same polynomial", {
  # Leroy Algorithm 2.13. The child is a reparametrisation, not an
  # approximation: evaluating the child's form at barycentric `mu` must equal
  # evaluating the parent's at the corresponding point of the parent simplex.
  # This is the test that catches a transposed or mis-levelled read-off.
  set.seed(13)
  for (K in 3:4) {
    lat <- bernstein_lattice(4L, K)
    coef <- stats::rnorm(lat$n_coef)
    box <- list(V = diag(K), coef = coef)
    lambda <- rand_lambda(K)
    kids <- subdivide(box, lat, lambda)
    expect_length(kids, K)
    for (kid in kids) {
      for (rep in 1:3) {
        mu <- rand_lambda(K)
        expect_equal(
          bern_eval(kid$coef, lat, mu),
          bern_eval(coef, lat, as.vector(kid$V %*% mu))
        )
      }
    }
  }
})

test_that("subdivide() drops the degenerate children", {
  # lambda_i = 0 puts the new point in the face opposite vertex i, so V^[i]
  # would be flat and its expansion meaningless.
  lat <- bernstein_lattice(3L, 4L)
  box <- list(V = diag(4L), coef = stats::rnorm(lat$n_coef))
  lambda <- c(0.5, 0.5, 0, 0)
  kids <- subdivide(box, lat, lambda)
  expect_length(kids, 2L)
})

test_that("bisect() returns the child with vertex `p` REPLACED first", {
  # Pinning the labelling, which is the opposite of "keeps vertex p".
  # Assuming the wrong order might cause us to certify over the wrong half.
  lat <- bernstein_lattice(3L, 3L)
  box <- list(V = diag(3L), coef = stats::rnorm(lat$n_coef))
  mid <- c(0.5, 0.5, 0)
  kids <- bisect(box, 1L, 2L, lat)

  expect_equal(kids[[1L]]$V[, 1L], mid)
  expect_equal(kids[[1L]]$V[, 2L], c(0, 1, 0))
  expect_equal(kids[[2L]]$V[, 2L], mid)
  expect_equal(kids[[2L]]$V[, 1L], c(1, 0, 0))
})

test_that("the two halves of a bisection cover the parent", {
  set.seed(14)
  lat <- bernstein_lattice(4L, 3L)
  coef <- stats::rnorm(lat$n_coef)
  kids <- bisect(list(V = diag(3L), coef = coef), 2L, 3L, lat)
  # Every point of the parent is in one half or the other, so the two child
  # ranges must between them cover the parent's range.
  parent_max <- max(vapply(
    seq_len(400L),
    function(i) bern_eval(coef, lat, rand_lambda(3L)),
    numeric(1L)
  ))
  expect_gte(max(vapply(kids, box_bound, numeric(1L))), parent_max)
})

# --- Enclosure ----------------------------------------------------------------

test_that("box_bound() encloses the polynomial from above", {
  # Garloff (1986) via the convex hull property, PBP 10.2 with 10.3 Remark 2.
  # Sampled rather than proved, so this is a smoke test for the sign of the
  # inequality rather than a proof of the enclosure.
  set.seed(15)
  lat <- bernstein_lattice(6L, 3L)
  coef <- stats::rnorm(lat$n_coef)
  box <- list(V = diag(3L), coef = coef)
  vals <- vapply(
    seq_len(2000L),
    function(i) bern_eval(coef, lat, rand_lambda(3L)),
    numeric(1L)
  )
  expect_gte(box_bound(box), max(vals))
  expect_lte(min(coef), min(vals))
})

test_that("vertex values are exact, so box_best() is attained", {
  set.seed(16)
  lat <- bernstein_lattice(5L, 4L)
  coef <- stats::rnorm(lat$n_coef)
  V <- diag(4L)
  V[, 1L] <- c(0.5, 0.5, 0, 0)
  box <- list(V = V, coef = reparametrise_to(coef, lat, V))

  best <- box_best(box, lat)
  # `theta` is in the *original* simplex's coordinates, so evaluate the
  # original form there.
  expect_equal(best$value, bern_eval(coef, lat, best$theta))
  expect_lte(best$value, box_bound(box))
})

test_that("boxes_best() picks the best box", {
  lat <- bernstein_lattice(2L, 3L)
  a <- list(V = diag(3L), coef = rep(1, lat$n_coef))
  b <- list(V = diag(3L), coef = rep(2, lat$n_coef))
  expect_equal(boxes_best(list(a, b), lat)$value, 2)
  expect_equal(boxes_best(list(b, a), lat)$value, 2)
})

test_that("longest_edge() returns the longest edge", {
  lat <- bernstein_lattice(2L, 3L)
  V <- diag(3L)
  V[, 1L] <- c(0.5, 0.5, 0) # shortens edges 1-2 and 1-3, leaving 2-3 longest
  expect_identical(longest_edge(V, lat$edges), c(2L, 3L))
})

# --- Reparametrisation --------------------------------------------------------

test_that("reparametrise_to() onto the standard simplex is the identity", {
  set.seed(17)
  lat <- bernstein_lattice(4L, 3L)
  coef <- stats::rnorm(lat$n_coef)
  expect_identical(reparametrise_to(coef, lat, diag(3L)), coef)
})

test_that("reparametrise_to() preserves the polynomial", {
  set.seed(18)
  lat <- bernstein_lattice(4L, 4L)
  coef <- stats::rnorm(lat$n_coef)
  V <- diag(4L)
  V[, 1L] <- c(0.5, 0.5, 0, 0)
  out <- reparametrise_to(coef, lat, V)
  for (rep in 1:5) {
    mu <- rand_lambda(4L)
    expect_equal(
      bern_eval(out, lat, mu),
      bern_eval(coef, lat, as.vector(V %*% mu))
    )
  }
})

test_that("reparametrise_to() agrees with subdivide() on a single replacement", {
  set.seed(19)
  lat <- bernstein_lattice(4L, 3L)
  coef <- stats::rnorm(lat$n_coef)
  lambda <- c(0.25, 0.5, 0.25)
  V <- diag(3L)
  V[, 2L] <- lambda
  kid <- subdivide(list(V = diag(3L), coef = coef), lat, lambda)[[2L]]
  expect_equal(reparametrise_to(coef, lat, V), kid$coef)
})

test_that("reparametrise_to() does not depend on the order of the vertices", {
  # A facet is a set. The plurality facet {theta_1 <= theta_2} spanned by
  # (tie, e_2, e_3) is the same simplex as (e_2, e_3, tie), and both must
  # describe the same polynomial. This needs no special handling: the polar
  # form is symmetric (PBP 11.2), so permuting the columns just permutes the
  # coefficients correspondingly.
  set.seed(28)
  lat <- bernstein_lattice(4L, 3L)
  coef <- stats::rnorm(lat$n_coef)
  tie <- c(0.5, 0.5, 0)
  matched <- cbind(tie, c(0, 1, 0), c(0, 0, 1))
  permuted <- cbind(c(0, 1, 0), c(0, 0, 1), tie)

  for (V in list(matched, permuted)) {
    out <- reparametrise_to(coef, lat, V)
    for (rep in 1:4) {
      mu <- rand_lambda(3L)
      expect_equal(
        bern_eval(out, lat, mu),
        bern_eval(coef, lat, as.vector(V %*% mu))
      )
    }
  }

  # The two are the same coefficients up to the relabelling of the vertices.
  expect_setequal(
    round(reparametrise_to(coef, lat, matched), 10),
    round(reparametrise_to(coef, lat, permuted), 10)
  )
})

test_that("reparametrise_to() reaches sub-simplices with no vertex in common", {
  # The case PBP Figure 11.4 rules out for repeated subdivision: the medial
  # triangle, every vertex interior to an edge of the original. No ordering of
  # single-vertex replacements reaches it, so this is the test that separates
  # the blossom construction (PBP 11.2) from the one it replaced.
  set.seed(29)
  lat <- bernstein_lattice(4L, 3L)
  coef <- stats::rnorm(lat$n_coef)
  medial <- cbind(c(0.5, 0.5, 0), c(0.5, 0, 0.5), c(0, 0.5, 0.5))

  out <- reparametrise_to(coef, lat, medial)
  for (rep in 1:8) {
    mu <- rand_lambda(3L)
    expect_equal(
      bern_eval(out, lat, mu),
      bern_eval(coef, lat, as.vector(medial %*% mu))
    )
  }

  # Exactly, not approximately: every step is a convex combination, so there is
  # no drift to absorb. A tolerance here would hide the difference between this
  # construction and the extrapolating one.
  expect_lt(
    max(abs(vapply(
      seq_len(50L),
      function(i) {
        mu <- rand_lambda(3L)
        bern_eval(out, lat, mu) - bern_eval(coef, lat, as.vector(medial %*% mu))
      },
      numeric(1L)
    ))),
    1e-13
  )
})

test_that("reparametrise_to() preserves the enclosure on an interior sub-simplex", {
  # Load-bearing for `certify()`: the bound is min/max of the coefficients, so
  # reparametrisation is only sound if the new coefficients still bracket the
  # polynomial over the new simplex.
  set.seed(30)
  for (K in 2:4) {
    lat <- bernstein_lattice(4L, as.integer(K))
    coef <- stats::rnorm(lat$n_coef)
    V <- vapply(
      seq_len(K),
      function(i) {
        z <- stats::rgamma(K, shape = 2)
        z / sum(z)
      },
      numeric(K)
    )
    out <- reparametrise_to(coef, lat, V)
    for (rep in 1:40) {
      mu <- rand_lambda(K)
      value <- bern_eval(coef, lat, as.vector(V %*% mu))
      expect_lte(value, max(out) + 1e-12)
      expect_gte(value, min(out) - 1e-12)
    }
  }
})

test_that("reparametrise_to() refuses vertices outside the simplex or degenerate", {
  # Outside means `dc_step()` extrapolates rather than interpolating, which is
  # what costs the stability and the enclosure both.
  lat <- bernstein_lattice(3L, 3L)
  coef <- rep(1, lat$n_coef)
  outside <- cbind(c(1.5, -0.5, 0), c(0, 1, 0), c(0, 0, 1))
  degenerate <- cbind(c(0.5, 0.5, 0), c(0.5, 0.5, 0), c(0, 0, 1))

  expect_error(reparametrise_to(coef, lat, outside), "standard simplex")
  expect_error(reparametrise_to(coef, lat, degenerate), "non-degenerate")
})

# --- Branch and bound ---------------------------------------------------------

# A fixed, reproducible seed box on the standard simplex.
seed_box <- function(lat, coef) list(V = diag(lat$K), coef = coef)

test_that("certify_sup() brackets the true supremum", {
  set.seed(20)
  lat <- bernstein_lattice(6L, 3L)
  coef <- stats::rnorm(lat$n_coef, sd = 2)
  res <- certify_sup(
    list(seed_box(lat, coef)),
    lat,
    tol = 1e-9,
    max_iter = 2000L
  )

  vals <- vapply(
    seq_len(5000L),
    function(i) bern_eval(coef, lat, rand_lambda(3L)),
    numeric(1L)
  )
  expect_gte(res$bound, max(vals))
  expect_gte(res$bound, res$incumbent)
  expect_equal(res$incumbent, bern_eval(coef, lat, res$theta))
})

test_that("the bound is valid at every iteration, not just at convergence", {
  # The whole claim of the method: refinement buys tightness, not validity. If
  # this fails then a run that hits `max_nodes` returns an invalid certificate.
  set.seed(21)
  lat <- bernstein_lattice(6L, 3L)
  coef <- stats::rnorm(lat$n_coef, sd = 2)
  truth <- certify_sup(
    list(seed_box(lat, coef)),
    lat,
    tol = 1e-12,
    max_iter = 5000L
  )$incumbent

  bounds <- vapply(
    c(0L, 1L, 2L, 5L, 20L, 100L),
    function(m) {
      certify_sup(list(seed_box(lat, coef)), lat, tol = 0, max_iter = m)$bound
    },
    numeric(1L)
  )
  expect_true(all(bounds >= truth))
  expect_false(is.unsorted(rev(bounds))) # non-increasing in the budget
})

test_that("the bound is non-increasing as the tolerance tightens", {
  set.seed(22)
  lat <- bernstein_lattice(20L, 3L)
  coef <- stats::rnorm(lat$n_coef, sd = 2)
  bounds <- vapply(
    c(1, 1e-1, 1e-3, 1e-6),
    function(tol) {
      certify_sup(
        list(seed_box(lat, coef)),
        lat,
        tol = tol,
        max_iter = 5000L
      )$bound
    },
    numeric(1L)
  )
  expect_false(is.unsorted(rev(bounds)))
})

test_that("certify_sup() stops within tol of the incumbent", {
  set.seed(23)
  lat <- bernstein_lattice(5L, 3L)
  coef <- stats::rnorm(lat$n_coef, sd = 2)
  tol <- 1e-4
  res <- certify_sup(
    list(seed_box(lat, coef)),
    lat,
    tol = tol,
    max_iter = 10000L
  )
  expect_lt(res$iterations, 10000L)
  # `eta` is added to the reported bound after the stopping test.
  expect_lte(res$bound - res$incumbent, tol + 1e-12)
})

test_that("the rounding slack is added and is non-negative", {
  set.seed(24)
  lat <- bernstein_lattice(5L, 3L)
  coef <- stats::rnorm(lat$n_coef, sd = 2)
  args <- list(list(seed_box(lat, coef)), lat, tol = 1e-9, max_iter = 2000L)
  with_eta <- do.call(certify_sup, c(args, list(round_slack = TRUE)))
  without <- do.call(certify_sup, c(args, list(round_slack = FALSE)))
  expect_gte(with_eta$bound, without$bound)
  expect_equal(with_eta$incumbent, without$incumbent)
})

test_that("pruning nodes via slack doesn't hurt the bound", {
  set.seed(25)
  lat <- bernstein_lattice(6L, 3L)
  coef <- stats::rnorm(lat$n_coef, sd = 2)
  truth <- certify_sup(
    list(seed_box(lat, coef)),
    lat,
    tol = 1e-12,
    max_iter = 5000L
  )$incumbent
  res <- certify_sup(
    list(seed_box(lat, coef)),
    lat,
    tol = 1e-9,
    max_iter = 2000L,
    slack = 0.5
  )
  expect_gte(res$bound, truth)
  expect_lte(res$bound, res$incumbent + 0.5 + 1e-9)
})

test_that("keep_argmax retains an enclosure of the maximiser", {
  set.seed(26)
  lat <- bernstein_lattice(6L, 3L)
  coef <- stats::rnorm(lat$n_coef, sd = 2)
  res <- certify_sup(
    list(seed_box(lat, coef)),
    lat,
    tol = 1e-6,
    max_iter = 2000L,
    keep_argmax = TRUE
  )
  expect_gt(length(res$active), 0L)
  # The incumbent's argmax must lie in some retained box, so some retained box
  # must bound it from above.
  expect_gte(
    max(vapply(res$active, box_bound, numeric(1L))),
    res$incumbent - 1e-9
  )
})

test_that("keep_argmax and slack are mutually exclusive", {
  lat <- bernstein_lattice(3L, 3L)
  expect_error(
    certify_sup(
      list(seed_box(lat, rep(1, lat$n_coef))),
      lat,
      slack = 0.1,
      keep_argmax = TRUE
    ),
    "`slack` must be 0"
  )
})

test_that("a constant polynomial is certified without any subdivision", {
  lat <- bernstein_lattice(4L, 3L)
  res <- certify_sup(list(seed_box(lat, rep(2.5, lat$n_coef))), lat, tol = 1e-9)
  expect_identical(res$iterations, 0L)
  expect_equal(res$incumbent, 2.5)
  expect_lt(res$bound - 2.5, 1e-12)
})

test_that("several seeds are certified as their union", {
  set.seed(27)
  lat <- bernstein_lattice(5L, 3L)
  coef <- stats::rnorm(lat$n_coef, sd = 2)
  halves <- bisect(seed_box(lat, coef), 1L, 2L, lat)
  joint <- certify_sup(halves, lat, tol = 1e-9, max_iter = 4000L)
  separate <- vapply(
    halves,
    function(b) certify_sup(list(b), lat, tol = 1e-9, max_iter = 4000L)$bound,
    numeric(1L)
  )
  expect_lte(joint$bound, max(separate) + 1e-9)
  expect_gte(joint$bound, joint$incumbent)
})

test_that("a run that prunes all nodes converges rather than running out", {
  # The other stopping path. Every other test converges by reaching tolerance.
  # Without this the empty-active-set branch never executes under test, and a
  # run that pruned its way to an answer would be reported as budget-limited,
  # which says the bound is loose when it is as tight as the method gets.
  set.seed(47)
  lat <- bernstein_lattice(6L, 3L)
  # Picking coefficients that are not maximised on the vertices:
  coef <- stats::rnorm(lat$n_coef, sd = 2)
  while (max(coef) == boxes_best(list(seed_box(lat, coef)), lat)$value) {
    coef[which.max(coef)] <- min(coef) - 1
  }
  seeds <- list(seed_box(lat, coef))
  res <- certify_sup(seeds, lat, tol = 0, max_iter = 500L, slack = 0.5)

  expect_length(res$active, 0L)
  expect_true(res$converged)
  expect_false(res$budget_hit)
  expect_lt(res$iterations, 500L)
})

test_that("pruning discards nodes that only tie the incumbent", {
  # `>` not `>=`: a node bounded by exactly the incumbent cannot improve on it.
  # Keeping it is conservative: same bound, more work, so no test that reads
  # the bound can tell the difference. `keep_argmax` deliberately inverts this,
  # which is the case that must not be broken while fixing the other.
  lat <- bernstein_lattice(2L, 3L)
  flat <- function(u) {
    node(list(V = diag(3L), coef = rep(u, lat$n_coef)), id = 1L)
  }
  nodes <- list(flat(1), flat(5), flat(9))

  expect_length(prune_active(nodes, 5, 0, 0, keep_argmax = FALSE)$keep, 1L)
  expect_length(prune_active(nodes, 5, 0, 0, keep_argmax = TRUE)$keep, 2L)
})

# --- Regression ---------------------------------------------------------------

test_that("certify_sup() returns the bound it has always returned", {
  # A pinned value on a fixed input. Everything else here tests properties,
  # which a systematically shifted bound would still satisfy.
  set.seed(55)
  lat <- bernstein_lattice(10L, 3L)
  coef <- abs(stats::rnorm(lat$n_coef, sd = 2))^2
  res <- certify_sup(
    list(list(V = diag(3L), coef = coef)),
    lat,
    tol = 1e-9,
    max_iter = 1000L
  )
  expect_equal(res$bound, 6.4082048901, tolerance = 1e-9)
  expect_equal(res$incumbent, 6.4082048895, tolerance = 1e-9)
  expect_identical(res$iterations, 83L)
  expect_lt(res$bound - res$incumbent, 1e-9)
})
