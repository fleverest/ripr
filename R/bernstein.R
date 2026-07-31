# Certified upper bounds on `sup_theta G` over a sub-simplex of the probability
# simplex, by branch and bound in the simplicial Bernstein basis.
#
# For a batch of `n` draws over `K` categories, `G(theta) = E_Q[p_theta / P_W]`
# expands in the multinomial basis, which *is* the degree-`n` Bernstein basis on
# the simplex. The Bernstein coefficients are then the realised likelihood
# ratios `q(x) / p_W(x)`, so the convex hull property gives a valid upper bound
# for free and de Casteljau subdivision refines it, quadratically in the
# sub-simplex diameter. See Leroy (2012), Reliable Computing 17(1), 11-21.
#
# Deliberately base R: this is a numerical inner loop over preallocated numeric
# vectors. Nothing here knows about faces, engines or families -- the seam to
# the rest of the package is oracle_bound.R.

# ---- lattice ---------------------------------------------------------------

#' Enumerate the degree-`n` tally lattice on `K` categories
#'
#' Precomputes the index structure the subdivision needs: the multi-index
#' matrix, the positions of the `K` vertex coefficients, the edge list, and --
#' for each edge -- the lattice lines parallel to it.
#'
#' The row order of `tally` is the coefficient order for every box built against
#' this lattice, and it is *not* assumed: `vertex` is found by key `match()` and
#' `rows` by `split()`, so any valid enumeration of the degree-`n` lattice
#' works. Callers integrating against an engine must pass that engine's own
#' outcome matrix, so the coefficients they build from `engine@log_q_mass` are
#' indexed consistently.
#'
#' @param n Degree (multinomial batch size). At least 1.
#' @param K Number of categories. At least 2.
#' @param tally Optional `(M, K)` count matrix giving the enumeration to use.
#'   `NULL` (default) uses the internal `compositions()` enumeration.
#' @return List with `n`, `K`, `tally`, `n_coef`, `pw` (base-`(n+1)` place
#'   values), `vertex` (positions of the `K` vertex coefficients), `edges` (the
#'   `2 x choose(K, 2)` edge list) and `rows`.
#' @keywords internal
#' @noRd
bernstein_lattice <- function(n, K, tally = NULL) {
  n <- as.integer(n)
  K <- as.integer(K)
  stopifnot(length(n) == 1L, length(K) == 1L, n >= 1L, K >= 2L)

  if (is.null(tally)) {
    tally <- compositions(n, K)
  } else {
    tally <- as.matrix(tally)
    storage.mode(tally) <- "integer"
    dimnames(tally) <- NULL
    # A wrong enumeration would mis-index every coefficient and yield a
    # plausible but invalid bound, so this is checked rather than assumed.
    # Non-negativity is part of the check because it is what makes the
    # base-(n+1) key injective.
    stopifnot(
      ncol(tally) == K,
      nrow(tally) == choose(n + K - 1L, K - 1L),
      all(is.finite(tally)),
      all(tally >= 0L),
      all(rowSums(tally) == n),
      !anyDuplicated(as.vector(tally %*% (n + 1)^(seq_len(K) - 1L)))
    )
  }

  # Tally n * e_j has base-(n+1) key n * (n+1)^(j-1).
  pw <- (n + 1)^(seq_len(K) - 1L)
  vertex <- match(n * pw, as.vector(tally %*% pw))

  edges <- utils::combn(K, 2L)
  rows <- apply(
    edges,
    2L,
    function(e) build_rows(tally, e[1L], e[2L], n),
    simplify = FALSE
  )
  names(rows) <- paste(edges[1L, ], edges[2L, ], sep = "-")

  list(
    n = n,
    K = K,
    tally = tally,
    n_coef = nrow(tally),
    pw = pw,
    vertex = vertex,
    edges = edges,
    rows = rows
  )
}

#' Every composition of `n` into `k` non-negative parts, lexicographically
#'
#' Prefixes `i = 0:n` and recurses, which emits the multi-indices in ascending
#' lexicographic order. Agrees row-for-row with the stars-and-bars enumeration
#' in `build_counts_matrix()`; `test-bernstein.R` pins that so a change to
#' either fails loudly rather than silently mis-indexing coefficients.
#'
#' @param n Total to be split.
#' @param k Number of parts.
#' @return `(choose(n + k - 1, k - 1), k)` integer matrix, no dimnames.
#' @keywords internal
#' @noRd
compositions <- function(n, k) {
  if (k == 1L) {
    return(matrix(as.integer(n), nrow = 1L, ncol = 1L))
  }
  out <- do.call(
    rbind,
    lapply(0:n, function(i) cbind(i, compositions(n - i, k - 1L)))
  )
  dimnames(out) <- NULL
  storage.mode(out) <- "integer"
  out
}

#' Positions of the lattice lines parallel to edge `(p, q)`
#'
#' Each line is returned as a position vector ordered by increasing power of
#' vertex `q`. Bisecting `(p, q)` uses barycentric weights supported on
#' `{p, q}`, so the other slots are spectators and the subdivision decouples
#' into one 1-D de Casteljau pass per line.
#'
#' @param tally `(M, K)` multi-index matrix.
#' @param p,q Edge endpoints (category indices).
#' @param n Degree.
#' @return List of integer position vectors.
#' @keywords internal
#' @noRd
build_rows <- function(tally, p, q, n) {
  spec <- tally[, -c(p, q), drop = FALSE]
  key <- if (ncol(spec) == 0L) {
    rep.int(0, nrow(tally))
  } else {
    as.vector(spec %*% (n + 1)^(seq_len(ncol(spec)) - 1L))
  }
  lapply(
    unname(split(seq_len(nrow(tally)), key)),
    function(i) i[order(tally[i, q])]
  )
}

# ---- primitives ------------------------------------------------------------

#' de Casteljau midpoint split of a 1-D Bernstein array
#'
#' Exact: every entry of the output is a convex combination of the input, so
#' this is unconditionally stable.
#'
#' @param b Numeric vector of degree-`length(b) - 1` Bernstein coefficients.
#' @return `list(left = , right = )`, each of the same length as `b`.
#' @keywords internal
#' @noRd
split_1d <- function(b) {
  m <- length(b) - 1L
  if (m <= 0L) {
    return(list(left = b, right = b))
  }
  left <- right <- numeric(m + 1L)
  left[1L] <- b[1L]
  right[m + 1L] <- b[m + 1L]
  cur <- b
  for (r in seq_len(m)) {
    cur <- 0.5 * (cur[-length(cur)] + cur[-1L])
    left[r + 1L] <- cur[1L]
    right[m + 1L - r] <- cur[length(cur)]
  }
  list(left = left, right = right)
}

#' Bisect a box's edge `(p, q)`, returning both children exactly
#'
#' A box is `list(V, coef)`: `V` is `K x K` with columns giving the
#' sub-simplex's vertices in barycentric coordinates of the original simplex
#' (which, for the probability simplex, are the parameter vectors themselves),
#' and `coef` are its Bernstein coefficients in `lat`'s row order.
#'
#' @param box A box.
#' @param p,q Edge endpoints.
#' @param lat A `bernstein_lattice()`.
#' @return List of the two child boxes: the one keeping vertex `p`, then the one
#'   keeping vertex `q`.
#' @keywords internal
#' @noRd
bisect <- function(box, p, q, lat) {
  cl <- cr <- box$coef
  for (line in lat$rows[[paste(p, q, sep = "-")]]) {
    lr <- split_1d(box$coef[line])
    cl[line] <- lr$left
    cr[line] <- lr$right
  }
  w <- 0.5 * (box$V[, p] + box$V[, q])
  v_left <- v_right <- box$V
  v_left[, q] <- w # child keeping vertex p
  v_right[, p] <- w # child keeping vertex q
  list(list(V = v_left, coef = cl), list(V = v_right, coef = cr))
}

#' The edge of a box with the greatest Euclidean length
#' @param V `K x K` vertex matrix.
#' @param edges `2 x E` edge list.
#' @return Length-2 integer vector of endpoints.
#' @keywords internal
#' @noRd
longest_edge <- function(V, edges) {
  d2 <- colSums((V[, edges[1L, ], drop = FALSE] -
    V[, edges[2L, ], drop = FALSE])^2)
  edges[, which.max(d2)]
}

# The convex hull property: G <= max coefficient over the box.
box_bound <- function(box) max(box$coef)

# Bernstein forms interpolate at the vertices, so the vertex coefficients are
# exact values of G -- hence a valid *lower* bound on the supremum.
box_values <- function(box, lat) box$coef[lat$vertex]

# The best vertex of a box: its value and the parameter vector attaining it.
box_best <- function(box, lat) {
  v <- box_values(box, lat)
  j <- which.max(v)
  list(value = v[[j]], theta = box$V[, j])
}

# The best vertex over a list of boxes.
boxes_best <- function(boxes, lat) {
  best <- list(value = -Inf, theta = NULL)
  for (b in boxes) {
    cand <- box_best(b, lat)
    if (cand$value > best$value) {
      best <- cand
    }
  }
  best
}

# ---- general reparametrisation ---------------------------------------------

#' Index structure for the intermediate degrees `0, ..., n`
#'
#' `up[[m + 1]]` is an `(n_coef(m - 1), K)` integer matrix whose `(beta, i)`
#' entry is the position of `beta + e_i` among the degree-`m` multi-indices, so
#' one de Casteljau step is a single matrix-vector product. Built once per
#' lattice and reused across every face.
#'
#' @param n Degree.
#' @param K Number of categories.
#' @return `list(key = , up = )`, both indexed by degree + 1.
#' @keywords internal
#' @noRd
degree_ladder <- function(n, K) {
  pw <- (n + 1)^(seq_len(K) - 1L)
  key <- vector("list", n + 1L)
  for (m in 0:n) {
    key[[m + 1L]] <- as.vector(compositions(m, K) %*% pw)
  }
  up <- vector("list", n + 1L)
  for (m in seq_len(n)) {
    idx <- matrix(0L, nrow = length(key[[m]]), ncol = K)
    for (i in seq_len(K)) {
      idx[, i] <- match(key[[m]] + pw[i], key[[m + 1L]])
    }
    up[[m + 1L]] <- idx
  }
  list(key = key, up = up)
}

# One de Casteljau step at barycentric weights `lambda`: degree m -> degree m-1.
dc_step <- function(cur, m, lambda, ladder) {
  idx <- ladder$up[[m + 1L]]
  as.vector(matrix(cur[idx], nrow = nrow(idx)) %*% lambda)
}

#' Reparametrise a Bernstein form onto an arbitrary sub-simplex
#'
#' Leroy (2012) section 2.3: `K` successive runs of de Casteljau (his Algorithm
#' 2.13) at the target vertices, with general barycentric weights. The
#' coefficient at multi-index `a` is the blossom of `G` evaluated at the new
#' vertices with multiplicities `a`, so the result is the *exact* Bernstein form
#' of the same polynomial over the new simplex.
#'
#' The row-decomposition shortcut `bisect()` uses does not apply here, because a
#' general vertex has no zero barycentric component: this is a full pyramid at
#' `O(K * n * choose(n + K - 1, K - 1))`. It is paid once per face at setup, not
#' once per bisection.
#'
#' @param coef Length-`n_coef` coefficient vector in `lat$tally` row order.
#' @param lat A `bernstein_lattice()`.
#' @param vertices `(K, K)` matrix whose columns are the new vertices in
#'   barycentric coordinates of the original simplex.
#' @return Coefficient vector over the new simplex, in `lat$tally` row order.
#' @keywords internal
#' @noRd
reparametrise_to <- function(coef, lat, vertices) {
  K <- lat$K
  n <- lat$n
  L <- as.matrix(vertices)
  coef <- as.numeric(coef)
  stopifnot(
    nrow(L) == K,
    ncol(L) == K,
    all(is.finite(L)),
    length(coef) == lat$n_coef
  )

  ladder <- degree_ladder(n, K)
  # `walk()` emits leaves in ascending lexicographic order of the multi-index,
  # which is exactly `compositions()`' order; `perm` carries that canonical
  # order to and from the (arbitrary) row order of `lat$tally`.
  perm <- match(ladder$key[[n + 1L]], as.vector(lat$tally %*% lat$pw))
  stopifnot(!anyNA(perm))

  res <- numeric(lat$n_coef)
  pos <- 0L

  # Replace the argument slots of the blossom one direction at a time: after
  # `a_j` steps in direction `j` the array holds b[v^alpha, w_1^a_1, ..., w_j^a_j].
  walk <- function(arr, m, j) {
    if (j == K) {
      cur <- arr
      mm <- m
      while (mm > 0L) {
        cur <- dc_step(cur, mm, L[, K], ladder)
        mm <- mm - 1L
      }
      pos <<- pos + 1L
      res[pos] <<- cur
      return(invisible(NULL))
    }
    cur <- arr
    mm <- m
    walk(cur, mm, j + 1L)
    while (mm > 0L) {
      cur <- dc_step(cur, mm, L[, j], ladder)
      mm <- mm - 1L
      walk(cur, mm, j + 1L)
    }
    invisible(NULL)
  }
  walk(coef[perm], n, 1L)

  out <- numeric(lat$n_coef)
  out[perm] <- res
  out
}

# ---- branch and bound ------------------------------------------------------

#' Certified upper bound on `sup G` over the union of the seed sub-simplices
#'
#' Validity does not depend on convergence: `bound` is a valid upper bound at
#' every iteration, so `G / bound <= 1` whenever you stop. Refinement buys
#' tightness, not validity.
#'
#' @param seeds List of boxes (see `bisect()`).
#' @param lat A `bernstein_lattice()`.
#' @param tol Stop once `bound - incumbent <= tol`.
#' @param max_iter Cap on bisections.
#' @param slack Prune with `U(S) <= incumbent + slack`. Cheaper, but the active
#'   set then encloses the slack-superlevel set rather than the argmax.
#' @param keep_argmax Prune with `U(S) < incumbent` instead, retaining ties, so
#'   that `active` is a certified enclosure of every maximiser. Forces
#'   `slack = 0`; the two modes are mutually exclusive.
#' @param round_slack Add `n * eps * max|coef|` to the reported bound. de
#'   Casteljau is all convex combinations and so is stable, but round-to-nearest
#'   can put the computed maximum coefficient marginally *below* the true one,
#'   which is the unsafe direction for a validity claim.
#' @return `list(bound, incumbent, theta, active, rejected, iterations,
#'   exhausted, trace)`.
#' @keywords internal
#' @noRd
certify_sup <- function(
  seeds,
  lat,
  tol = 1e-3,
  max_iter = 500L,
  slack = 0,
  keep_argmax = FALSE,
  round_slack = TRUE
) {
  if (keep_argmax && slack > 0) {
    stop("`slack` must be 0 when keep_argmax = TRUE")
  }
  max_iter <- as.integer(max_iter)
  active <- seeds
  best <- boxes_best(active, lat)
  incumbent <- best$value
  # Children's coefficients are convex combinations of their parent's, so the
  # seed maximum bounds |coef| over every box the run will ever create.
  eta <- if (round_slack) {
    lat$n *
      .Machine$double.eps *
      max(vapply(seeds, function(b) max(abs(b$coef)), numeric(1L)))
  } else {
    0
  }
  rejected <- list()
  trace <- numeric(max_iter)
  it <- 0L
  exhausted <- FALSE

  repeat {
    bounds <- vapply(active, box_bound, numeric(1L))
    u <- if (length(bounds)) max(bounds) else -Inf
    # Pruned boxes were all bounded by `incumbent + slack`, so the maximum over
    # the whole domain is bounded by the larger of that and the active maximum.
    bound <- max(incumbent + slack, u) + eta

    if (!length(active)) {
      exhausted <- TRUE
      break
    }
    if (bound - eta - incumbent <= tol || it >= max_iter) {
      break
    }

    it <- it + 1L
    j <- which.max(bounds)
    e <- longest_edge(active[[j]]$V, lat$edges)
    kids <- bisect(active[[j]], e[1L], e[2L], lat)
    active <- c(active[-j], kids)
    kid_best <- boxes_best(kids, lat)
    if (kid_best$value > best$value) {
      best <- kid_best
    }
    incumbent <- best$value

    bounds <- vapply(active, box_bound, numeric(1L))
    keep <- if (keep_argmax) {
      # Retain ties, generously: dropping a box that attains the maximum would
      # break the enclosure claim, whereas keeping a spare one only costs work.
      bounds >= incumbent - eta - 8 * .Machine$double.eps * max(1, abs(incumbent))
    } else {
      bounds > incumbent + slack
    }
    rejected <- c(rejected, active[!keep])
    active <- active[keep]
    trace[it] <- max(
      incumbent + slack,
      if (any(keep)) max(bounds[keep]) else -Inf
    ) +
      eta
  }

  list(
    bound = bound,
    incumbent = incumbent,
    theta = best$theta,
    active = active,
    rejected = rejected,
    iterations = it,
    exhausted = exhausted,
    trace = trace[seq_len(it)]
  )
}

# ---- a priori cost model ---------------------------------------------------

#' Largest mixed second difference of a box's Bernstein coefficients
#'
#' Definition 2.8 of Leroy (2012); the index `e_0` wraps to `e_K`. Adding
#' `e_a + e_b` to a tally shifts its base-`(n+1)` key by `pw[a] + pw[b]` with no
#' carrying, since `|gamma| = n - 2` keeps every slot at most `n`.
#'
#' @param box A box.
#' @param lat A `bernstein_lattice()`.
#' @return Numeric scalar.
#' @keywords internal
#' @noRd
second_difference_norm <- function(box, lat) {
  K <- lat$K
  n <- lat$n
  pw <- lat$pw
  key <- as.vector(lat$tally %*% pw)
  gkey <- as.vector(compositions(n - 2L, K) %*% pw)
  wrap <- function(i) if (i == 0L) K else i
  coef_at <- function(kk) box$coef[match(kk, key)]

  pairs <- utils::combn(K, 2L)
  diffs <- apply(pairs, 2L, function(ij) {
    i <- ij[1L]
    j <- ij[2L]
    max(abs(
      coef_at(gkey + pw[i] + pw[wrap(j - 1L)]) +
        coef_at(gkey + pw[wrap(i - 1L)] + pw[j]) -
        coef_at(gkey + pw[wrap(i - 1L)] + pw[wrap(j - 1L)]) -
        coef_at(gkey + pw[i] + pw[j])
    ))
  })
  max(diffs)
}

#' Leroy (2012) Theorem 3.6: worst-case subdivision effort for accuracy `tol`
#'
#' For a priori capacity planning only. It ignores the cut-off test entirely
#' (his Remark 3.7) and ran roughly 70x pessimistic against measured expansions
#' on the reference examples, so user-facing cost figures should come from the
#' measured `iterations` instead.
#'
#' `ratio` is the required `shrink^(-2N)`; `mesh` follows from Theorem 2.17 with
#' `m(Delta) = sqrt(2)`, and the sub-simplex count from covering at that mesh.
#'
#' @param lat A `bernstein_lattice()`.
#' @param sd_norm A `second_difference_norm()`.
#' @param tol Target accuracy.
#' @param shrink Per-step diameter shrink factor.
#' @return `list(dimension_constant, steps, mesh, subsimplices)`.
#' @keywords internal
#' @noRd
leroy_cost <- function(lat, sd_norm, tol, shrink = 0.5) {
  k <- lat$K - 1L
  const <- k^2 * (k + 1) * (k + 2)^2 * (k + 3) / 288
  ratio <- lat$n * const * sd_norm / tol
  list(
    dimension_constant = const,
    steps = log(ratio) / (2 * log(1 / shrink)),
    mesh = sqrt(2 / ratio),
    subsimplices = ratio^(k / 2)
  )
}
