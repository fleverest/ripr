#' @include numerics.R
NULL

# Bernstein enclosure over a simplex, and branch and bound on top of it.
#
# For a multinomial random variable `X`, `E_theta[X]` is a polynomial in
# `theta` of degree `n_trials`, and over a simplex the multinomial basis
# *is* the degree-`n` Bernstein basis. So the coefficients are the realised
# values of `X` over the lattice support, the convex hull property bounds the
# polynomial by their range for free, and de Casteljau subdivision tightens it
# quadratically (Leroy, 2012) in the sub-simplex diameter. That is  what makes
# a *proven* upper bound available at all, as against the oracle's multi-start
# gradient ascent search, which only ever gives a lower bound.
#
# `dc_step()` runs one step of the de Casteljau algorithm
# (Prautzsch-Boehm-Paluszny 10.4, from the Bernstein recursion in 10.1);
# `dc_pyramid()` iterates it; `dc_child()` reads a subsimplex expansion off the
# pyramid (Leroy Algorithm 2.13, step 4). Subdivision and reparametrisation are
# both callers -- see PBP 11.3, which treats them together.

# --- Size guard for Bernstein -------------------------------------------------

#' For a multinomial random variable `X`, `E_theta[X]` takes the form of a
#' polynomial in `theta` of degree `n_trials`. Its Bernstein form has one
#' coefficient per point of the multinomial sample space.
#' @keywords internal
#' @noRd
bernstein_size <- function(n_trials, k) choose(n_trials + k - 1, k - 1)


#' The largest batch (`n_trials`) that would fit within a budget
#'
#' Reported in the refusal, since the batch size is one thing a caller can
#' actually change.
#' @keywords internal
#' @noRd
largest_batch <- function(k, max_coefficients) {
  # Binary search over n_trials over the support size
  if (bernstein_size(1L, k) > max_coefficients) {
    return(0L)
  }
  hi <- 1L
  while (bernstein_size(2L * hi, k) <= max_coefficients) {
    hi <- 2L * hi
  }
  lo <- hi
  hi <- 2L * hi
  while (lo < hi) {
    mid <- (lo + hi + 1L) %/% 2L
    if (bernstein_size(mid, k) <= max_coefficients) {
      lo <- mid
    } else {
      hi <- mid - 1L
    }
  }
  lo
}


#' Refuse a certification too large to attempt
#'
#' Checked before anything is built, since the lattice is where the memory and
#' most of the setup time go. Exact arithmetic on `choose()`, so the guard itself
#' costs nothing.
#'
#' This is a resource limit rather than a correctness one: raising it costs time
#' and memory and nothing else. The other guards, e.g. checking that the family
#' permits a bound on the expectation at all, are a correctness concern and can
#' not be overridden.
#' @keywords internal
#' @noRd
check_bernstein_size <- function(n_trials, k, max_coefficients) {
  size <- bernstein_size(n_trials, k)
  if (size <= max_coefficients) {
    return(invisible(size))
  }
  stop(
    "certifying would need ",
    format(size, big.mark = ",", scientific = FALSE),
    " Bernstein coefficients for n_trials = ",
    n_trials,
    " with k = ",
    k,
    ", above `max_coefficients` (",
    format(max_coefficients, big.mark = ",", scientific = FALSE),
    ").\n",
    "The count is choose(n_trials + k - 1, k - 1), so it is the batch size ",
    "that drives it: ",
    largest_batch(k, max_coefficients),
    " would fit.\n",
    "Raise `max_coefficients` to attempt it anyway.",
    call. = FALSE
  )
}


# ---- lattice ---------------------------------------------------------------

#' Enumerate the degree-`n` tally lattice on `K` categories
#'
#' Precomputes every index the de Casteljau routines need: the multi-index
#' matrix, the `K` vertex positions (PBP 10.2: `b(a_j) = b_{n e_j}`), the edge
#' list, the degree ladder used by `dc_step()`, and the read-off map used by
#' `dc_child()`.
#'
#' The row order of `tally` is the coefficient order from `compositions()`, so
#' be careful of the order of the coefficients upstream when implementing with
#' a new family.
#' @param n Degree (multinomial batch size). At least 1.
#' @param K Number of categories. At least 2.
#' @return List with `n`, `K`, `tally`, `n_coef`, `pw` (base-`(n+1)` place
#'   values), `vertex` (positions of the `K` vertex coefficients), `edges` (the
#'   `2 x choose(K, 2)` edge list) `up` and `readoff`.
#' @keywords internal
#' @noRd
bernstein_lattice <- function(n, K) {
  n <- as.integer(n)
  K <- as.integer(K)
  stopifnot(length(n) == 1L, length(K) == 1L, n >= 1L, K >= 2L)

  tally <- compositions(n, K)
  pw <- (n + 1)^(seq_len(K) - 1L)
  # Tally n * e_j has base-(n+1) key n * (n+1)^(j-1).
  vertex <- match(n * pw, as.vector(tally %*% pw))

  edges <- utils::combn(K, 2L)
  # `up[[m + 1]][beta, i]` is the position of `beta + e_i` among the degree-`m`
  # multi-indices, so one de Casteljau step is a gather and a matrix product.
  key <- lapply(0:n, function(m) as.vector(compositions(m, K) %*% pw))
  up <- vector("list", n + 1L)
  for (m in seq_len(n)) {
    idx <- matrix(0L, nrow = length(key[[m]]), ncol = K)
    for (i in seq_len(K)) {
      idx[, i] <- match(key[[m]] + pw[i], key[[m + 1L]])
    }
    up[[m + 1L]] <- idx
  }

  # Leroy Algorithm 2.13 step 4: `b_alpha(V^[i]) = b^(alpha_i)_{alpha-hat-i}`.
  # Level `l` of the pyramid holds degree `n - l`, and zeroing slot `i` leaves
  # degree `n - alpha_i`, so the level to read is `alpha_i`. Grouped by level so
  # the read-off is a handful of vectorised gathers.
  readoff <- lapply(seq_len(K), function(i) {
    lev <- tally[, i]
    hat <- tally
    hat[, i] <- 0L
    hkey <- as.vector(hat %*% pw)
    lapply(sort(unique(lev)), function(l) {
      rows <- which(lev == l)
      list(
        level = l + 1L,
        rows = rows,
        pos = match(hkey[rows], key[[n - l + 1L]])
      )
    })
  })

  list(
    n = n,
    K = K,
    tally = tally,
    n_coef = nrow(tally),
    pw = pw,
    vertex = vertex,
    edges = edges,
    up = up,
    readoff = readoff
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


# ---- primitives ------------------------------------------------------------

#' One de Casteljau step at barycentric weights `lambda`: degree `m` -> `m - 1`
#'
#' PBP 10.4: `b_i <- [b_{i+e_0} + ... + b_{i+e_d}] u`, which follows from the
#' Bernstein recursion `B^n_i = u_0 B^{n-1}_{i-e_0} + ... + u_d B^{n-1}_{i-e_d}`
#' in 10.1. Every entry is an affine combination of its parents, and a convex
#' one when `lambda >= 0`.
#' @keywords internal
#' @noRd
dc_step <- function(cur, m, lambda, lat) {
  idx <- lat$up[[m + 1L]]
  as.vector(matrix(cur[idx], nrow = nrow(idx)) %*% lambda)
}


#' The full de Casteljau pyramid
#'
#' PBP 10.4: `n` steps reduce the degree-`n` array to the single value at
#' `lambda`. Level `l` holds the degree-`n - l` intermediates `b^(l)`; the book
#' calls the collection a tetrahedral array. All the levels are kept, because
#' the subsimplex expansions are read off them.
#' @keywords internal
#' @noRd
dc_pyramid <- function(coef, lat, lambda) {
  levels <- vector("list", lat$n + 1L)
  levels[[1L]] <- coef
  for (l in seq_len(lat$n)) {
    levels[[l + 1L]] <- dc_step(levels[[l]], lat$n - l + 1L, lambda, lat)
  }
  levels
}


#' Read the expansion over `V^[i]` off a pyramid
#'
#' Leroy Algorithm 2.13, step 4. `V^[i]` is the simplex with vertex `i` replaced
#' by the point `lambda` was taken at, and its coefficients are
#' `b_alpha(V^[i]) = b^(alpha_i)_{alpha-hat_i}` -- one entry of one pyramid
#' level per output coefficient, no arithmetic.
#' @keywords internal
#' @noRd
dc_child <- function(pyr, lat, i) {
  out <- numeric(lat$n_coef)
  for (g in lat$readoff[[i]]) {
    out[g$rows] <- pyr[[g$level]][g$pos]
  }
  out
}


#' Split a box at the point with barycentric weights `lambda`
#'
#' PBP 11.3. Returns one child per non-degenerate `V^[i]`: `lambda_i = 0` puts
#' the new point in the face opposite vertex `i`, so `V^[i]` would be flat.
#' *All* the children come from one pyramid -- the value `b^(n)_0` sits in every
#' one of them, and its dependency cone is the whole pyramid, so computing one
#' child costs exactly what computing them all costs.
#'
#' A box is `list(V, coef)`: `V` is `K x K` with columns giving the
#' sub-simplex's vertices in barycentric coordinates of the original simplex
#' (which, for the probability simplex, are the parameter vectors themselves),
#' and `coef` are its Bernstein coefficients in `lat`'s row order.
#' @keywords internal
#' @noRd
subdivide <- function(box, lat, lambda) {
  pyr <- dc_pyramid(box$coef, lat, lambda)
  point <- box$V %*% lambda
  lapply(which(lambda != 0), function(i) {
    V <- box$V
    V[, i] <- point
    list(V = V, coef = dc_child(pyr, lat, i))
  })
}


#' Bisect a box's edge `(p, q)`, returning both children exactly
#'
#' Leroy Example 2.15: binary splitting at the midpoint of an edge. Only
#' `lambda_p` and `lambda_q` are non-zero, so exactly two children are
#' non-degenerate. Midpoints rather than arbitrary edge points because that is
#' what bounds the shrinking factor, and hence the subdivision count (Leroy
#' Lemma 2.16, Theorem 3.6).
#'
#' @param box A box.
#' @param p,q Edge endpoints.
#' @param lat A `bernstein_lattice()`.
#' @return The two child boxes, `V^[p]` then `V^[q]` -- that is, the one with
#'   vertex `p` *replaced* first. Note this is the opposite labelling to
#'   "keeps vertex `p`".
#' @keywords internal
#' @noRd
bisect <- function(box, p, q, lat) {
  lambda <- numeric(lat$K)
  lambda[c(p, q)] <- 0.5
  subdivide(box, lat, lambda)
}

#' The edge of a box with the greatest Euclidean length
#' @param V `K x K` vertex matrix.
#' @param edges `2 x E` edge list.
#' @return Length-2 integer vector of endpoints.
#' @keywords internal
#' @noRd
longest_edge <- function(V, edges) {
  d2 <- colSums(
    (V[, edges[1L, ], drop = FALSE] -
      V[, edges[2L, ], drop = FALSE])^2
  )
  edges[, which.max(d2)]
}

# PBP 10.2 (convex hull property) with 10.3 Remark 2 (functional surface, so
# coefficients are Bezier *ordinates*): i.e. G <= max coefficient over the box.
box_bound <- function(box) max(box$coef)

# PBP 10.2: `b(a_0) = b_{n0...0}, ..., b(a_d) = b_{0...0n}`. The vertex
# coefficients are exact values of G, hence a valid *lower* bound.
vertex_values <- function(box, lat) box$coef[lat$vertex]

# The best vertex of a box: its value and the parameter vector attaining it.
box_best <- function(box, lat) {
  v <- vertex_values(box, lat)
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

#' Reparametrise a Bernstein form onto an arbitrary sub-simplex
#'
#' PBP 11.2: the polar form `b[x_1 ... x_n]` is the unique symmetric multiaffine
#' map with `b[x ... x] = b(x)`, and its values at the vertex arguments
#'
#'   `b_alpha = b[v_1 ... v_1 v_2 ... v_2 ... v_K ... v_K]`, `v_j` taken
#'   `alpha_j` times,
#'
#' *are* the Bezier coefficients over `conv(v_1, ..., v_K)`. The recursion
#' consuming one argument per step is PBP 11.2 (1), which is `dc_step()`; when
#' every argument is the same point it collapses to de Casteljau's algorithm,
#' so a single moved vertex agrees exactly with plain subdivision.
#'
#' Cost is `choose(n + K, K)` `dc_step()` calls. This is called once per sunull
#' at the start before `certify_sup()`, never inside the branch-and-bound loop,
#' which subdivides with `bisect()` instead of full reparametrisation.
#'
#' @param coef Length-`n_coef` coefficient vector in `lat$tally` row order.
#' @param lat A `bernstein_lattice()`.
#' @param vertices `(K, K)` matrix whose columns are the new vertices in
#'   barycentric coordinates of the original simplex.
#' @return Coefficient vector over the new simplex, in `lat$tally` row order.
#' @keywords internal
#' @noRd
reparametrise_to <- function(coef, lat, vertices) {
  V <- as.matrix(vertices)
  stopifnot(
    nrow(V) == lat$K,
    ncol(V) == lat$K,
    all(is.finite(V)),
    length(coef) == lat$n_coef,
    "vertices must lie in the standard simplex" = all(V >= -1e-12) &&
      max(abs(colSums(V) - 1)) < 1e-9,
    "vertices must span a non-degenerate simplex" = abs(det(V)) > 1e-12
  )

  # `compositions()` loops the first coordinate ascending and recurses on the
  # rest, which is this recursion, so the sub-results concatenate directly into
  # coefficient order.
  fill <- function(cur, j, remaining) {
    if (j == lat$K) {
      for (deg in rev(seq_len(remaining))) {
        cur <- dc_step(cur, deg, V[, j], lat)
      }
      return(cur)
    }
    parts <- vector("list", remaining + 1L)
    for (a in 0:remaining) {
      if (a > 0L) {
        cur <- dc_step(cur, remaining - a + 1L, V[, j], lat)
      }
      parts[[a + 1L]] <- fill(cur, j + 1L, remaining - a)
    }
    unlist(parts, use.names = FALSE)
  }

  fill(as.numeric(coef), 1L, lat$n)
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
#'   converged, budget_hit, trace)`. `converged` and `budget_hit` are mutually
#'   exclusive and exactly one is `TRUE`: the search either pruned or closed the
#'   gap, or it ran out of `max_iter`. Only the second qualifies the bound.
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
  converged <- FALSE
  budget_hit <- FALSE

  repeat {
    bounds <- vapply(active, box_bound, numeric(1L))
    u <- if (length(bounds)) max(bounds) else -Inf
    # Pruned boxes were all bounded by `incumbent + slack`, so the maximum over
    # the whole domain is bounded by the larger of that and the active maximum.
    bound <- max(incumbent + slack, u) + eta

    if (!length(active)) {
      converged <- TRUE
      break
    }
    if (bound - eta - incumbent <= tol) {
      converged <- TRUE
      break
    }
    if (it >= max_iter) {
      budget_hit <- TRUE
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
      bounds >=
        incumbent - eta - 8 * .Machine$double.eps * max(1, abs(incumbent))
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
    converged = converged,
    budget_hit = budget_hit,
    trace = trace[seq_len(it)]
  )
}
