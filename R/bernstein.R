#' @include numerics.R
NULL

# Bernstein enclosure over a simplex, and branch and bound on top of it.
#
# `E_theta[X]` is a degree-`n_trials` polynomial in `theta`, and over a simplex
# the multinomial basis *is* the Bernstein basis: the coefficients are the
# values of `X` on the lattice, their range bounds the polynomial over the
# simplex, and de Casteljau subdivision tightens that quadratically in terms of
# the diameter of the sub-simplices.
# Leroy (2012), Reliable Computing 17(1), 11-21.
# "PBP" is shorthand for Prautzsch, Boehm and Paluszny (2002).

# --- Size guard for Bernstein -------------------------------------------------

# Number of Bernstein coefficients: one per point of the multinomial sample
# space.
bernstein_size <- function(n_trials, k) choose(n_trials + k - 1, k - 1)


# Reject a certification too large to attempt, before the lattice is built. A
# resource limit only: raising `max_coefficients` costs a lot of time and
# memory.
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
    "that drives it: reduce `n_trials`, or raise `max_coefficients` to ",
    "attempt it anyway.",
    call. = FALSE
  )
}


# ---- lattice ---------------------------------------------------------------

# Enumerate the degree-`n` tally lattice on `K` categories, with every index
# the de Casteljau routines need: `vertex` (PBP 10.2), `edges` (2 x E), `up`
# for `dc_step()`, `readoff` for `dc_child()`, and `pw` (base-`(n+1)` place
# values keying multi-indices). `tally` rows follow `enumerate_counts()`, the
# same order `enumerate_space()` lists a `count_space` in, so `x` on the
# multinomial sample space is already a coefficient vector. A new family must
# supply coefficients in this order.
bernstein_lattice <- function(n, K) {
  n <- as.integer(n)
  K <- as.integer(K)
  stopifnot(length(n) == 1L, length(K) == 1L, n >= 1L, K >= 2L)

  tally <- enumerate_counts(n, K)
  pw <- (n + 1)^(seq_len(K) - 1L)
  # Tally n * e_j has base-(n+1) key n * (n+1)^(j-1).
  vertex <- match(n * pw, as.vector(tally %*% pw))

  edges <- utils::combn(K, 2L)
  # `up[[m + 1]][beta, i]` is the position of `beta + e_i` among the degree-`m`
  # multi-indices, so one de Casteljau step is a gather and a matrix product.
  key <- lapply(0:n, function(m) as.vector(enumerate_counts(m, K) %*% pw))
  up <- vector("list", n + 1L)
  for (m in seq_len(n)) {
    idx <- matrix(0L, nrow = length(key[[m]]), ncol = K)
    for (i in seq_len(K)) {
      idx[, i] <- match(key[[m]] + pw[i], key[[m + 1L]])
    }
    up[[m + 1L]] <- idx
  }

  # Leroy Algorithm 2.13 step 4: `b_alpha(V^[i]) = b^(alpha_i)_{alpha-hat-i}`,
  # read from pyramid level `alpha_i`. Grouped by level for vectorised gathers.
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

# ---- primitives ------------------------------------------------------------

# One de Casteljau step (PBP 10.4) at barycentric weights `lambda`, degree `m`
# to `m - 1`. Convex combinations when `lambda >= 0`.
dc_step <- function(cur, m, lambda, lat) {
  idx <- lat$up[[m + 1L]]
  as.vector(matrix(cur[idx], nrow = nrow(idx)) %*% lambda)
}


# All `n` de Casteljau levels; level `l` holds degree `n - l`. Every level is
# kept because the subsimplex expansions are read off them.
dc_pyramid <- function(coef, lat, lambda) {
  levels <- vector("list", lat$n + 1L)
  levels[[1L]] <- coef
  for (l in seq_len(lat$n)) {
    levels[[l + 1L]] <- dc_step(levels[[l]], lat$n - l + 1L, lambda, lat)
  }
  levels
}


# The expansion over `V^[i]` (vertex `i` replaced by the split point), read
# off the pyramid without arithmetic (Leroy Algorithm 2.13, step 4).
dc_child <- function(pyr, lat, i) {
  out <- numeric(lat$n_coef)
  for (g in lat$readoff[[i]]) {
    out[g$rows] <- pyr[[g$level]][g$pos]
  }
  out
}


# Split a box at barycentric weights `lambda` (PBP 11.3), one child per
# `lambda_i != 0` (the others would be flat).
#
# A box is `list(V, coef)`: `V` is `K x K`, rows the sub-simplex's vertices in
# barycentric coordinates of the original (for the probability simplex, the
# parameter vectors themselves); `coef` in `lat`'s row order.
subdivide <- function(box, lat, lambda) {
  pyr <- dc_pyramid(box$coef, lat, lambda)
  point <- drop(lambda %*% box$V)
  lapply(which(lambda != 0), function(i) {
    V <- box$V
    V[i, ] <- point
    list(V = V, coef = dc_child(pyr, lat, i))
  })
}


# Bisect edge `(p, q)` at its midpoint (Leroy Example 2.15); midpoints bound
# the shrinking factor and hence the subdivision count (Lemma 2.16, Thm 3.6).
# Returns `V^[p]` then `V^[q]`: the child with vertex `p` *replaced* first.
bisect <- function(box, p, q, lat) {
  lambda <- numeric(lat$K)
  lambda[c(p, q)] <- 0.5
  subdivide(box, lat, lambda)
}

# Endpoints of the longest edge of vertex matrix `V`.
longest_edge <- function(V, edges) {
  d2 <- rowSums(
    (V[edges[1L, ], , drop = FALSE] -
      V[edges[2L, ], , drop = FALSE])^2
  )
  edges[, which.max(d2)]
}

# Convex hull property (PBP 10.2, 10.3 Remark 2): G <= max coefficient.
box_bound <- function(box) max(box$coef)

# Vertex coefficients are exact values of G (PBP 10.2), hence a lower bound.
vertex_values <- function(box, lat) box$coef[lat$vertex]

box_best <- function(box, lat) {
  v <- vertex_values(box, lat)
  j <- which.max(v)
  list(value = v[[j]], theta = box$V[j, ])
}

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

# How the rows of `v` leave the standard simplex, as a clause for an error
# message, or `NULL` if they don't. The one membership test for vertex
# matrices; negativity is reported before the sum.
simplex_departure <- function(v, neg_tol = 1e-12, sum_tol = 1e-9) {
  if (any(v < -neg_tol)) {
    return(paste0("the smallest coordinate is ", format(min(v))))
  }
  sums <- rowSums(v)
  off <- abs(sums - 1)
  if (max(off) >= sum_tol) {
    return(paste0(
      "the coordinates of one sum to ",
      format(sums[which.max(off)]),
      " rather than 1"
    ))
  }
  NULL
}

# A lower-dimensional simplex as `K` vertex rows, by repeating its last vertex.
pad_vertices <- function(V, K) {
  V[c(seq_len(nrow(V)), rep(nrow(V), K - nrow(V))), , drop = FALSE]
}

# Bernstein coefficients over the sub-simplex with (barycentric) vertex rows
# `vertices`, via the polar form (PBP 11.2): `b_alpha` is the blossom with
# `v_j` taken `alpha_j` times, each argument consumed by one `dc_step()`.
reparametrise_to <- function(coef, lat, vertices) {
  V <- as.matrix(vertices)
  stopifnot(
    nrow(V) == lat$K,
    ncol(V) == lat$K,
    all(is.finite(V)),
    length(coef) == lat$n_coef,
    "vertices must lie in the standard simplex" = is.null(
      simplex_departure(V)
    )
  )

  # This recursion visits multi-indices in `enumerate_counts()` order, so the
  # parts concatenate straight into coefficient order.
  fill <- function(cur, j, remaining) {
    if (j == lat$K) {
      for (deg in rev(seq_len(remaining))) {
        cur <- dc_step(cur, deg, V[j, ], lat)
      }
      return(cur)
    }
    parts <- vector("list", remaining + 1L)
    for (a in 0:remaining) {
      if (a > 0L) {
        cur <- dc_step(cur, remaining - a + 1L, V[j, ], lat)
      }
      parts[[a + 1L]] <- fill(cur, j + 1L, remaining - a)
    }
    unlist(parts, use.names = FALSE)
  }

  fill(as.numeric(coef), 1L, lat$n)
}

# ---- branch and bound ------------------------------------------------------

# A branch-and-bound node: a box plus its cached upper bound and lineage.
node <- function(box, id, parent = NA_integer_, depth = 0L, born = 0L) {
  box$ub <- box_bound(box)
  box$id <- id
  box$parent <- parent
  box$depth <- depth
  # With `parent` and the record's `retired`, makes the run replayable.
  box$born <- born
  box
}

# A node without its coefficients, for the history.
node_stub <- function(box, retired, fate) {
  list(
    id = box$id,
    parent = box$parent,
    depth = box$depth,
    born = box$born,
    retired = retired,
    fate = fate,
    ub = box$ub,
    V = box$V
  )
}

node_ubs <- function(nodes) vapply(nodes, function(b) b$ub, numeric(1L))


# Why the search should stop, or `NULL`. Only "budget_hit" qualifies the
# result.
stop_reason <- function(n_active, gap, tol, it, max_iter) {
  if (n_active == 0L) {
    return("converged")
  }
  if (gap <= tol) {
    return("converged")
  }
  if (it >= max_iter) {
    return("budget_hit")
  }
  NULL
}


# Global upper bound on `sup G` over one sub-simplex.
#
# `bound` is valid (up to rounding) at every iteration. It includes the largest
# bound ever pruned, since a pruned node may still hold the supremum.
# `shared_incumbent` is a value attained elsewhere (e.g. another cell of the
# same null), used only for pruning and stopping. Exactly one of `converged` and
# `budget_hit` is `TRUE`.
certify_sup <- function(
  box,
  lat,
  tol = 1e-3,
  max_iter = 500L,
  shared_incumbent = -Inf
) {
  max_iter <- as.integer(max_iter)

  active <- list(node(box, id = 1L))
  ubs <- node_ubs(active)
  next_id <- 2L
  best <- box_best(box, lat)
  incumbent_trace <- numeric(max_iter)
  retired_nodes <- list()
  trace <- numeric(max_iter)
  it <- 0L
  reason <- NULL

  pruned_ub <- -Inf

  repeat {
    incumbent <- max(best$value, shared_incumbent)
    bound <- max(best$value, ubs, pruned_ub)

    reason <- stop_reason(
      length(active),
      bound - incumbent,
      tol,
      it,
      max_iter
    )
    if (!is.null(reason)) {
      break
    }

    it <- it + 1L

    j <- which.max(ubs)
    parent <- active[[j]]
    e <- longest_edge(parent$V, lat$edges)
    kids <- bisect(parent, e[1L], e[2L], lat)
    kids <- lapply(seq_along(kids), function(i) {
      node(
        kids[[i]],
        id = next_id + i - 1L,
        parent = parent$id,
        depth = parent$depth + 1L,
        born = it
      )
    })
    next_id <- next_id + length(kids)
    retired_nodes[[length(retired_nodes) + 1L]] <-
      node_stub(parent, it, "split")
    active <- c(active[-j], kids)
    ubs <- c(ubs[-j], node_ubs(kids))
    kid_best <- boxes_best(kids, lat)
    if (kid_best$value > best$value) {
      best <- kid_best
    }

    # Leroy Lemma 3.2 cut-off: drop nodes whose bound can't beat the incumbent
    # (ties too).
    incumbent <- max(best$value, shared_incumbent)
    keep <- ubs > incumbent
    if (!all(keep)) {
      pruned_ub <- max(pruned_ub, ubs[!keep])
      for (b in active[!keep]) {
        retired_nodes[[length(retired_nodes) + 1L]] <-
          node_stub(b, it, "pruned")
      }
      active <- active[keep]
      ubs <- ubs[keep]
    }
    trace[it] <- max(best$value, ubs, pruned_ub)
    incumbent_trace[it] <- best$value
  }

  list(
    bound = bound,
    incumbent = best$value,
    theta = best$theta,
    iterations = it,
    history = c(
      retired_nodes,
      lapply(active, node_stub, retired = NA_integer_, fate = "active")
    ),
    converged = identical(reason, "converged"),
    budget_hit = identical(reason, "budget_hit"),
    trace = trace[seq_len(it)],
    incumbent_trace = incumbent_trace[seq_len(it)]
  )
}
