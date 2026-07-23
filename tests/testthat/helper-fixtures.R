# Tiny shared fixtures and finite-difference helpers. Everything is kept small
# (K = 3, single-digit n) so the whole suite runs in a couple of seconds.
#
# The plurality null is built here, in the tests, purely from the general
# polytope_face API -- it is an *exemplar* of a union-structured null, not part
# of the package. Face j (j in 2:K) is the closed half-space {theta_j >=
# theta_1}, whose vertices are the basis points e_k (k != 1) plus the tie point
# (e_1 + e_j) / 2.
plurality_faces <- function(K) {
  lapply(2:K, function(j) {
    basis <- lapply(setdiff(seq_len(K), 1L), function(k) {
      v <- numeric(K)
      v[k] <- 1
      v
    })
    tie <- numeric(K)
    tie[c(1L, j)] <- 0.5
    polytope_face(
      vertices = do.call(cbind, c(basis, list(tie))),
      face_index = j
    )
  })
}

# A K = 3 finite-mixture alternative with candidate 1 leading (strictly inside
# the plurality alternative), expressed with the general mixture_dist.
small_alt <- function() {
  mixture_dist(
    components = matrix(c(0.6, 0.25, 0.15, 0.5, 0.3, 0.2), nrow = 3),
    weights = c(0.5, 0.5)
  )
}

# Mean parameter of a mixture_dist, a convenient reference for init points.
alt_mean <- function(alt) as.vector(alt@components %*% alt@weights)

# Central finite-difference gradient of a scalar function f: R^k -> R.
fd_grad <- function(f, x, eps = 1e-6) {
  vapply(seq_along(x), function(i) {
    xp <- x
    xm <- x
    xp[i] <- xp[i] + eps
    xm[i] <- xm[i] - eps
    (f(xp) - f(xm)) / (2 * eps)
  }, numeric(1L))
}

# Brute-force E_Q[p_theta / P] over an enumerated multinomial support: a plain
# weighted sum, independent of the engine's log-space arithmetic.
brute_expect_ratio <- function(alt, fam, theta, log_P) {
  X <- support(fam)
  log_q <- dist_log_density(alt, fam, X)
  log_p <- log_density(fam, theta)
  sum(exp(log_q) * exp(log_p - log_P))
}

# Init atoms: project a reference point onto each face.
init_on_faces <- function(faces, ref) {
  do.call(cbind, lapply(faces, function(f) init_point(f, ref)))
}
