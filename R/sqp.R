# Sequential quadratic programming solver
#
# A general exact-Hessian SQP (SLSQP from nloptr uses a quasi-Newton
# approximation) for smooth convex objectives subject to x >= 0.

#' Minimise a smooth convex function over `x >= 0` by exact-Hessian SQP
#'
#' Each iteration minimises the quadratic model
#' `f(x) + g'(z - x) + (z - x)' H (z - x) / 2` over `z >= 0` with
#' `qp_nonneg()`, warm-started from `x`, then backtracks (Armijo) along the
#' segment from `x` to that minimiser, which stays feasible since the orthant
#' is convex. Once the active set settles the steps are Newton steps on the
#' free coordinates, so convergence is locally quadratic.
#'
#' It stops when the KKT residual `max(|g_i|` for `x_i > 0`, `-g_i` for
#' `x_i = 0)` is at most `tol`, after `max_iter` iterations, or when no step
#' decreases `fn`. Near the optimum `fn` changes below rounding, so a full
#' step is also accepted when its predicted decrease is itself at rounding
#' level and `fn` does not rise beyond rounding.
#'
#' `x` must be feasible with `fn(x)` finite; `fn` may return `Inf` outside its
#' domain, which the line search then avoids. `hess(x)` should be positive
#' definite (add a small ridge if it may be singular), though `qp_nonneg()`
#' falls back to a least-squares solve on a singular free block. Returns
#' `list(par, value, kkt, iterations)`; `kkt` above `tol` means it stopped
#' short.
#'
#' References: Nocedal & Wright (2006). Specifically:
#'  - §16.5 Active-set methods for convex QPs, and
#'  - §18 Sequential Quadratic Programming
#'
#' Written by Claude Opus 5.5 (Anthropic), checked by Floyd against the above
#' @keywords internal
#' @noRd
sqp_nonneg <- function(x, fn, gr, hess, tol = 1e-12, max_iter = 500L) {
  # KKT residuals for x >= 0 (see N&W §12.3, Thm 12.1 for definition)
  # |g| where free, the positive part of -g at the bound.
  residual <- function(x, g) {
    on <- x > 0
    max(c(abs(g[on]), -g[!on], 0))
  }

  f <- fn(x)
  iterations <- 0L
  repeat {
    # Stop when converged, or out of iterations
    g <- gr(x)
    kkt <- residual(x, g)
    if (kkt <= tol || iterations >= max_iter) {
      break
    }

    # SQP subproblem (see N&W §18): minimise the quadratic model over
    # z >= 0 warm-started at x, by the primal active-set method (N&W §16.5).
    h <- hess(x)
    step <- qp_nonneg(h, g - as.vector(h %*% x), x) - x

    # Check that the directional derivative is negative (and finite)
    descent <- sum(g * step)
    if (!is.finite(descent) || descent >= 0) {
      break
    }

    # Backtracking line search (N&W §3.1, Alg. 3.1) along the feasible segment
    # from x to the QP minimiser.
    rounding <- 1e-14 * max(1, abs(f))
    t <- 1
    accepted <- FALSE
    step_size <- max(abs(step))
    scale <- max(1, abs(x))
    while (t * step_size > .Machine$double.eps * scale) {
      # Loop stops when the step has become negligible relative to x
      x_new <- pmax(x + t * step, 0) # clip rounding negatives to exact zeros
      f_new <- fn(x_new)
      # Armijo sufficient decrease (N&W eq. 3.4), or a full step whose
      # decrease is below rounding (the end game).
      sufficient <- f_new <= f + 1e-4 * t * descent
      negligible <- t == 1 && -descent < rounding && f_new <= f + rounding
      if (is.finite(f_new) && (sufficient || negligible)) {
        accepted <- TRUE
        break
      }
      t <- t / 2
    }
    if (!accepted) {
      break
    }

    # Take the step.
    x <- x_new
    f <- f_new
    iterations <- iterations + 1L
  }
  list(par = x, value = f, kkt = kkt, iterations = iterations)
}


#' Minimise `z' H z / 2 + c' z` over `z >= 0`, from a feasible `z`
#'
#' A primal active-set method: solve on the free set, step back to the
#' boundary if that leaves the orthant, otherwise free the bound coordinate
#' with the most negative gradient, until none is negative. Finite.
#'
#' It is very similar to using `quadprog::solve.QP`, but that requires a
#' p.d. Hessian. Our hessians will often be nearly singular (when atoms are
#' close to one another).
#' @keywords internal
#' @noRd
qp_nonneg <- function(h, c, z, max_iter = 10L * length(c) + 50L) {
  n <- length(c)
  free <- which(z > 0)
  for (it in seq_len(max_iter)) {
    if (length(free)) {
      h_free <- h[free, free, drop = FALSE]
      z_free <- tryCatch(
        solve(h_free, -c[free]),
        error = function(e) {
          sol <- qr.coef(qr(h_free), -c[free])
          ifelse(is.na(sol), 0, sol)
        }
      )
      if (all(z_free > 0)) {
        z[] <- 0
        z[free] <- z_free
      } else {
        p <- z_free - z[free]
        neg <- p < 0
        ratio <- -z[free][neg] / p[neg]
        alpha <- min(1, ratio)
        z[free] <- pmax(z[free] + alpha * p, 0)
        z[free][neg][ratio <= alpha] <- 0
        free <- free[z[free] > 0]
        next
      }
    }
    grad <- as.vector(h %*% z) + c
    bound <- setdiff(seq_len(n), free)
    if (!length(bound)) {
      break
    }
    j <- bound[which.min(grad[bound])]
    if (grad[j] >= -1e-15 * max(1, abs(c))) {
      break
    }
    free <- c(free, j)
  }
  z
}
