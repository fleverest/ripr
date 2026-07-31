# The certified counterpart to `oracle()`. Where `oracle()` is a heuristic
# global search returning a *lower* bound on the face maximum, `oracle_bound()`
# returns a proven *upper* bound on `sup_theta G` over the face -- which is the
# quantity a duality-gap certificate actually needs.
#
# The seam is deliberately narrow: everything geometric and numerical lives in
# bernstein.R, everything about which problems qualify lives in
# `certifiable()`, and nothing here touches `oracle()` or the optimisation loop.

#' Tuning for the Bernstein branch-and-bound bound
#'
#' @param tol Stop subdividing once `bound - incumbent <= tol`. The bound is
#'   valid at every iteration, so this trades tightness for time, never
#'   validity.
#' @param max_iter Cap on bisections per face. Reaching it is not an error: the
#'   bound returned is still valid, merely looser.
#' @param max_coef Refuse to certify problems with more Bernstein coefficients
#'   than this, so an oversized batch is a planned fallback to the heuristic
#'   oracle rather than an out-of-memory kill. Coefficients per sub-simplex are
#'   `choose(n + K - 1, K - 1)`: 10,626 at `n = 20, K = 5`; 4.6M at
#'   `n = 100, K = 5`. Certification caps the usable batch size more tightly
#'   than optimisation does.
#' @param keep_argmax Prune with `U(S) < incumbent` rather than
#'   `U(S) <= incumbent`, retaining ties, so that the returned `active` set is a
#'   certified enclosure of every maximiser. Costs a little work and some
#'   memory; `active` is `NULL` without it.
#' @param round_slack Add `n * .Machine$double.eps * max(abs(coef))` to the
#'   reported bound. de Casteljau is all convex combinations and so is stable,
#'   but round-to-nearest can put the computed maximum coefficient marginally
#'   *below* the true one, which is the unsafe direction for a validity claim.
#' @return A named list of the above.
#' @export
bnb_control <- function(
  tol = 1e-3,
  max_iter = 500L,
  max_coef = 2e5,
  keep_argmax = FALSE,
  round_slack = TRUE
) {
  stopifnot(
    is.numeric(tol), length(tol) == 1L, is.finite(tol), tol >= 0,
    is.numeric(max_iter), length(max_iter) == 1L, max_iter >= 0,
    is.numeric(max_coef), length(max_coef) == 1L, max_coef >= 0,
    is.logical(keep_argmax), length(keep_argmax) == 1L, !is.na(keep_argmax),
    is.logical(round_slack), length(round_slack) == 1L, !is.na(round_slack)
  )
  list(
    tol = as.numeric(tol),
    max_iter = as.integer(max_iter),
    max_coef = as.numeric(max_coef),
    keep_argmax = keep_argmax,
    round_slack = round_slack
  )
}

# Is `vertices` a (possibly degenerate) sub-simplex of the probability simplex?
#
# Validity of the bound needs the barycentric weights to be non-negative and to
# sum to one -- that is what makes the Bernstein basis over the sub-simplex
# non-negative and a partition of unity. Non-degeneracy is *not* required: de
# Casteljau evaluates the blossom regardless, and the convex hull property holds
# for the hull of whatever vertices are supplied.
face_simplex_reason <- function(vertices, K) {
  if (!is.matrix(vertices) || !all(is.finite(vertices))) {
    return("its vertex matrix is not a finite numeric matrix")
  }
  if (nrow(vertices) != K) {
    return(sprintf(
      "its vertices have %d coordinates, not the family's %d",
      nrow(vertices), K
    ))
  }
  if (ncol(vertices) != K) {
    return(sprintf(
      paste0(
        "it has %d vertices, not %d: the Bernstein reparametrisation is ",
        "defined for simplices only"
      ),
      ncol(vertices), K
    ))
  }
  if (any(vertices < 0)) {
    return("it has a vertex outside the probability simplex (negative entry)")
  }
  if (max(abs(colSums(vertices) - 1)) > 1e-9) {
    return("it has a vertex whose coordinates do not sum to 1")
  }
  NULL
}

#' Why can this (face, engine, family) triple not be certified?
#'
#' Returns `NULL` when it can, and otherwise a character vector of
#' human-readable clauses for [certify()] to fold into its `require_certified`
#' error. Every applicable blocker is reported, not just the first: a user told
#' only that their engine is stochastic would swap it out and then discover
#' that their family has no Bernstein form either.
#' @keywords internal
#' @noRd
certifiable_reason <- function(face, engine, family, control = bnb_control()) {
  why <- character(0L)
  if (!deterministic(engine)) {
    why <- c(why, sprintf(
      "the %s engine is stochastic, so its `G` is an estimate rather than a polynomial",
      S7_class(engine)@name
    ))
  }
  bf <- bernstein_form(family)
  if (is.null(bf)) {
    why <- c(why, sprintf(
      "the %s family has no Bernstein form (its density is not a Bernstein basis function of theta)",
      S7_class(family)@name
    ))
  } else if (n_coefficients(bf) > control$max_coef) {
    why <- c(why, sprintf(
      "the Bernstein form needs %.0f coefficients, above `max_coef = %.0f`",
      n_coefficients(bf), control$max_coef
    ))
  }
  if (!S7_inherits(face, polytope_face)) {
    why <- c(why, sprintf(
      "face %s is a %s: unbounded geometries have no simplex to subdivide",
      face_label(face), S7_class(face)@name
    ))
  } else if (!is.null(bf)) {
    bad <- face_simplex_reason(face@vertices, bf$K)
    if (!is.null(bad)) {
      why <- c(why, sprintf("face %s cannot be used: %s", face_label(face), bad))
    }
  }
  if (length(why) == 0L) NULL else why
}

# A short label for a face in diagnostics: its own index when it has one.
face_label <- function(face) {
  fi <- tryCatch(face@face_index, error = function(e) NA_real_)
  if (length(fi) == 1L && !is.na(fi)) as.character(fi) else "<unlabelled>"
}

#' Can `sup_theta G` be bounded deterministically over this face?
#'
#' A plain predicate rather than a generic: it must be cheap, total, and
#' testable, and it has to consult three objects at once. Including the cost
#' budget lets callers plan rather than discover failure mid-run.
#'
#' All four conditions are necessary. Being exact is not sufficient -- the
#' Bernstein identity is a property of the *multinomial* pmf, so
#' [bernstein_form()] is consulted on the family, not the engine. The face must
#' be a bounded simplex inside the probability simplex, since that is what the
#' reparametrisation subdivides.
#'
#' @param face A `face`.
#' @param engine An `expectation_engine`.
#' @param family A `family`.
#' @param control A [bnb_control()] list, consulted for `max_coef`.
#' @return `TRUE` or `FALSE`.
#' @export
certifiable <- function(face, engine, family, control = bnb_control()) {
  is.null(certifiable_reason(face, engine, family, control))
}

#' Certified upper bound on `sup_theta G` over a face
#'
#' The upper-bound counterpart to [oracle()].
#'
#' Contract: returns `list(bound =, incumbent =, theta =, active =,
#' iterations =, method =)` where **`bound` is a valid upper bound on
#' `sup_theta G` over the face** -- unlike [oracle()], whose `value` is a lower
#' bound. So `incumbent <= sup G <= bound`. `theta` is the best point found, and
#' lies on the face. `active` is a list of sub-simplices whose union contains
#' every maximiser, or `NULL` when the run was not asked to track them
#' (`keep_argmax = FALSE`, the default).
#'
#' `bound` is valid at *every* iteration, not only at convergence: `tol` and
#' `max_iter` govern how tight it is, never whether it holds.
#'
#' # Floating point
#'
#' The validity claim is exact in real arithmetic. In floating point it holds to
#' a relative accuracy of order `1e-14`, and `round_slack` (see [bnb_control()])
#' is a calibrated pad of the right order rather than a proven error envelope.
#' The loss comes from dynamic range: de Casteljau averages coefficients whose
#' maximum can be far larger than `sup G` -- a mixture that gives some outcome
#' little mass makes `q(x) / p_W(x)` enormous there -- and averaging numbers
#' `R` times larger than the result costs about `log10(R)` digits. Measured on a
#' deliberately ill-conditioned case (`max(coef) / sup G` about 3500), the pad
#' absorbed all but 4e-15 of the relative error.
#'
#' The practical consequence is only that `bound` and an independent evaluation
#' of `G` -- [expect_ratio()]'s log-space sum, say -- agree to about `1e-14`
#' relative, not exactly, and either may be the larger. Nothing statistical
#' turns on it: a rescaling that is off by a factor of `1 + 1e-14` is not a
#' failure mode of an e-variable. Compare the two with a relative tolerance.
#'
#' Only `list(polytope_face, exact_engine)` has a method. `halfspace_face` is
#' deliberately left without one -- it is unbounded, so there is no simplex to
#' subdivide, and that is a permanent gap rather than a temporary one. Use
#' [certifiable()] to test before dispatching.
#'
#' @param face A `face`.
#' @param engine An `expectation_engine`.
#' @param log_Pw Length-M log-density of the candidate mixture over the engine's
#'   outcome set -- the same vector [oracle_step()] computes.
#' @param family The `family` supplying the [bernstein_form()].
#' @param control A [bnb_control()] list.
#' @param ... Ignored.
#' @return See the contract above.
#' @export
oracle_bound <- new_generic(
  "oracle_bound",
  c("face", "engine"),
  function(face, engine, log_Pw, family, control = bnb_control(), ...) {
    S7::S7_dispatch()
  }
)

method(oracle_bound, list(polytope_face, exact_engine)) <- function(
  face,
  engine,
  log_Pw,
  family,
  control = bnb_control(),
  ...
) {
  bf <- bernstein_form(family)
  if (is.null(bf)) {
    stop(
      "`oracle_bound()` needs a family with a Bernstein form; ",
      "test with `certifiable()` first"
    )
  }
  why <- face_simplex_reason(face@vertices, bf$K)
  if (!is.null(why)) {
    stop(sprintf("this face cannot be certified: %s", why), call. = FALSE)
  }

  lat <- bernstein_lattice(bf$n, bf$K, tally = engine@outcomes)

  # The whole objective, with no re-derivation: G(theta) = sum_x (q(x)/p_W(x))
  # P_theta(x), and P_theta(x) is the Bernstein basis function at x. NaN comes
  # only from -Inf - -Inf, i.e. an outcome carrying no Q-mass under a mixture
  # that also gives it none; it contributes nothing, matching `expect_ratio()`.
  coef <- exp(nan_to_neginf(engine@log_q_mass - log_Pw))
  if (!all(is.finite(coef))) {
    # Not defensive padding: a non-finite coefficient means Q puts mass where
    # P_W has none, so the candidate is not an e-variable at *any* rescaling.
    stop(
      "the candidate mixture has no support where Q has mass: ",
      "`G` is unbounded and no rescaling makes this an e-variable",
      call. = FALSE
    )
  }

  # Divide out any residual drift in the column sums so the barycentric weights
  # sum to exactly one; `face_simplex_reason()` has already bounded the drift.
  V <- sweep(face@vertices, 2L, colSums(face@vertices), "/")
  box <- list(V = V, coef = reparametrise_to(coef, lat, V))

  res <- certify_sup(
    list(box),
    lat,
    tol = control$tol,
    max_iter = control$max_iter,
    keep_argmax = control$keep_argmax,
    round_slack = control$round_slack
  )
  list(
    bound = res$bound,
    incumbent = res$incumbent,
    theta = res$theta,
    # Without `keep_argmax` the pruning rule discards boxes that may attain the
    # maximum, so the active set is not an enclosure and must not be presented
    # as one.
    active = if (control$keep_argmax) res$active else NULL,
    iterations = res$iterations,
    method = "bernstein_bnb"
  )
}
