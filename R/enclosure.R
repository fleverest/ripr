#' @include region.R family.R
NULL

# The interface certify() bounds expectations through, and the branch and
# bound that drives it. A bounding method is an `enclosure` subclass: the
# driver owns the simplicial geometry (which edge to bisect, the node tree,
# pruning and stopping), and the method supplies upper bounds and attained
# values on each sub-simplex. `bounding_method()` decides which method, if any,
# covers a cell. The built-in methods are in `bounding_methods.R`.

# --- The interface ------------------------------------------------------------

#' Enclosure: bounding expectations over a cell
#'
#' [certify()] is, in essence, just an interface to `enclosure`. In `ripr`,
#' an `enclosure` is a method of bounding
#' \eqn{\sup_{\theta \in C} E_\theta[X]}{sup_{theta in C} E_theta[X]} over a
#' simplex \eqn{C}{C} of the null. [certify()] runs branch and bound over
#' the cells, keeping a tree of sub-simplices ("boxes"). At each iteration, it
#' splits the boxes with the largest upper bound by cutting its longest edge, and
#' prunes the boxes whose bound cannot beat the best value ("incumbent") attained
#' anywhere. The enclosure provides the interface to do a single step: it only
#' bounds and evaluates the incumbent on a single box, and provides the
#' interface for splitting.
#'
#' `enclosure` is abstract. A new bounding method subclasses it and implements:
#'
#' \describe{
#'   \item{`bounding_method(family, cell)`}{Registered on the pairs of
#'   [parametric_family] and [convex_region] the method covers. Returns an
#'   unprepared instance of the method's class, or a string saying why the
#'   cell cannot be bounded. The default simply raises an error.}
#'   \item{`prepare_enclosure(enclosure, x, family, control)`}{Precomputes
#'   whatever the method needs once per [certify()] call, shared by every cell
#'   that uses it (e.g. the Bernstein lattice for multinomial random variables,
#'   along with `x` evaluated on the entire sample space). Returns the prepared
#'   enclosure. The default returns it unchanged.}
#'   \item{`enclosure_root(enclosure, cell)`}{The first box for a cell: a list
#'   whose `V` is the vertex matrix, one vertex per row, plus anything else the
#'   method keeps per box.}
#'   \item{`enclosure_upper(enclosure, box)`}{An upper bound on
#'   \eqn{E_\theta[X]}{E_theta[X]} over the box.}
#'   \item{`enclosure_best(enclosure, box)`}{`list(value, theta)`: a value of
#'   \eqn{E_\theta[X]}{E_theta[X]} attained at a point `theta` of the box.}
#'   \item{`enclosure_split(enclosure, box, edge)`}{Returns the two boxes
#'   covering `box` resulting from cutting the edge between vertex rows
#'   `edge[1]` and `edge[2]`.}
#'   \item{`enclosure_edge(enclosure, box)`}{Optional: selects the edge to
#'   bisect. The default is to select the longest (in Euclidean distance)
#'   edge.}
#' }
#'
#' An effective enclosure should shrink with the size of the box. For instance,
#' Bernstein branch and bound converges at a rate of O(h^2), where h is the box
#' diameter.
#'
#' @param label Name of the method, reported in a certificate's `method`.
#' @return An `enclosure`.
#' @seealso [certify()]
#' @export
enclosure <- new_class(
  "enclosure",
  abstract = TRUE,
  properties = list(label = class_character)
)


#' @rdname enclosure
#' @param family A [parametric_family].
#' @param cell A [convex_region], one cell of a null.
#' @export
bounding_method <- new_generic(
  "bounding_method",
  c("family", "cell"),
  function(family, cell) S7::S7_dispatch()
)


#' @rdname enclosure
#' @param enclosure An `enclosure`: unprepared for `prepare_enclosure()`,
#'   prepared for the rest.
#' @param x A [random_variable].
#' @param control A list of options, e.g. `max_coefficients` for the enclosure
#'   method.
#' @export
prepare_enclosure <- new_generic(
  "prepare_enclosure",
  "enclosure",
  function(enclosure, x, family, control) S7::S7_dispatch()
)

method(prepare_enclosure, enclosure) <- function(
  enclosure,
  x,
  family,
  control
) {
  enclosure
}


#' @rdname enclosure
#' @export
enclosure_root <- new_generic(
  "enclosure_root",
  "enclosure",
  function(enclosure, cell) S7::S7_dispatch()
)


#' @rdname enclosure
#' @param box A box, as returned by `enclosure_root()` or `enclosure_split()`.
#' @export
enclosure_upper <- new_generic(
  "enclosure_upper",
  "enclosure",
  function(enclosure, box) S7::S7_dispatch()
)


#' @rdname enclosure
#' @export
enclosure_best <- new_generic(
  "enclosure_best",
  "enclosure",
  function(enclosure, box) S7::S7_dispatch()
)


#' @rdname enclosure
#' @param edge The two vertex rows of `box$V` whose edge is bisected.
#' @export
enclosure_split <- new_generic(
  "enclosure_split",
  "enclosure",
  function(enclosure, box, edge) S7::S7_dispatch()
)


#' @rdname enclosure
#' @export
enclosure_edge <- new_generic(
  "enclosure_edge",
  "enclosure",
  function(enclosure, box) S7::S7_dispatch()
)

method(enclosure_edge, enclosure) <- function(enclosure, box) {
  V <- box$V
  longest_edge(V, utils::combn(nrow(V), 2L))
}


#' Endpoints of the longest edge of vertex matrix `V`, among `edges` (2 x E)
#' @keywords internal
#' @noRd
longest_edge <- function(V, edges) {
  d2 <- rowSums(
    (V[edges[1L, ], , drop = FALSE] -
      V[edges[2L, ], , drop = FALSE])^2
  )
  edges[, which.max(d2)]
}


# --- Branch and bound ---------------------------------------------------------

#' A branch-and-bound node: a box plus its cached upper bound and lineage
#' @keywords internal
#' @noRd
node <- function(enc, box, id, parent = NA_integer_, depth = 0L, born = 0L) {
  box$ub <- enclosure_upper(enc, box)
  box$id <- id
  box$parent <- parent
  box$depth <- depth
  # With `parent` and the record's `retired`, makes the run replayable.
  box$born <- born
  box
}


#' A node without the method's data, for the history
#' @keywords internal
#' @noRd
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


#' @keywords internal
#' @noRd
node_ubs <- function(nodes) vapply(nodes, function(b) b$ub, numeric(1L))


#' The best value attained over several boxes
#' @keywords internal
#' @noRd
best_of <- function(enc, boxes) {
  best <- list(value = -Inf, theta = NULL)
  for (b in boxes) {
    cand <- enclosure_best(enc, b)
    if (cand$value > best$value) {
      best <- cand
    }
  }
  best
}


#' Why the search should stop, or `NULL`. Only "budget_hit" qualifies the
#' result.
#' @keywords internal
#' @noRd
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


#' Yields a global upper bound on `sup E_theta[X]` over a cell
#'
#' `bound` is valid (up to rounding) at every iteration. It includes the
#' largest bound ever pruned, since a pruned node may still hold the supremum.
#' `shared_incumbent` is a value attained elsewhere (e.g. another cell of the
#' same null), used only for pruning and stopping. Exactly one of `converged`
#' and `budget_hit` is `TRUE`.
#' @keywords internal
#' @noRd
branch_and_bound <- function(
  enc,
  box,
  tol = 1e-3,
  max_iter = 500L,
  shared_incumbent = -Inf
) {
  max_iter <- as.integer(max_iter)

  active <- list(node(enc, box, id = 1L))
  ubs <- node_ubs(active)
  next_id <- 2L
  best <- enclosure_best(enc, box)
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
    kids <- enclosure_split(enc, parent, enclosure_edge(enc, parent))
    kids <- lapply(seq_along(kids), function(i) {
      node(
        enc,
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
    kid_best <- best_of(enc, kids)
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
