#' Bundle a RIPr problem: family + null + alternative + engine
#'
#' The problem object every algorithm step consumes. Normalises each component
#' and builds an exact engine by default when none is supplied.
#'
#' @param family Sampling model, in any form accepted by [as_family()].
#' @param null Null-hypothesis faces, in any form accepted by [as_faces()].
#' @param alternative The numerator Q as an [outcome_distribution] over the
#'   sample space (e.g. `as_marginal(point_mixing(...), family)`).
#' @param engine Optional `expectation_engine`; defaults to
#'   `exact_engine(family, alternative)`.
#' @return List with components `family`, `null` (list of faces), `alternative`,
#'   `engine`.
#' @export
ripr_problem <- function(family, null, alternative, engine = NULL) {
  family <- as_family(family)
  null <- as_faces(null)
  alternative <- as_outcome_distribution(alternative)
  if (is.null(engine)) {
    engine <- exact_engine(family = family, alternative = alternative)
  }
  list(
    family = family,
    null = null,
    alternative = alternative,
    engine = engine
  )
}
