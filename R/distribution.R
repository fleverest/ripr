#' Probability distribution (abstract root)
#'
#' The umbrella for every probability law the package handles. It has two
#' branches: a [mixing] measure over a family's *parameter* space Theta, and an
#' [outcome_distribution] over its *sample* space X. The two are bridged by
#' [as_marginal()], which pushes a mixing through a family's kernel to obtain the
#' induced outcome distribution (a [marginal]).
#'
#' This root carries no interface of its own -- it exists so the two branches
#' share a common type. The evaluable interface (`dist_log_density`,
#' `dist_sample`) lives on [outcome_distribution], the sample-space branch the
#' RIPr problem actually consumes.
#' @export
distribution <- new_class("distribution", abstract = TRUE)
