# The headline guarantee: the certified, rescaled e-variable has null
# expectation at most 1 everywhere on the null hypothesis.

test_that("certified e-variable satisfies E_{P_theta}[e] <= 1 on the null", {
  n <- 8
  alt <- point_dist(theta_star = c(0.6, 0.25, 0.15))
  fam <- multinomial_family(n, 3)
  faces <- plurality_faces(3)
  prob <- ripr_problem(fam, faces, alt)

  set.seed(30)
  res <- run_ripr(
    prob,
    init_atoms = init_on_faces(faces, alt@theta_star),
    init_atom_faces = seq_along(faces),
    fw_iters = 20, em_iters = 5, n_seeds = 120, verbose = FALSE
  )

  X <- support(fam)
  # The rescaled e-variable e(x) = (Q / P*)(x) / (1 + gap), straight from the
  # fitted e_variable object.
  e <- e_value(res$e_variable, X)

  E_e <- function(theta) sum(exp(log_density(fam, theta)) * e)

  # Boundary null points: theta_1 = theta_j and candidate j dominates.
  boundary <- list(
    c(0.5, 0.5, 0.0),
    c(0.45, 0.45, 0.10),
    c(0.4, 0.4, 0.2),
    c(1 / 3, 1 / 3, 1 / 3)
  )
  # Interior null points: candidate 1 strictly trails.
  interior <- list(
    c(0.2, 0.5, 0.3),
    c(0.1, 0.6, 0.3),
    c(0.25, 0.4, 0.35),
    c(0.0, 0.7, 0.3)
  )
  for (theta in c(boundary, interior)) {
    expect_lte(E_e(theta), 1 + 1e-6)
  }
})
