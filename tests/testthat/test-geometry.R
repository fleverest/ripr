# Face geometry contracts: projection lands in the face and is idempotent,
# contains() agrees with project() fixed points, and the oracle returns a point
# on the face.

test_that("polytope face projection lands in the face and is idempotent", {
  faces <- plurality_faces(4)
  set.seed(20)
  for (f in faces) {
    for (rep in 1:3) {
      theta <- as.vector(rmultinom(1, 30, rep(1 / 4, 4))) / 30
      p <- project(f, theta)
      expect_true(contains(f, p, tol = 1e-6))
      # Idempotent: projecting an on-face point returns it.
      expect_equal(project(f, p), p, tolerance = 1e-6)
    }
  }
})

test_that("contains() agrees with project() fixed points", {
  # Boundary face for K = 3: theta_1 = theta_2 with candidate 2 dominating.
  f <- polytope_face(
    vertices = matrix(c(0.5, 0.5, 0, 1 / 3, 1 / 3, 1 / 3), nrow = 3),
    face_index = 2
  )
  on_face <- project(f, c(0.5, 0.5, 0))
  expect_true(contains(f, on_face))

  off_face <- c(0.8, 0.1, 0.1) # candidate 1 strictly leads: not on the tie face
  expect_false(contains(f, off_face, tol = 1e-6))
})

test_that("polytope oracle returns a point on the face", {
  n <- 6
  alt <- small_alt()
  fam <- multinomial_family(n, 3)
  eng <- exact_engine(fam, as_marginal(alt, fam))
  M <- nrow(support(fam))
  obj <- ripr:::fw_objective(eng, fam, log(rep(1 / M, M)))

  f <- plurality_faces(3)[[1L]]
  set.seed(21)
  res <- oracle(f, obj, n_seeds = 60L)
  expect_true(contains(f, res$theta, tol = 1e-6))
  expect_length(res$theta, 3L)
})

test_that("halfspace face projection and membership are closed form", {
  f <- halfspace_face(v = c(1, -1), c = 0, face_index = 1) # {theta_1 <= theta_2}
  inside <- c(0.2, 0.8)
  expect_true(contains(f, inside))
  expect_equal(project(f, inside), inside)

  outside <- c(0.9, 0.1)
  p <- project(f, outside)
  expect_true(contains(f, p, tol = 1e-8))
  expect_equal(project(f, p), p, tolerance = 1e-8)
})
