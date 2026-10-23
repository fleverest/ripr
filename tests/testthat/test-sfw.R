# Properties of sfw_step().
#
# Each step is a fully-corrective FW step followed by a slide, so it must
# record one row per step, never end above the KL a plain corrective step
# would reach from the same state, keep atoms in their own parts, and reach a
# small gap quickly on a problem fw_step() alone converges on slowly.

sfw_fixture <- function(n = 8L, q = c(0.4, 0.35, 0.25)) {
  fam <- multinomial_family(n_trials = n, k = 3L)
  parts <- lapply(2:3, function(j) {
    vertices <- diag(3)
    vertices[1L, ] <- replace(numeric(3), c(1L, j), 0.5)
    simplex_region(vertices = vertices)
  })
  ripr_init(
    fam(q),
    null_model(fam, parts),
    exact_engine(),
    control = ripr_control(n_seeds = 30L, n_restarts = 4L)
  )
}

test_that("sfw_step records one sfw row per step, with exact KL", {
  set.seed(1)
  st <- sfw_step(sfw_fixture(), times = 3L)
  expect_identical(st@trace$phase, c("init", "sfw", "sfw", "sfw"))
  expect_equal(utils::tail(st@trace$kl, 1L), exact_kl(st), tolerance = 1e-10)
  expect_equal(sum(weights(st@mixing)), 1)
  expect_true(in_own_part(st))
})

test_that("the slide never leaves KL above the corrective step alone", {
  set.seed(1)
  st <- fw_step(sfw_fixture(), times = 2L, correct = TRUE)
  set.seed(2)
  fw <- fw_step(st, correct = TRUE)
  set.seed(2)
  sfw <- sfw_step(st)
  expect_lte(
    utils::tail(sfw@trace$kl, 1L),
    utils::tail(fw@trace$kl, 1L) + 1e-12
  )
})

test_that("sfw_step reaches a small gap in few steps", {
  set.seed(1)
  st <- sfw_step(sfw_fixture(), times = 15L, until = gap_below(1e-8))
  expect_lt(st@oracle$value - 1, 1e-8)
  expect_lte(sum(st@trace$phase == "sfw"), 10L)
  expect_lte(n_atoms(st@mixing), 4L)
})

test_that("sfw_step honours until and checks max_eval", {
  st <- sfw_fixture()
  none <- sfw_step(st, times = 5L, until = function(state, candidate) TRUE)
  expect_identical(nrow(none@trace), nrow(st@trace))
  expect_error(sfw_step(st, max_eval = 0))
  expect_no_error(sfw_step(st, max_eval = 5L))
})

test_that("merge_coincident folds same-part duplicates only", {
  atoms <- rbind(
    c(0.2, 0.3, 0.5),
    c(0.2, 0.3, 0.5),
    c(0.2, 0.3, 0.5),
    c(0.6, 0.2, 0.2)
  )
  m <- ripr:::merge_coincident(
    atoms,
    c(0.1, 0.2, 0.3, 0.4),
    part = c(1L, 1L, 2L, 1L)
  )
  expect_identical(nrow(m$atoms), 3L)
  expect_equal(m$weights, c(0.3, 0.3, 0.4))
  expect_identical(m$part, c(1L, 2L, 1L))
})


test_that("sfw_step slides a state that came from another verb before inserting", {
  set.seed(1)
  st <- fw_step(sfw_fixture(), times = 2L)
  slides <- 0L
  slide <- ripr:::slide_state
  local_mocked_bindings(
    slide_state = function(...) {
      slides <<- slides + 1L
      slide(...)
    },
    .package = "ripr"
  )
  stepped <- sfw_step(st)
  expect_identical(slides, 2L) # once before inserting, once after
  slides <- 0L
  sfw_step(stepped)
  expect_identical(slides, 1L) # a state from sfw_step() is not slid first
})

test_that("merge_coincident merges at the weight-averaged position, within tolerance", {
  a <- c(0.2, 0.3, 0.5)
  b <- a + c(1e-10, -1e-10, 0)
  m <- ripr:::merge_coincident(rbind(a, b), c(0.25, 0.75), part = c(1L, 1L))
  expect_identical(nrow(m$atoms), 1L)
  expect_equal(as.vector(m$atoms), 0.25 * a + 0.75 * b, tolerance = 0)
  expect_equal(m$weights, 1)
  # Beyond atol + rtol * max(|a|, |b|), nothing merges.
  far <- ripr:::merge_coincident(rbind(a, a + c(1e-6, -1e-6, 0)), c(0.5, 0.5), part = c(1L, 1L))
  expect_identical(nrow(far$atoms), 2L)
  # rtol scales with the coordinates' size.
  big <- c(1e6, 0)
  near <- rbind(big, big + c(1e-4, 0))
  expect_identical(nrow(ripr:::merge_coincident(near, c(0.5, 0.5), c(1L, 1L), atol = 1e-9, rtol = 0)$atoms), 2L)
  expect_identical(nrow(ripr:::merge_coincident(near, c(0.5, 0.5), c(1L, 1L), atol = 1e-9, rtol = 1e-9)$atoms), 1L)
})

test_that("ripr_control checks the merge tolerances", {
  expect_identical(ripr_control()$merge_atol, 1e-9)
  expect_identical(ripr_control(merge_rtol = 0)$merge_rtol, 0)
  expect_error(ripr_control(merge_atol = -1), "larger than or equal to 0")
})
