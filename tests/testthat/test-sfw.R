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
  st <- fw_step(sfw_fixture(), times = 2L, correct_weights = TRUE)
  set.seed(2)
  fw <- fw_step(st, correct_weights = TRUE)
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

test_that("sfw_step honours until and a tiny slide budget", {
  st <- sfw_fixture()
  none <- sfw_step(st, times = 5L, until = function(state, candidate) TRUE)
  expect_identical(nrow(none@trace), nrow(st@trace)) # No further iterations
  set.seed(1)
  long <- sfw_step(st) # No max_eval
  set.seed(1)
  st@control$slide_max_eval <- 2L
  expect_warning(short <- sfw_step(st), class = "ripr_slide_warning")
  expect_lte(
    utils::tail(long@trace$kl, 1L),
    utils::tail(short@trace$kl, 1L)
  )
})

test_that("the searches use the control's tolerance and budget", {
  st <- sfw_fixture()
  st@control$search_xtol <- 1e-5
  st@control$search_max_eval <- 7L
  seen <- list()
  maximise <- ripr:::maximise_over
  local_mocked_bindings(
    maximise_over = function(..., xtol, max_eval) {
      seen[[length(seen) + 1L]] <<- c(xtol, max_eval)
      maximise(..., xtol = xtol, max_eval = max_eval)
    },
    .package = "ripr"
  )
  set.seed(1)
  fw_step(st)
  em_step(st)
  expect_gt(length(seen), 0L)
  for (s in seen) {
    expect_identical(s, c(1e-5, 7))
  }
})

test_that("near_tol decides which maxima are one", {
  # Two end points 1e-8 apart are one maximum at the default tolerance and
  # two at a zero one. The search output is synthetic: insert_candidates()
  # reads only the end points' positions, their parts, and the order of their
  # values (all above one).
  st <- ripr_init(
    sfw_fixture()@alternative,
    sfw_fixture()@null,
    exact_engine(),
    atoms = list(rbind(c(0.3, 0.4, 0.3)), rbind(c(0.3, 0.3, 0.4)))
  )
  p <- c(0.1, 0.5, 0.4)
  found <- list(
    theta = c(0.1, 0.6, 0.3),
    value = 1.3,
    part = 1L,
    maxima = list(
      theta = rbind(p, p + c(0, 1e-8, -1e-8)),
      value = c(1.2, 1.1),
      part = c(1L, 1L)
    )
  ) # found$maxima$theta has two rows that are close
  expect_length(ripr:::insert_candidates(st, found), 2L)
  st@control$near_tol <- 0
  expect_length(ripr:::insert_candidates(st, found), 3L)
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
  # Default tolerances
  atol <- ripr_control()$merge_atol
  rtol <- ripr_control()$merge_rtol
  # Construct a, b
  a <- c(0.2, 0.3, 0.5)
  delta <- atol + rtol * norm(a, type = "2")
  b <- a + c(delta / 2, -delta / 2, 0) # closer than atol + rtol * max(|a|, |b|)
  # Check that they merge
  m <- ripr:::merge_coincident(rbind(a, b), c(0.25, 0.75), part = c(1L, 1L))
  expect_identical(nrow(m$atoms), 1L)
  expect_equal(as.vector(m$atoms), 0.25 * a + 0.75 * b, tolerance = 0)
  expect_equal(m$weights, 1)

  # Beyond atol + rtol * max(|a|, |b|), nothing merges.
  b <- a + c(delta, -delta, 0) # further than atol + rtol * max(|a|, |b|)
  m <- ripr:::merge_coincident(rbind(a, b), c(0.25, 0.75), part = c(1L, 1L))
  expect_identical(nrow(m$atoms), 2L)
  # rtol scales with the coordinates' size.
  big <- c(1e6, 0)
  delta <- atol + rtol * 1e6
  near <- rbind(big, big + c(delta, 0))
  expect_identical(
    nrow(
      ripr:::merge_coincident(
        near,
        c(0.5, 0.5),
        c(1L, 1L),
        atol = 1e-9,
        rtol = 0
      )$atoms
    ),
    2L
  )
  expect_identical(
    nrow(
      ripr:::merge_coincident(
        near,
        c(0.5, 0.5),
        c(1L, 1L),
        atol = 1e-9,
        rtol = 1e-9
      )$atoms
    ),
    1L
  )
})

test_that("ripr_control checks the merge tolerances", {
  expect_identical(ripr_control()$merge_atol, 1e-9)
  expect_identical(ripr_control(merge_rtol = 0)$merge_rtol, 0)
  expect_error(ripr_control(merge_atol = -1), "larger than or equal to 0")
})

test_that("insert_candidates merges a maximum found from two parts into one", {
  # (0, 0.5, 0.5) lies on the face both parts share, so each part's search
  # can end there; it must be inserted once, in the part with the higher G.
  # The search output is synthetic, as in the near_tol test above.
  st <- ripr_init(
    sfw_fixture()@alternative,
    sfw_fixture()@null,
    exact_engine(),
    atoms = list(rbind(c(0.3, 0.4, 0.3)), rbind(c(0.3, 0.3, 0.4)))
  )
  shared <- c(0, 0.5, 0.5)
  found <- list(
    theta = c(0.1, 0.6, 0.3),
    value = 1.3,
    part = 1L,
    maxima = list(
      theta = rbind(shared, shared + c(0, 1e-8, -1e-8), c(0.3, 0.4, 0.3)),
      value = c(1.1, 1.2, 1.05),
      part = c(1L, 2L, 1L)
    )
  )
  cands <- ripr:::insert_candidates(st, found)
  expect_length(cands, 2L)
  expect_equal(cands[[1L]]$theta, found$theta)
  expect_equal(cands[[2L]]$theta, shared, tolerance = 1e-7)
  expect_identical(cands[[2L]]$part, 2L)
  expect_true(contains(parts(st@null@region)[[2L]], cands[[2L]]$theta))
  # merge_coincident() itself still keeps parts apart unless told otherwise.
  m <- ripr:::merge_coincident(rbind(shared, shared), c(1, 1), c(1L, 2L))
  expect_identical(nrow(m$atoms), 2L)
})

test_that("the slide merges nearly coincident pairs it would leave apart", {
  # Start with a bunch of random atoms and run sliding Frank--Wolfe;
  # we want to make sure that atoms are merged appropriately.
  near_pairs <- function(st, tol) {
    a <- atoms(st@mixing)
    same <- outer(st@part, st@part, "==") & upper.tri(diag(nrow(a)))
    d <- as.matrix(stats::dist(a, method = "maximum"))[same]
    sum(d > 1e-9 & d <= tol)
  }
  stalled <- 0L
  for (seed in 1:6) {
    set.seed(seed)
    st <- ripr_init(
      sfw_fixture()@alternative,
      sfw_fixture()@null,
      exact_engine(),
      atoms = lapply(parts(sfw_fixture()@null@region), function(p) {
        ch <- chart(p)
        ch$to_theta_batch(ch$seed(40L))
      })
    )
    ld <- ripr:::compile_engine(st@engine)
    dg <- ripr:::compile_engine_grad(st@engine)
    apart <- ripr:::slide_state(st, ld, dg, tol = 0)$state
    merged <- ripr:::slide_state(st, ld, dg)$state
    stalled <- stalled + near_pairs(apart, st@control$near_tol)
    expect_identical(near_pairs(merged, st@control$near_tol), 0L)
    expect_lte(exact_kl(merged), exact_kl(apart) + 1e-12)
  }
  skip_if(stalled == 0L, "no slide stalled with a pair apart on this platform")
  expect_gt(stalled, 0L)
})

test_that("merge_coincident merges to a fixed point", {
  # The second and third points merge at their midpoint, which lands within
  # tolerance of the first: one pass would leave that pair, so passes repeat.
  tol <- 1e-3
  pts <- rbind(c(0, 0), c(1.2, 0.4) * tol, c(0.4, 1.2) * tol)
  m <- ripr:::merge_coincident(
    pts,
    c(1, 1, 1),
    c(1L, 1L, 1L),
    atol = tol,
    rtol = 0
  )
  expect_identical(nrow(m$atoms), 1L)
  expect_equal(m$weights, 3)
})

test_that("the slide's merge never costs KL", {
  # With a tolerance wide enough to catch genuinely distinct atoms, merging
  # them would raise KL, so the slide keeps them apart.
  fam <- multinomial_family(n_trials = 8L, k = 3L)
  medial <- simplex_region(
    vertices = rbind(c(0.5, 0.5, 0), c(0.5, 0, 0.5), c(0, 0.5, 0.5))
  )
  q <- fam(finite_dist(
    atoms = rbind(c(0.7, 0.15, 0.15), c(0.15, 0.7, 0.15)),
    weights = c(0.5, 0.5)
  ))
  set.seed(1)
  st <- sfw_step(
    ripr_init(q, null_model(fam, medial), exact_engine()),
    times = 10L,
    until = gap_below(1e-9)
  )
  expect_gt(n_atoms(st@mixing), 1L)
  ld <- ripr:::compile_engine(st@engine)
  dg <- ripr:::compile_engine_grad(st@engine)
  wide <- ripr:::slide_state(st, ld, dg, tol = 1)$state
  expect_identical(n_atoms(wide@mixing), n_atoms(st@mixing))
  expect_lte(exact_kl(wide), exact_kl(st) + 1e-12)
})
