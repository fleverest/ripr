
<!-- README.md is generated from README.Rmd. Edit the .Rmd and run devtools::build_readme(). -->

# ripr

`ripr` computes the **reverse information projection** (RIPr) `P*` of an
alternative distribution `Q` onto a (possibly non-convex,
union-structured) null hypothesis, together with the **duality-gap
certificate** for the rescaled e-variable `(Q / P*) / (1 + gap)`, whose
expectation under every null distribution is at most 1.

The projection is found by Frank–Wolfe and EM over a mixture of atoms on
the null faces. It requires the specification of a sampling `family`
(the model `p_theta`), the null geometry (a union of `face`s), and a
surrogate `alternative` `Q`. The alternative is an
`outcome_distribution` — a law over the sample space — which you can
supply directly or build from a `mixing` measure over the parameter
space (e.g. `point_mixing`, `finite_mixing`), marginalised through the
family. These components are templated as
[S7](https://rconsortium.github.io/S7/) base classes, so new sample
spaces drop in without touching the algorithmic core. The algorithms
require the ability to compute KL divergences and gradients, which may
be implemented either via **exact** numerical methods or by **Monte
Carlo** sampling. These two routes are encapsulated in the `engine`
class.

``` r
# install from github
#renv::install("fleverest/ripr")
```

## Example 1: a one-sided binomial null

Let `X ~ Binomial(10, p)`, written as a 2-category multinomial with
`theta = (p, 1 - p)`. Consider testing `H_0: p <= 1/2` against
`H_1: p > 1/2`, adopting the simple alternative `Q = Binomial(10, 3/4)`.
The RIPr of `Q` onto the convex null is a **point mass at `p = 1/2`**.
We demonstrate this using both the exact and Monte Carlo engines.

``` r
#library(ripr)
devtools::load_all()
n <- 10
fam <- multinomial_family(n_trials = n, k = 2)
# `ripr_problem` takes Q as an `outcome_distribution` over the sample space.
# Build one by marginalising a `mixing` measure over the parameter theta (here
# it is just a point mass at theta = (3/4, 1/4)) through the family, giving the
# outcome law Q = Binomial(10, 3/4).
Q <- as_marginal(point_mixing(theta_star = c(0.75, 0.25)), fam)

# H_0 = {p <= 1/2} is the simplex segment from p = 0 to the tie p = 1/2.
face <- polytope_face(vertices = cbind(c(0, 1), c(0.5, 0.5)), face_index = 1)
null <- null_region(faces = list(face))
```

**Exact engine** (enumerates all 11 outcomes):

``` r
set.seed(1)
prob <- ripr_problem(fam, null, Q)
res <- run_ripr(
  prob,
  # Deliberately initialising somewhere sub-optimal for illustrative purposes
  init_atoms = matrix(c(0.25, 0.75), ncol = 1),
  init_atom_faces = 1L,
  fw_iters = 8,
  em_iters = 5,
  prune_threshold = 1e-6, # drop the numerically-dead atoms for output
  verbose = FALSE
)

cbind(
  p = res$projection@mixing@components[1, ],
  weight = res$projection@mixing@weights
)
##              p weight
## [1,] 0.4999973      1
```

`res$projection` is a `marginal` — the fitted `P*` as a distribution
over outcomes — and its `@mixing` is the `finite_mixing` holding the
atoms and weights on the null. `run_ripr()` returns the fitted
e-variable directly; call `e_value()` on new data. Observing `(8, 2)`
yields evidence against `H_0`, while a tie `(5, 5)` yields evidence in
favour of it:

``` r
res$e_variable
## <e_variable>  e(x) = (Q / P*) / (1 + gap)
##   numerator  Q  : point_mixing
##   projection P* : finite_mixing (1 atom)
##   gap           : 0   (correction 1 + gap = 1)
e_value(res$e_variable, rbind(c(8, 2), c(5, 5)))
## [1] 6.4074351 0.2373047
```

**Monte Carlo engine** (sample 4000 draws from `Q`) reaches the same
projection:

``` r
eng_mc <- mc_engine(fam, Q, n_draws = 1000)
res_mc <- run_ripr(
  ripr_problem(fam, null, Q, engine = eng_mc),
  init_atoms = matrix(c(0.25, 0.75), ncol = 1),
  init_atom_faces = 1L,
  fw_iters = 10,
  em_iters = 10,
  gap_tol = 1e-3,
  verbose = FALSE
)

cbind(
  p = res_mc$projection@mixing@components[1, ],
  weight = res_mc$projection@mixing@weights
)
##              p weight
## [1,] 0.4999897      1
```

## Example 2: a 2-D Gaussian cone null (union of half-spaces)

Now let `X ~ N(mu, I_2)` and test the **cone alternative**
`H_1: mu_1 > |mu_2|` — equivalently `mu_1 > mu_2` **and** `mu_1 > -mu_2`
— with the point alternative `mu = (2, 1)`. Its complement, the null
`H_0 = {mu_1 <= mu_2} ∪ {mu_1 <= -mu_2}`, is a union of two half-spaces.

``` r
fam <- gaussian_family(dim = 2) # standard: sigma = I
mu <- c(2, 1)
Q <- as_marginal(point_mixing(theta_star = mu), fam)

faces <- list(
  halfspace_face(v = c(1, -1), c = 0, face_index = 1), # mu_1 <= mu_2
  halfspace_face(v = c(1, 1), c = 0, face_index = 2) # mu_1 <= -mu_2
)
null <- null_region(faces = faces)
```

The sample space is continuous here, so we use the Monte Carlo engine:
1000 draws to fit the mixture, and fresh `certify_draws = 5000` samples,
independent of the fit, to certify the result (drawn twice by default —
see below).

``` r
set.seed(2)
eng <- mc_engine(fam, Q, n_draws = 1000)
prob <- ripr_problem(fam, null, Q, engine = eng)
init <- cbind(
  init_point(faces[[1]], mu),
  init_point(faces[[2]], mu)
)

res <- run_ripr(
  prob,
  init_atoms = init,
  init_atom_faces = c(1, 2),
  fw_iters = 25,
  em_iters = 5,
  n_seeds = 50,
  gap_tol = 1e-3,
  certify_draws = 5000,
  verbose = FALSE
)

M <- res$projection@mixing@components
cbind(
  mu1 = M[1, ],
  mu2 = M[2, ],
  weight = round(res$projection@mixing@weights, 3)
)
##           mu1       mu2 weight
## [1,] 1.591870  1.591870  0.951
## [2,] 1.168331 -1.168331  0.048
## [3,] 2.751958 -2.752034  0.001
## [4,] 2.779545 -2.779711  0.000
```

The projection concentrates on the `mu_1 = mu_2` tie-point — the null
boundary nearest the alternative `mu = (2, 1)` — with a little mass on
the `mu_1 = -mu_2` face. The certificate reports the guaranteed e-value
growth rate and, for a Monte Carlo engine, the standard error of the
estimated gap:

``` r
c(
  gap = res$gap,
  gap_se = res$certificate$gap_se,
  growth_rate = res$certificate$growth_rate
)
##         gap      gap_se growth_rate 
##  0.10495159  0.04916562  0.05122194
```

### Re-certifying on a larger sample

The Monte Carlo certificate above is a stochastic estimate: `certify()`
resamples fresh draws from `Q` (it does *not* reuse the fit draws) and,
by default, split-samples — one sample to locate the worst-case
direction `theta*`, a second to estimate the gap at it, so the estimate
carries no selection bias. Its standard error shrinks like
`1 / sqrt(n)`, so re-certifying the *same* fitted projection on more
draws tightens it, with no re-fitting — just raise `n_draws`:

``` r
set.seed(9)
cert <- certify(res$projection, prob, n_draws = 50000)
c(new_gap = cert$gap, new_se = cert$gap_se)
##    new_gap     new_se 
## 0.08703831 0.01361741
```

`certify()` takes the fitted `projection` (or any `finite_mixing` of
atoms on the null) and sweeps the oracle over the null faces. It
resamples the problem’s own engine by default; pass `engine =` to
certify against a different one (say an `exact_engine`), and
`split = FALSE` for a single-sample estimate.
