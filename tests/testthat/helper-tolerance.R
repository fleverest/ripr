# Certification is mathematically exact, but floating point violates a tiny bit.
# This just defines a consistent relative tolerance for unit tests.
rounding_tol <- function(x) 1e-12 * max(1, abs(x))
