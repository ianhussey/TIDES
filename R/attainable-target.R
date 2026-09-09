# Targeted certification for brimmest(): reduce ONE reported tuple to the
# handful of integer states it could have come from, without enumerating the
# whole attainable lattice.
#
# .attainable_lattice() answers "which (mean, sd) pairs are reachable?" and
# then brimmest() asks whether the report is among them. For a single report
# that is a membership question answered by a reachability-of-everything
# computation, and on a wide scale the full lattice is not merely slow but
# refused outright by its own size guard: a 0-63 inventory at n = 50 needs a
# 3151 x 24801 state table.
#
# The report pins both coordinates almost completely, which is what this file
# exploits:
#
#   * the mean's rounding interval admits only integer sums S in
#     [n*(m_lo - l), n*(m_hi - l)] -- usually one to a few values, not
#     thousands (this is the GRIM condition, read as a constraint rather than
#     a test);
#   * each candidate S then pins the sum of squares to a narrow window,
#     because Q = S^2/n + (n - 1) * sd^2 and sd is known to its reporting
#     precision.
#
# .target_states() below builds those (S, Q-window) targets, inverse-rounding
# the reported values once rather than forward-rounding the lattice cell by
# cell for every rounding rule. R/certify-sandwich.R then decides whether n
# integers in [0, W] can hit any of them.
#
# The rest of the file is the routing arithmetic brimmest() needs to choose
# between this route and the full lattice before committing to either.
#
# The verdict is identical to lattice membership; only the work differs. That
# equivalence is asserted cell for cell in tests/testthat/test-brimmest.R.

# Internal: how many cells would the full lattice need? Mirrors the arithmetic
# in .attainable_lattice() without allocating anything, so brimmest() can
# choose a route before committing to one. NA when the design is not on the
# integer grid at all.
.lattice_cells <- function(l, u, n, mg) {
  W <- mg * (u - l)
  if (abs(W - round(W)) > 1e-9) {
    return(NA_real_)
  }
  W <- as.integer(round(W))
  if (W < 1L) {
    return(NA_real_)
  }
  ys <- 0:W
  pr <- ys * (W - ys)
  g <- Reduce(.gcd2, pr[pr > 0])
  if (!length(g) || is.na(g) || g < 1) {
    g <- 1
  }
  (n * W + 1) * (n * max(pr) / g + 1)
}

# Internal: the full lattice, memoised per design. Enumerating it is the
# dominant cost of the lattice route, and repeated calls on one design are
# the common case when screening a table of reports.
.lattice_cache <- new.env(parent = emptyenv())

.lattice_cached <- function(l, u, n, mg, max_cells = 2e7) {
  key <- paste(l, u, n, mg, sep = "\r")
  hit <- get0(key, envir = .lattice_cache, inherits = FALSE)
  if (!is.null(hit)) {
    return(hit)
  }
  res <- .attainable_lattice(l, u, n, mg, max_cells = max_cells)
  assign(key, res, envir = .lattice_cache)
  res
}

# Internal: map the forward rounding vocabulary used by the lattice
# (.round_reported) onto the inverse-rounding vocabulary of unround_interval().
.rounding_inverse <- function(rounding) {
  switch(
    rounding,
    half_up = "up",
    half_down = "down",
    native = "even",
    ceiling = "ceiling",
    floor = "floor",
    trunc = "trunc",
    anti_trunc = "anti_trunc",
    stop("unknown rounding rule: ", rounding)
  )
}

# Internal: the (S, Q) targets implied by one reported (mean, sd) under one
# rounding rule. Returns NULL when the report admits no integer sum at all
# (the GRIM condition failing), which is already a certificate of
# impossibility. `mg` is the granularity multiplier from .scoring_geometry().
.target_states <- function(
  l,
  u,
  n,
  mg,
  mean,
  sd,
  mean_digits,
  sd_digits,
  rounding
) {
  inv <- .rounding_inverse(rounding)
  W <- as.integer(round(mg * (u - l)))
  iv_m <- unround_interval(mean, mean_digits, inv)
  iv_s <- unround_interval(sd, sd_digits, inv)

  # S = n * mg * (mean - l) must be a whole number of grid steps. Endpoint
  # inclusion matters: a half-up report owns its lower endpoint but not its
  # upper, and an S sitting exactly on an excluded endpoint is not admissible.
  s_lo <- n * mg * (iv_m$lo - l)
  s_hi <- n * mg * (iv_m$hi - l)
  s_from <- as.integer(ceiling(s_lo - 1e-9))
  s_to <- as.integer(floor(s_hi + 1e-9))
  # seq.int() counts DOWN when from > to, so an empty range has to be caught
  # explicitly; otherwise a mean admitting no integer sum yields two phantoms
  if (s_to < s_from) {
    return(NULL)
  }
  S <- seq.int(s_from, s_to)
  if (!iv_m$lo_incl) {
    S <- S[S > s_lo + 1e-9]
  }
  if (!iv_m$hi_incl) {
    S <- S[S < s_hi - 1e-9]
  }
  S <- S[S >= 0L & S <= n * W]
  if (!length(S)) {
    return(NULL)
  }

  # Q = S^2/n + (n - 1) * (mg * sd)^2, and the SD interval carries its own
  # endpoint inclusion. The interval must be clamped to non-negative first: a
  # reported SD of 0 unrounds to something like [-0.005, 0.005), and squaring
  # that negative endpoint would put a positive floor under ss and wrongly
  # exclude the zero-variance sample. Clamping also makes the lower endpoint
  # attainable, so its exclusion flag no longer applies.
  sd_lo <- max(0, iv_s$lo)
  sd_hi <- max(0, iv_s$hi)
  if (sd_hi <= 0 && !iv_s$hi_incl) {
    return(NULL)
  }
  lo_incl <- iv_s$lo_incl || iv_s$lo < 0
  ss_lo <- (n - 1) * (mg * sd_lo)^2
  ss_hi <- (n - 1) * (mg * sd_hi)^2
  base <- S^2 / n
  data.frame(
    S = S,
    Q_lo = base + ss_lo,
    Q_hi = base + ss_hi,
    lo_incl = lo_incl,
    hi_incl = iv_s$hi_incl
  )
}
