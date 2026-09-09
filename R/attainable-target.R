# Targeted certification for brimmest(): reduce ONE reported tuple to the
# handful of integer states it could have come from, without enumerating the
# whole attainable lattice.
#
# attainable_lattice() answers "which (mean, sd) pairs are reachable?" and then
# brimmest() asks whether the report is among them. For a single report that is
# a membership question answered by a reachability-of-everything computation,
# and on a wide scale the full lattice is not merely slow but refused outright
# by its own size guard: a 0-63 inventory at n = 50 needs a 3151 x 24801 state
# table.
#
# The report pins both coordinates almost completely, which is what this file
# exploits: the mean's rounding interval admits only a few integer sums (the
# GRIM condition, read as a constraint rather than a test), and each candidate
# sum then pins the sum of squares to a narrow window, because
# Q = S^2/n + (n - 1) * sd^2 and sd is known to its reporting precision.
#
# target_states() builds those (S, Q-window) targets by inverse-rounding the
# reported values once, rather than forward-rounding the lattice cell by cell.
# R/certify-sandwich.R then decides whether n integers in [0, W] can hit any of
# them. The verdict is identical to lattice membership; only the work differs,
# and that equivalence is asserted cell for cell in
# tests/testthat/test-brimmest.R.

# Internal: how many cells would the full lattice need? Mirrors the arithmetic
# in attainable_lattice() without allocating anything, so brimmest() can choose
# a route before committing to one. NA when the design is not on the integer
# grid at all.
lattice_cells <- function(l, u, n, multiplier) {
  span <- multiplier * (u - l)
  if (!dplyr::near(span, round(span), tol = TOLERANCE)) {
    return(NA_real_)
  }

  span <- as.integer(round(span))
  if (span < 1L) {
    return(NA_real_)
  }

  values <- 0:span
  contributions <- values * (span - values)
  divisor <- Reduce(gcd2, contributions[contributions > 0])
  if (!length(divisor) || is.na(divisor) || divisor < 1) {
    divisor <- 1
  }

  (n * span + 1) * (n * max(contributions) / divisor + 1)
}

# Internal: the full lattice, memoised per design. Enumerating it is the
# dominant cost of the lattice route, and repeated calls on one design are the
# common case when screening a table of reports.
lattice_cache <- new.env(parent = emptyenv())

lattice_cached <- function(l, u, n, multiplier, max_cells = MAX_LATTICE_CELLS) {
  key <- paste(l, u, n, multiplier, sep = "\r")
  hit <- get0(key, envir = lattice_cache, inherits = FALSE)
  if (!is.null(hit)) {
    return(hit)
  }

  out <- attainable_lattice(l, u, n, multiplier, max_cells = max_cells)
  assign(key, out, envir = lattice_cache)

  out
}

# Internal: map the forward rounding vocabulary used by the lattice
# (round_reported) onto the inverse-rounding vocabulary of unround_interval().
rounding_inverse <- function(rounding) {
  switch(
    rounding,
    half_up = "up",
    half_down = "down",
    native = "even",
    ceiling = "ceiling",
    floor = "floor",
    trunc = "trunc",
    anti_trunc = "anti_trunc",
    cli::cli_abort("Unknown rounding rule: {.val {rounding}}.")
  )
}

# Internal: the (S, Q) targets implied by one reported (mean, sd) under one
# rounding rule. Returns NULL when the report admits no integer sum at all (the
# GRIM condition failing), which is already a certificate of impossibility.
target_states <- function(
  l,
  u,
  n,
  multiplier,
  mean,
  sd,
  digits_mean,
  digits_sd,
  rounding
) {
  inverse <- rounding_inverse(rounding)
  span <- as.integer(round(multiplier * (u - l)))
  interval_mean <- unround_interval(mean, digits_mean, inverse)
  interval_sd <- unround_interval(sd, digits_sd, inverse)

  # S = n * multiplier * (mean - l) must be a whole number of grid steps.
  # Endpoint inclusion matters: a half-up report owns its lower endpoint but
  # not its upper, and an S sitting exactly on an excluded endpoint is not
  # admissible.
  sum_lo <- n * multiplier * (interval_mean$lo - l)
  sum_hi <- n * multiplier * (interval_mean$hi - l)
  sum_from <- as.integer(ceiling(sum_lo - TOLERANCE))
  sum_to <- as.integer(floor(sum_hi + TOLERANCE))
  # seq.int() counts DOWN when from > to, so an empty range has to be caught
  # explicitly; otherwise a mean admitting no integer sum yields two phantoms
  if (sum_to < sum_from) {
    return(NULL)
  }

  sums <- seq.int(sum_from, sum_to)
  if (!interval_mean$lo_incl) {
    sums <- sums[sums > sum_lo + TOLERANCE]
  }
  if (!interval_mean$hi_incl) {
    sums <- sums[sums < sum_hi - TOLERANCE]
  }
  sums <- sums[sums >= 0L & sums <= n * span]
  if (!length(sums)) {
    return(NULL)
  }

  # Q = S^2/n + (n - 1) * (multiplier * sd)^2, and the SD interval carries its
  # own endpoint inclusion. The interval must be clamped to non-negative first:
  # a reported SD of 0 unrounds to something like [-0.005, 0.005), and squaring
  # that negative endpoint would put a positive floor under the sum of squares
  # and wrongly exclude the zero-variance sample. Clamping also makes the lower
  # endpoint attainable, so its exclusion flag no longer applies.
  sd_lo <- max(0, interval_sd$lo)
  sd_hi <- max(0, interval_sd$hi)
  if (sd_hi <= 0 && !interval_sd$hi_incl) {
    return(NULL)
  }

  base <- sums^2 / n

  tibble::tibble(
    S = sums,
    Q_lo = base + (n - 1) * (multiplier * sd_lo)^2,
    Q_hi = base + (n - 1) * (multiplier * sd_hi)^2,
    lo_incl = interval_sd$lo_incl || interval_sd$lo < 0,
    hi_incl = interval_sd$hi_incl
  )
}
