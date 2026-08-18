# Targeted reachability for brimmest(): certify ONE reported tuple without
# enumerating the whole attainable lattice.
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
# So the question becomes: can n integers in [0, W] hit one of a handful of
# (S, Q-window) targets? Three consequences:
#
#  (1) The reported values are inverse-rounded to intervals once, instead of
#      the lattice being forward-rounded cell by cell for every rounding rule.
#  (2) States that cannot reach any target are never visited. A partial sum s
#      after t items is viable only if the remaining k = n - t items can still
#      carry it to some target sum, and likewise on the R axis, which confines
#      the DP to a corridor around the targets rather than the full rectangle.
#  (3) A hit can be declared the moment it appears, at any layer: an item with
#      value y = 0 adds nothing to either coordinate, so a target reached
#      after t < n items pads out to exactly n with scores at the scale
#      minimum. Possibility is therefore often near-instant, while
#      impossibility still costs a full sweep of the corridor -- asymmetric,
#      but in the useful direction.
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

# Internal: can n integers in [0, W] realise any of the given (S, Q-window)
# targets? Exact in both directions. Returns TRUE / FALSE, or NA when the
# corridor would still exceed max_cells (caller falls back to the lattice).
.attainable_target <- function(W, n, tg, max_cells = 2e7) {
  if (is.null(tg) || !nrow(tg)) {
    return(FALSE)
  }
  ys <- 0:W
  pr <- ys * (W - ys) # each item's contribution to R
  g <- Reduce(.gcd2, pr[pr > 0])
  if (!length(g) || is.na(g) || g < 1) {
    g <- 1
  }
  pm <- max(pr) / g
  drs <- pr / g

  # R = W*S - Q, on the same reduced axis the full lattice uses. Every term
  # y(W - y) is divisible by g, so R always is.
  r_lo <- ceiling((W * tg$S - tg$Q_hi) / g - 1e-9)
  r_hi <- floor((W * tg$S - tg$Q_lo) / g + 1e-9)
  # a strictly excluded SD endpoint cannot supply a target state exactly on it
  ex_hi <- !tg$hi_incl & abs((W * tg$S - tg$Q_hi) / g - r_lo) < 1e-9
  ex_lo <- !tg$lo_incl & abs((W * tg$S - tg$Q_lo) / g - r_hi) < 1e-9
  r_lo <- r_lo + ex_hi
  r_hi <- r_hi - ex_lo
  keep <- r_hi >= 0 & r_lo <= n * pm & r_lo <= r_hi
  if (!any(keep)) {
    return(FALSE)
  }
  tg <- tg[keep, , drop = FALSE]
  r_lo <- pmax(0, r_lo[keep])
  r_hi <- pmin(n * pm, r_hi[keep])

  S_min <- min(tg$S)
  S_max <- max(tg$S)
  R_min <- min(r_lo)
  R_max <- max(r_hi)

  # Per-layer corridor: after t items, a partial (s, r) is viable only if the
  # remaining k = n - t items can still bridge to some target.
  s_win <- function(t) c(max(0L, S_min - (n - t) * W), min(t * W, S_max))
  r_win <- function(t) c(max(0, R_min - (n - t) * pm), min(t * pm, R_max))

  cells <- max(vapply(
    seq_len(n),
    function(t) {
      sw <- s_win(t)
      rw <- r_win(t)
      if (sw[2] < sw[1] || rw[2] < rw[1]) {
        return(0)
      }
      (sw[2] - sw[1] + 1) * (rw[2] - rw[1] + 1)
    },
    numeric(1)
  ))
  if (cells > max_cells) {
    return(NA)
  }

  hit <- function(sa, ra, M) {
    for (i in seq_len(nrow(tg))) {
      si <- tg$S[i] - sa + 1L
      if (si < 1L || si > nrow(M)) {
        next
      }
      lo <- max(r_lo[i], ra) - ra + 1L
      hi <- min(r_hi[i], ra + ncol(M) - 1L) - ra + 1L
      if (hi < lo) {
        next
      }
      if (any(M[si, lo:hi] != as.raw(0))) return(TRUE)
    }
    FALSE
  }

  z <- as.raw(0)
  sa <- 0L
  ra <- 0 # window origin of the live layer
  A <- matrix(z, 1L, 1L)
  A[1L, 1L] <- as.raw(1) # empty sample: (0, 0)

  for (t in seq_len(n)) {
    sw <- s_win(t)
    rw <- r_win(t)
    if (sw[2] < sw[1] || rw[2] < rw[1]) {
      return(FALSE)
    }
    nsa <- sw[1]
    nra <- rw[1]
    B <- matrix(z, sw[2] - sw[1] + 1L, rw[2] - rw[1] + 1L)
    for (i in seq_along(ys)) {
      v <- ys[i]
      dr <- drs[i]
      # source rows/cols in A that land inside B after the (v, dr) shift
      s0 <- max(sa, nsa - v)
      s1 <- min(sa + nrow(A) - 1L, sw[2] - v)
      if (s1 < s0) {
        next
      }
      r0 <- max(ra, nra - dr)
      r1 <- min(ra + ncol(A) - 1L, rw[2] - dr)
      if (r1 < r0) {
        next
      }
      bi <- (s0 + v - nsa + 1L):(s1 + v - nsa + 1L)
      bj <- (r0 + dr - nra + 1L):(r1 + dr - nra + 1L)
      B[bi, bj] <- B[bi, bj] |
        A[
          (s0 - sa + 1L):(s1 - sa + 1L),
          (r0 - ra + 1L):(r1 - ra + 1L),
          drop = FALSE
        ]
    }
    A <- B
    sa <- nsa
    ra <- nra
    # a target met now pads to exactly n with scores at the scale minimum
    if (hit(sa, ra, A)) return(TRUE)
  }
  FALSE
}
