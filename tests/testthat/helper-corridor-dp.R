# The corridor DP that brimmest() used to fall back on, kept here as the
# independent oracle for .certify_fast().
#
# It answers exactly the question .certify_fast() answers -- can n integers in
# [0, W] realise any of these (S, Q-window) targets? -- by a different method:
# a dynamic program over reachable (sum, sum-of-squares) states, confined to a
# corridor around the targets. Two implementations agreeing cell for cell over
# a whole reporting grid is the real check on the fast path, so the DP earns
# its keep in the suite even though the package no longer ships it.
#
# It left R/ because it was dead weight there: the fast path settled every one
# of 2,719,364 grid cells across six designs without ever exhausting its node
# budget, so the fallback branch never ran in production. See
# paper/method-brimmest/certification-fast-path-plan.md.
#
# Depends on strait-internal .gcd2(); everything else is self-contained.

# Internal: can n integers in [0, W] realise any of the given (S, Q-window)
# targets? Exact in both directions. Returns TRUE / FALSE, or NA when the
# corridor would still exceed max_cells (caller falls back to the lattice).
.attainable_target <- function(W, n, tg, max_cells = 2e7) {
  if (is.null(tg) || !nrow(tg)) {
    return(FALSE)
  }
  ys <- 0:W
  pr <- ys * (W - ys) # each item's contribution to R
  g <- Reduce(strait:::.gcd2, pr[pr > 0])
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
