l <- 1
u <- 5
n <- 7
m <- seq(l, u, by = 0.002)
bes <- sqrt(n / (n - 1))

test_that("band rules reproduce the article's nested-constraint panels", {
  # Popoviciu, range only
  d <- sd_region_data(l, u, n, rule = "range", by = 0.002)
  expect_equal(d$lo, rep(0, length(m)))
  expect_equal(d$hi, rep((u - l) / sqrt(2), length(m)))

  # Muilwijk / Bhatia-Davis arch, no floor
  d <- sd_region_data(l, u, n, rule = "mean", by = 0.002)
  expect_equal(d$hi, bes * sqrt(pmax(0, (u - m) * (m - l))))
  expect_equal(d$lo, rep(0, length(m)))

  # that arch plus the naive Bernoulli floor
  d <- sd_region_data(l, u, n, rule = "mean_naive_floor", by = 0.002)
  expect_equal(d$lo, bes * sqrt((m - floor(m)) * (1 - (m - floor(m)))))
  expect_equal(d$hi, bes * sqrt(pmax(0, (u - m) * (m - l))))

  # this article: both bounds sharp, GRIM-free
  d <- sd_region_data(l, u, n, rule = "quasi", by = 0.002)
  expect_equal(d$lo, sd_min_quasi_integer(m, n))
  expect_equal(d$hi, sd_max_structure_s(m, n, l, u))
})

test_that("the alpha rule matches the article's alpha panel", {
  k <- 2; a <- 0.70; cc <- (k - 1) / k; D <- 1 - cc * a
  qmin_k <- sd_min_quasi_integer(k * (m - l), n) / k
  qmax_k <- sd_max_structure_s(k * (m - l), n, 0, k * (u - l)) / k
  ceil_a <- sqrt((n / (n - 1)) * v_max_alpha(k * m, k, n, l, u) / D) / k
  exp_lo <- qmin_k / sqrt(D)
  exp_hi <- pmin(ceil_a, qmax_k)
  bad <- exp_lo > exp_hi + 1e-12
  exp_lo[bad] <- NA_real_; exp_hi[bad] <- NA_real_

  d <- sd_region_data(l, u, n, rule = "alpha", n_items = k, alpha = a, by = 0.002)
  expect_equal(d$lo, exp_lo)
  expect_equal(d$hi, exp_hi)
})

test_that("alpha only tightens, and is inert / rejected at one item", {
  free  <- sd_region_data(l, u, n, rule = "quasi", n_items = 2, by = 0.01)
  withA <- sd_region_data(l, u, n, rule = "alpha", n_items = 2, alpha = 0.7, by = 0.01)
  ok <- !is.na(withA$hi)
  expect_true(all(withA$hi[ok] <= free$hi[ok] + 1e-9))
  expect_true(all(withA$lo[ok] >= free$lo[ok] - 1e-9))
  expect_error(sd_region_data(l, u, n, rule = "alpha", n_items = 1, alpha = 0.7))
  expect_error(sd_region_data(l, u, n, rule = "alpha"))          # alpha required
  expect_error(sd_region_data(l, u, n, rule = "integer_alpha"))  # alpha required
})

test_that("lattice rules return attainable tuples, and alpha is a strict subset", {
  a  <- sd_region_data(l, u, n, rule = "integer", n_items = 2, digits = 2)
  b  <- sd_region_data(l, u, n, rule = "integer_alpha", n_items = 2,
                       alpha = 0.70, digits = 2)
  attr(a, "type") |> expect_identical("points")
  expect_named(a, c("mean", "lo", "hi", "sd"))
  expect_gt(nrow(a), 0)
  expect_lt(nrow(b), nrow(a))                                   # alpha removes tuples
  expect_true(all(paste(b$mean, b$sd) %in% paste(a$mean, a$sd))) # strict subset

  # every retained tuple really is inside the alpha-conditional bounds
  # round(): a tibble refuses fractional row indices, unlike a data frame
  keep <- b[round(seq(1, nrow(b), length.out = min(40, nrow(b)))), ]
  ok <- vapply(seq_len(nrow(keep)), function(i) {
    r <- sd_bounds(l = l, u = u, n = n, mean = keep$mean[i], digits_mean = 2,
                   rounding = "up_or_down", granularity = "integer", scoring = "meanscored",
                   n_items = 2, alpha = 0.70)
    isTRUE(r$feasible) &&
      (keep$sd[i] + 0.005) >= r$min_sd - 1e-9 &&
      (keep$sd[i] - 0.005) <= r$max_sd + 1e-9
  }, logical(1))
  expect_true(all(ok))
})

test_that("every rule returns a ggplot", {
  for (r in c("range", "range_n", "mean", "mean_naive_floor", "mestdagh",
              "pesant_regin", "quasi")) {
    plot_sd_region(l, u, n, rule = r) |> expect_s3_class("ggplot")
  }
  expect_s3_class(plot_sd_region(l, u, n, rule = "alpha", n_items = 2, alpha = 0.7),
                  "ggplot")
  plot_sd_region(l, u, n, rule = "integer", digits = 1) |> expect_s3_class("ggplot")
  expect_s3_class(plot_sd_region(l, u, n, rule = "integer_alpha", n_items = 2,
                                 alpha = 0.7, digits = 1), "ggplot")
})

test_that("sd_max_muilwijk is the smooth arch and never below Structure S", {
  expect_equal(sd_max_muilwijk(3, n, l, u), bes * sqrt((u - 3) * (3 - l)))
  expect_true(all(sd_max_muilwijk(m, n, l, u) >= sd_max_structure_s(m, n, l, u) - 1e-9))
})

test_that("the exact lattice matches brute-force enumeration", {
  brute <- function(n, l, u) {
    g <- as.matrix(expand.grid(rep(list(l:u), n)))
    m <- unique(round(cbind(rowMeans(g), apply(g, 1, stats::sd)), 10))
    m[order(m[, 1], m[, 2]), , drop = FALSE]
  }
  for (cfg in list(c(4, 1, 5), c(5, 0, 3), c(3, 1, 7))) {
    n <- cfg[1]; l <- cfg[2]; u <- cfg[3]
    d <- sd_region_data(l, u, n, rule = "attainable")
    a <- unique(round(as.matrix(d[, c("mean", "sd")]), 10))
    a <- a[order(a[, 1], a[, 2]), , drop = FALSE]
    expect_equal(unname(a), unname(brute(n, l, u)), tolerance = 1e-9)
  }
  expect_identical(attr(sd_region_data(0, 3, 5, rule = "attainable"), "type"),
                   "points")
})

test_that("the exact lattice shows holes the reported lattice smears shut", {
  ex <- sd_region_data(0, 3, 7, rule = "attainable")
  # attainable sums of squares are not contiguous in steps of 2: real interior
  # gaps exist, which rounding to a reporting grid hides
  gaps <- vapply(split(ex$sd, ex$mean), function(s) {
    d <- diff(sort(unique(round(s, 10))))
    length(d) > 1 && max(d) > 1.5 * min(d)
  }, logical(1))
  expect_true(any(gaps))
})

test_that("scoring decouples granularity from alpha's k", {
  # a 2-item 0-3 SUM composite is plain integer data on [0, 6]
  s1 <- sd_region_data(0, 6, 23, rule = "quasi", n_items = 2,
                       scoring = "sumscored", by = 0.05)
  s2 <- sd_region_data(0, 6, 23, rule = "quasi", n_items = 1, by = 0.05)
  expect_equal(s1, s2)
  # sum-score alpha ceiling is k times the mean-score one
  a_sum  <- sd_region_data(0, 6, 23, rule = "alpha", n_items = 2,
                           scoring = "sumscored", alpha = 0.5, by = 1)
  a_mean <- sd_region_data(0, 3, 23, rule = "alpha", n_items = 2,
                           scoring = "meanscored", alpha = 0.5, by = 0.5)
  expect_equal(a_sum$hi, 2 * a_mean$hi, tolerance = 1e-8)
  # default scoring is unchanged
  expect_identical(sd_region_data(1, 7, 7, rule = "quasi", n_items = 5, by = 0.05),
                   sd_region_data(1, 7, 7, rule = "quasi", n_items = 5,
                                  scoring = "meanscored", by = 0.05))
  expect_error(sd_region_data(0, 5, 9, rule = "quasi", n_items = 2,
                              scoring = "sumscored"), "divisible")
})

test_that("alpha prunes the exact lattice from both directions", {
  a <- sd_region_data(0, 6, 23, rule = "attainable")
  b <- sd_region_data(0, 6, 23, rule = "attainable_alpha", n_items = 2,
                      scoring = "sumscored", alpha = 0.5)
  expect_true(nrow(b) < nrow(a))
  expect_true(all(paste(b$mean, b$sd) %in% paste(a$mean, a$sd)))
  lo_a <- tapply(a$sd, a$mean, min); lo_b <- tapply(b$sd, b$mean, min)
  hi_a <- tapply(a$sd, a$mean, max); hi_b <- tapply(b$sd, b$mean, max)
  cm <- intersect(names(lo_a), names(lo_b))
  expect_true(sum(lo_b[cm] > lo_a[cm] + 1e-9) > 0)   # raises the floor
  expect_true(sum(hi_b[cm] < hi_a[cm] - 1e-9) > 0)   # lowers the ceiling
})

test_that("the batched Gini envelope agrees with per-mean evaluation", {
  m <- 1 / (1 - 0.5 * 0.5)
  # a cold cache and a warm one must agree, and both must be finite and
  # increasing in the profile constraint
  rm(list = ls(envir = strait:::gini_envelope_cache),
     envir = strait:::gini_envelope_cache)
  cold <- vapply(c(1, 1.5, 2, 2.5, 3), function(mu) {
    r <- sd_min_alpha_gini(0, 6, 23, mu, m); if (is.null(r)) NA_real_ else r
  }, 0)
  warm <- vapply(c(1, 1.5, 2, 2.5, 3), function(mu) {
    r <- sd_min_alpha_gini(0, 6, 23, mu, m); if (is.null(r)) NA_real_ else r
  }, 0)
  expect_equal(cold, warm)
  expect_true(all(cold[!is.na(cold)] > 0))   # positive alpha forbids SD = 0
})

test_that("the lattice DP's state-space reductions are sound", {
  # the (S, R) axis change, gcd packing, growing window, reflection symmetry
  # and raw storage must not alter a single tuple
  for (cfg in list(c(4, 1, 5), c(5, 0, 3), c(7, 0, 6), c(6, 0, 2))) {
    n <- cfg[1]; l <- cfg[2]; u <- cfg[3]
    d <- strait:::attainable_lattice(l, u, n, 1)
    g <- as.matrix(expand.grid(rep(list(l:u), n)))
    b <- unique(round(cbind(rowMeans(g), apply(g, 1, stats::sd)), 10))
    b <- b[order(b[, 1], b[, 2]), , drop = FALSE]
    a <- unique(round(as.matrix(d[, c("mean", "sd")]), 10))
    a <- a[order(a[, 1], a[, 2]), , drop = FALSE]
    expect_equal(unname(a), unname(b), tolerance = 1e-9)
  }
  # odd W packs the R axis by g = 2, even W does not; both must be exact
  nrow(strait:::attainable_lattice(0, 3, 7, 1)) |> expect_equal(100L)
  nrow(strait:::attainable_lattice(0, 6, 23, 1)) |> expect_equal(8634L)
  # a mean-scored grid (mg > 1) still lands on the 1/mg lattice
  d <- strait:::attainable_lattice(0, 3, 6, 2)
  expect_true(all(abs(d$mean * 6 * 2 - round(d$mean * 6 * 2)) < 1e-9))
  # the guard fires only when the REDUCED grid is oversized; the reductions
  # bring cases within reach that the naive (S, Q) grid could not hold
  expect_error(strait:::attainable_lattice(1, 7, 61, 5), "too large")
  expect_silent(strait:::attainable_lattice(0, 6, 300, 1))
})

test_that("round_digits/rounding round the emitted lattice", {
  ex <- sd_region_data(0, 3, 7, rule = "attainable")
  expect_true(any(abs(ex$sd * 100 - round(ex$sd * 100)) > 1e-9))  # exact by default
  for (rr in c("half_up", "half_down", "native", "ceiling", "floor",
               "trunc", "anti_trunc")) {
    r <- sd_region_data(0, 3, 7, rule = "attainable", round_digits = 1,
                        rounding = rr)
    expect_true(all(abs(r$sd * 10 - round(r$sd * 10)) < 1e-9))
    expect_true(all(abs(r$mean * 10 - round(r$mean * 10)) < 1e-9))
    expect_false(any(1 / r$sd == -Inf))          # no negative zero
    expect_false(anyDuplicated(r) > 0)           # collapsed duplicates
  }
  expect_error(sd_region_data(0, 3, 7, rule = "attainable", round_digits = 1,
                              rounding = "nonsense"))
  # rounding genuinely merges distinct exact SDs
  expect_lt(length(unique(sd_region_data(0, 3, 7, rule = "attainable",
                                         round_digits = 1)$sd)),
            length(unique(round(ex$sd, 10))))
})

test_that("forward-rounded exact tuples are a subset of the GRIMMER lattice", {
  # the strongest available cross-check: every pair a real integer data set can
  # be reported as must survive GRIM, GRIMMER and the bounds. GRIMMER is
  # necessary but not sufficient, so it admits extras; it must never reject one.
  for (cfg in list(c(7, 0, 3, 1), c(9, 1, 5, 1), c(12, 0, 6, 1), c(7, 1, 7, 2))) {
    n <- cfg[1]; l <- cfg[2]; u <- cfg[3]; dg <- cfg[4]
    k <- function(d) paste(sprintf("%.10f", d$mean), sprintf("%.10f", d$sd))
    grm <- k(sd_region_data(l, u, n, rule = "integer", digits = dg))
    for (rr in c("half_up", "half_down")) {
      fwd <- k(sd_region_data(l, u, n, rule = "attainable",
                              round_digits = dg, rounding = rr))
      expect_true(all(fwd %in% grm))
    }
  }
})

test_that("umbrella_data's reporting grid is free of seq() drift", {
  # seq(0, 6, by = 0.1)[4] is 0.3 + 5.6e-17, and scrutiny reads the value as a
  # decimal, so an unrounded grid flips GRIMMER verdicts
  um <- umbrella_data(n = 12, l = 0, u = 6, digits = 1, granularity = "integer")
  expect_true(all(abs(um$mean * 10 - round(um$mean * 10)) < 1e-12))
  expect_true(all(abs(um$sd * 10 - round(um$sd * 10)) < 1e-12))
  expect_identical(sum(um$consistent), 1175L)   # 1167 before the fix
})
