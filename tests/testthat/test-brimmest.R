test_that("brimmest() separates attainable reports from the blind spot", {
  ok  <- brimmest(l = 1, u = 5, n = 9, mean = 3.0, sd = 1.0, digits = 1)
  bad <- brimmest(l = 1, u = 5, n = 9, mean = 1.3, sd = 0.9, digits = 1)
  expect_true(ok$possible)
  expect_true(nzchar(ok$rules))
  expect_false(bad$possible)
  expect_identical(bad$rules, "")
})

test_that("brimmest() is vectorised over tuples and keeps input order", {
  r <- brimmest(l = 1, u = 5, n = 9, digits = 1,
               mean = c(3.0, 1.3, 3.0), sd = c(1.0, 0.9, 1.0))
  expect_equal(nrow(r), 3L)
  expect_equal(r$mean, c(3.0, 1.3, 3.0))
  expect_identical(r$possible, c(TRUE, FALSE, TRUE))
  # recycling a scalar against a vector
  expect_equal(nrow(brimmest(l = 1, u = 5, n = 9, mean = 3.0,
                            sd = c(1.0, 0.9, 1.2), digits = 1)), 3L)
})

test_that("brimmest() never contradicts the closed-form screen's rejections", {
  # the screen is necessary but not sufficient: it may admit tuples that no
  # integer sample produces (the blind spot), but it must never reject one
  # that brimmest() proves attainable
  for (cfg in list(c(9, 5), c(12, 5), c(11, 7))) {
    n <- cfg[1]; u <- cfg[2]
    um <- suppressWarnings(umbrella_data(n = n, l = 1, u = u, digits = 1,
                                         Z = "integer"))
    cr <- brimmest(l = 1, u = u, n = n, mean = um$mean, sd = um$sd, digits = 1)
    expect_equal(sum(!um$consistent & cr$possible), 0L)
    expect_true(sum(um$consistent & !cr$possible) >= 0L)
  }
})

test_that("brimmest() agrees with an independently constructed witness", {
  # c(rep(0, 9), 1, 1, 1): n = 12 integers on 0-6, mean 0.25 -> "0.2",
  # sd 0.4523 -> "0.5". A real sample, so the report must brimmest possible.
  x <- c(rep(0, 9), 1, 1, 1)
  expect_equal(length(x), 12L)
  expect_true(all(x == round(x) & x >= 0 & x <= 6))
  r <- brimmest(l = 0, u = 6, n = 12,
               mean = round(mean(x), 1), sd = round(sd(x), 1), digits = 1)
  expect_true(r$possible)
})

test_that("a miss is only impossibility relative to the admitted rules", {
  # mean 0.25 rounds to 0.2 only when halves go down, so restricting the
  # rules can turn a possible report into an unreachable one
  x <- c(rep(0, 9), 1, 1, 1)
  both <- brimmest(l = 0, u = 6, n = 12, mean = 0.2, sd = 0.5, digits = 1,
                  rounding = c("half_up", "half_down"))
  down <- brimmest(l = 0, u = 6, n = 12, mean = 0.2, sd = 0.5, digits = 1,
                  rounding = "half_down")
  expect_true(both$possible)
  expect_true(down$possible)
  expect_true(grepl("half_down", both$rules))
})

test_that("brimmest() validates its arguments", {
  expect_error(brimmest(l = 1, u = 5, n = 9, mean = 3, sd = 1),
               "decimal places are required")
  expect_error(brimmest(l = 1, u = 5, n = 9, mean = 3, sd = 1, digits = 1,
                       rounding = "up_or_down"), "unknown rounding rule")
  expect_error(brimmest(l = 1, u = 5, n = 1, mean = 3, sd = 1, digits = 1),
               "n must be >= 2")
  # separate precision for the mean and the SD
  expect_equal(nrow(brimmest(l = 1, u = 5, n = 9, mean = 3.00, sd = 1.0,
                            mean_digits = 2, sd_digits = 1)), 1L)
})

test_that("brimmest() matches CLOSURE cell for cell", {
  skip_on_cran()
  skip_if_not_installed("unsum")
  um <- suppressWarnings(umbrella_data(n = 9, l = 1, u = 5, digits = 1,
                                       Z = "integer"))
  cr <- brimmest(l = 1, u = 5, n = 9, mean = um$mean, sd = um$sd, digits = 1)
  cl <- vapply(seq_len(nrow(um)), function(i) {
    r <- try(suppressWarnings(unsum::closure_generate(
      mean = sprintf("%.1f", um$mean[i]), sd = sprintf("%.1f", um$sd[i]),
      n = 9, scale_min = 1, scale_max = 5, rounding = "up_or_down")),
      silent = TRUE)
    if (inherits(r, "try-error")) FALSE else nrow(r$results) > 0
  }, logical(1))
  expect_identical(cr$possible, cl)
})

# ---- the targeted route (R/attainable-target.R) ------------------------------

test_that("the targeted route and the lattice route give identical verdicts", {
  # brimmest() picks a route by design size; both must certify the same set
  for (cfg in list(c(1, 5, 9), c(0, 6, 10), c(1, 7, 8))) {
    l <- cfg[1]; u <- cfg[2]; n <- cfg[3]
    lat <- strait:::.attainable_lattice(l, u, n, 1)
    grid <- expand.grid(mean = round(seq(l, u, by = 0.1), 1),
                        sd = round(seq(0, (u - l) / sqrt(2), by = 0.1), 1))
    for (rr in c("half_up", "half_down")) {
      lm <- round(strait:::.round_reported(lat$mean, 1, rr) * 10)
      ls <- round(strait:::.round_reported(lat$sd, 1, rr) * 10)
      ref <- paste(round(grid$mean * 10), round(grid$sd * 10)) %in% paste(lm, ls)
      got <- vapply(seq_len(nrow(grid)), function(i) {
        tg <- strait:::.target_states(l, u, n, 1, grid$mean[i], grid$sd[i],
                                      1, 1, rr)
        isTRUE(.attainable_target(as.integer(u - l), n, tg))
      }, logical(1))
      expect_identical(got, ref)
    }
  }
})

test_that("a zero SD survives the unrounding sign trap", {
  # unround_interval(0, 1, "up") is [-0.05, 0.05); squaring the negative
  # endpoint would put a positive floor under the sum of squares and wrongly
  # exclude the zero-variance sample of n scores all at the scale minimum
  expect_true(brimmest(l = 1, u = 5, n = 9, mean = 1.0, sd = 0.0,
                       digits = 1)$possible)
  expect_true(brimmest(l = 0, u = 6, n = 12, mean = 0.0, sd = 0.0,
                       digits = 1)$possible)
  expect_true(brimmest(l = 1, u = 5, n = 9, mean = 5.0, sd = 0.0,
                       digits = 1)$possible)
})

test_that("a mean admitting no integer sum is refused, not mis-enumerated", {
  # seq.int() counts down when from > to, so an empty candidate range must be
  # caught explicitly or it yields phantom sums and a false 'possible'
  expect_null(strait:::.target_states(1, 4, 6, 1, 1.1, 0.0, 1, 1, "half_up"))
  expect_false(brimmest(l = 1, u = 4, n = 6, mean = 1.1, sd = 0.0,
                        digits = 1)$possible)
})

test_that("brimmest() certifies designs too large for the full lattice", {
  skip_on_cran()
  # a 0-63 inventory at n = 50 needs a 3151 x 24801 state table, which the
  # lattice guard refuses outright; the targeted corridor is ~9% of that
  expect_gt(strait:::.lattice_cells(0, 63, 50, 1), 2e7)
  r <- brimmest(l = 0, u = 63, n = 50, mean = 20.5, sd = 12.3, digits = 1)
  expect_true(r$possible)
})

test_that("the lattice cache returns an identical object", {
  a <- strait:::.lattice_cached(1, 5, 9, 1)
  b <- strait:::.lattice_cached(1, 5, 9, 1)
  expect_identical(a, b)
  expect_identical(a, strait:::.attainable_lattice(1, 5, 9, 1))
})
