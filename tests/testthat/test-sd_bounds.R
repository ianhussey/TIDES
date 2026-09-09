test_that("sd_bounds() returns the documented one-row data frame", {
  b <- sd_bounds(l = 1, u = 5, n = 9, mean = 2)
  expect_s3_class(b, "data.frame")
  nrow(b) |> expect_equal(1L)
  expect_true(all(c("min_sd", "max_sd", "feasible", "min_rule", "max_rule",
                    "grim", "grimmer", "sd_in_bounds", "note") %in% names(b)))
})

test_that("with no constraints the only bound is s >= 0", {
  b <- sd_bounds()
  expect_equal(b$min_sd, 0)
  expect_true(is.infinite(b$max_sd))
})

test_that("nested constraints only tighten the ceiling", {
  span   <- sd_bounds(l = 1, u = 5)$max_sd
  parity <- sd_bounds(l = 1, u = 5, n = 9)$max_sd
  sharp  <- sd_bounds(l = 1, u = 5, n = 9, mean = 2)$max_sd
  expect_lte(parity, span + 1e-9)
  expect_lte(sharp, parity + 1e-9)
})

test_that("the floor never exceeds the ceiling on feasible cells", {
  b <- sd_bounds(l = 1, u = 7, n = 25, mean = 4.12, granularity = "quasiinteger")
  expect_lte(b$min_sd, b$max_sd + 1e-9)
})

test_that("bounds match exhaustive integer enumeration (walls)", {
  bad <- 0L
  for (n in 3:5) {
    g <- as.matrix(expand.grid(rep(list(1:5), n)))
    S <- rowSums(g); sdv <- apply(g, 1, sd)
    for (s in unique(S)) {
      i <- S == s
      r <- sd_bounds(l = 1, u = 5, n = n, mean = s / n, granularity = "integer")
      if (!(r$feasible &&
            abs(r$min_sd - min(sdv[i])) < 1e-9 &&
            abs(r$max_sd - max(sdv[i])) < 1e-9)) bad <- bad + 1L
    }
  }
  expect_equal(bad, 0L)
})

test_that("attained extremes create a nonzero floor and a feasible band", {
  r <- sd_bounds(a = 2, b = 4, n = 9, mean = 3, granularity = "integer")
  expect_gt(r$min_sd, 0)
  # a mean outside the [a + W/n, b - W/n] band is infeasible
  bad <- sd_bounds(a = 2, b = 4, n = 9, mean = 2.05, granularity = "integer")
  expect_false(bad$feasible)
})

test_that("meanscored granularity equals the affine reduction, / n_items", {
  m <- 3
  r  <- sd_bounds(l = 1, u = 5, n = 9, mean = 3, granularity = "integer",
                  scoring = "meanscored", n_items = m)
  r2 <- sd_bounds(l = 1 * m, u = 5 * m, n = 9, mean = 3 * m, granularity = "integer")
  expect_equal(r$min_sd, r2$min_sd / m)
  expect_equal(r$max_sd, r2$max_sd / m)
})

test_that("meanscored + alpha reduces to sumscored + alpha, / n_items", {
  m <- 4
  ma <- sd_bounds(l = 1, u = 5, n = 25, mean = 3.2, granularity = "integer",
                  scoring = "meanscored", n_items = m, alpha = 0.75)
  sa <- sd_bounds(l = 1 * m, u = 5 * m, n = 25, mean = 3.2 * m, granularity = "integer",
                  scoring = "sumscored", n_items = m, alpha = 0.75)
  expect_equal(ma$min_sd, sa$min_sd / m)
  expect_equal(ma$max_sd, sa$max_sd / m)
})

test_that("the granularity floor is robust to scaling dust (no spurious floor)", {
  # A rounded mean scaled onto the integer grid can land a hair off an integer
  # (round(4/3, 10) * 3 = 3.9999999999). The floor must still be exactly 0
  # (an all-equal integer sample exists), not a spurious ~1e-6.
  r <- sd_bounds(l = 1, u = 5, n = 9, mean = round(4 / 3, 10), granularity = "integer",
                 scoring = "meanscored", n_items = 3)
  expect_true(r$feasible)
  expect_equal(r$min_sd, 0)
})

test_that("input guards fire", {
  expect_error(sd_bounds(scoring = "singleitem", n_items = 2))
  expect_error(sd_bounds(l = 1, u = 5, n = 9, mean = 3, alpha = 0.7,
                         scoring = "singleitem"))
  expect_error(sd_bounds(l = 5, u = 1))
  expect_error(sd_bounds(n = 1, mean = 3))
})
