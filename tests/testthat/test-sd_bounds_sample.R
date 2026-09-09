test_that("sd_bounds_sample() attains the bound and the mean", {
  for (cfg in list(c(3, 30, 1, 7), c(89/30, 30, 1, 7), c(4.12, 25, 1, 7), c(2.44, 9, 1, 5))) {
    mean <- cfg[1]; n <- cfg[2]; l <- cfg[3]; u <- cfg[4]
    smax <- sd_bounds_sample(l, u, n, mean, "max")
    smin <- sd_bounds_sample(l, u, n, mean, "min")
    expect_length(smax, n); expect_length(smin, n)
    expect_equal(mean(smax), mean); expect_equal(mean(smin), mean)
    expect_equal(sd(smax), sd_bounds(l = l, u = u, n = n, mean = mean)$max_sd, tolerance = 1e-6)
    expect_equal(sd(smin), sd_bounds(l = l, u = u, n = n, mean = mean, granularity = "quasiinteger")$min_sd,
                 tolerance = 1e-6)
  }
})

test_that("granularity = 'integer' yields an all-integer sample and needs a GRIM mean", {
  x <- sd_bounds_sample(1, 7, 30, 89/30, "max", granularity = "integer")
  expect_true(all(x == round(x)))
  expect_error(sd_bounds_sample(1, 7, 30, 2.97, "max", granularity = "integer"))  # not GRIM
})

test_that("edge cases: mean at a limit gives SD 0", {
  sd(sd_bounds_sample(1, 7, 10, 1, "max")) |> expect_equal(0)
  expect_equal(sd(sd_bounds_sample(1, 7, 10, 4, "min")), 0)  # whole-number mean
})
