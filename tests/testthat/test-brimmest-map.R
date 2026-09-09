# brimmest_map() must reach exactly the verdicts brimmest() would reach
# row by row. Everything it adds -- design grouping, tuple de-duplication --
# is an efficiency, so agreement with the per-row result is the whole contract.

test_that("grouped certification equals row-by-row certification", {
  set.seed(2)
  d <- tibble::tibble(mean = round(runif(60, 1, 5), 1),
                  sd   = round(runif(60, 0, 2), 1),
                  n    = sample(c(9, 12, 15), 60, replace = TRUE))
  got <- brimmest_map(d, l = 1, u = 5, digits_mean = 1, digits_sd = 1, include_inputs = FALSE)
  ref <- lapply(seq_len(nrow(d)), function(i)
    brimmest(l = 1, u = 5, n = d$n[i], mean = d$mean[i], sd = d$sd[i],
             digits_mean = 1, digits_sd = 1))
  expect_identical(got$possible, vapply(ref, function(x) x$possible, logical(1)))
  expect_identical(got$rules, vapply(ref, function(x) x$rules, character(1)))
  nrow(got) |> expect_identical(nrow(d))
})

test_that("repeated tuples and row order survive de-duplication", {
  d <- tibble::tibble(mean = c(3.0, 1.3, 3.0, 2.5, 1.3),
                  sd   = c(1.0, 0.9, 1.0, 1.2, 0.9))
  r <- brimmest_map(d, l = 1, u = 5, n = 9, digits_mean = 1, digits_sd = 1,
                         include_inputs = FALSE)
  expect_identical(r$possible, c(TRUE, FALSE, TRUE, FALSE, FALSE))
  # identical inputs must give identical outputs wherever they appear
  expect_identical(r$possible[1], r$possible[3])
  expect_identical(r$rules[2], r$rules[5])
})

test_that("include_inputs controls whether the inputs come back", {
  d <- tibble::tibble(mean = c(3.0, 1.3), sd = c(1.0, 0.9))
  with_in <- brimmest_map(d, l = 1, u = 5, n = 9, digits_mean = 1, digits_sd = 1)
  without <- brimmest_map(d, l = 1, u = 5, n = 9, digits_mean = 1, digits_sd = 1,
                               include_inputs = FALSE)
  names(with_in) |> expect_identical(c("mean", "sd", "possible", "rules"))
  names(without) |> expect_identical(c("possible", "rules"))
  expect_identical(with_in$possible, without$possible)
})

test_that("rounding is a set for the call, not a value per row", {
  # mean 0.25 rounds to 0.2 only when halves go down, so the admitted set
  # decides the verdict -- and it must be passed through unmangled rather
  # than recycled across rows the way a per-row constant would be
  d <- tibble::tibble(mean = 0.2, sd = 0.5)
  up <- brimmest_map(d, l = 0, u = 6, n = 12, digits_mean = 1, digits_sd = 1,
                          rounding = "half_up", include_inputs = FALSE)
  both <- brimmest_map(d, l = 0, u = 6, n = 12, digits_mean = 1, digits_sd = 1,
                            rounding = c("half_up", "half_down"),
                            include_inputs = FALSE)
  expect_false(up$possible)
  expect_true(both$possible)
  expect_identical(both$rules, "half_down")
  # supplying it per row is refused rather than silently read
  expect_error(brimmest_map(cbind(d, rounding = "half_up"),
                                 l = 0, u = 6, n = 12, digits_mean = 1, digits_sd = 1),
               "whole call")
})

test_that("brimmest_map() validates its inputs", {
  d <- tibble::tibble(mean = c(3.0, 1.3), sd = c(1.0, 0.9))
  expect_error(brimmest_map(list(mean = 1, sd = 1)), "data frame")
  expect_error(brimmest_map(d, mean = 3, l = 1, u = 5, n = 9, digits_mean = 1, digits_sd = 1),
               "both a column and a constant")
  expect_error(brimmest_map(d, l = 1, u = 5, digits_mean = 1, digits_sd = 1), "`n` is required")
  expect_error(brimmest_map(d, l = 1, u = 5, n = 9), "digits_mean")
  expect_error(brimmest_map(d, l = 1, u = 5, n = 9, digits_mean = 1, digits_sd = 1,
                                 granularity = "integer"), "[Uu]nknown constant")
  expect_error(brimmest_map(d[0, ], l = 1, u = 5, n = 9, digits_mean = 1, digits_sd = 1),
               "no rows")
})

test_that("the closed-form screen never rejects what brimmest_map() allows", {
  # the same nesting brimmer() and brimmest() obey, exercised through the two
  # batch entry points against each other
  d <- tidyr::expand_grid(mean = seq(1, 5, by = 0.1), sd = seq(0, 2, by = 0.1))
  screen <- brimmer_map(d, l = 1, u = 5, n = 9, digits_mean = 1,
                             digits_sd = 1, granularity = "integer",
                             include_inputs = FALSE)
  exact <- brimmest_map(d, l = 1, u = 5, n = 9, digits_mean = 1, digits_sd = 1,
                             include_inputs = FALSE)
  sum(!screen$consistent & exact$possible) |> expect_equal(0L)
  # and it does admit cells the exact test rejects: the blind spot
  expect_true(sum(screen$consistent & !exact$possible) > 0L)
})
