test_that("sd_bounds_curve() is hole-free under quasi-integer", {
  cv <- sd_bounds_curve(l = 1, u = 7, n = 30, granularity = "quasiinteger")
  expect_true(all(c("mean", "min_sd", "max_sd", "feasible", "pomp_mean",
                    "parity_max", "ceil_parity", "floor_parity") %in% names(cv)))
  expect_false(any(is.na(cv$min_sd[cv$feasible])))
  expect_false(any(is.na(cv$max_sd[cv$feasible])))
})

test_that("umbrella_data() passing set matches an independent GRIMMER + bounds screen", {
  # kept small deliberately: CRAN scrutiny (0.6.1) evaluates GRIMMER about 25x
  # slower than 0.6.2, so a wide grid here costs minutes on the CRAN path
  um <- umbrella_data(n = 12, l = 1, u = 3, digits = 1)
  expect_true(all(c("mean", "sd", "min_sd", "max_sd", "in_bounds",
                    "grimmer", "consistent") %in% names(um)))
  # all three verdict categories must be present, or this proves little
  expect_true(any(um$consistent))
  expect_true(any(um$in_bounds & !is.na(um$grimmer) & !um$grimmer))
  expect_true(any(!um$in_bounds))
  indep <- with(um, {
    gm <- as.logical(scrutiny::grimmer(x = mean, sd = sd, n = 12,
                                       digits_x = 1, digits_sd = 1,
                                       items = 1, rounding = "up_or_down"))
    ib <- (sd + 0.05) >= min_sd - 1e-9 & (sd - 0.05) <= max_sd + 1e-9
    gm & ib
  })
  expect_equal(um$consistent, indep)
})

test_that("the plot functions return ggplot objects", {
  cv  <- sd_bounds_curve(l = 1, u = 7, n = 30, granularity = "quasiinteger")
  pts <- brimmer_map(
    tibble::tibble(mean = c(3, 4, 5), sd = c(2.5, 0.2, 1.8)),
    l = 1, u = 7, n = 30, digits_mean = 1, digits_sd = 1,
    granularity = "quasiinteger")
  plot_sd_bounds(cv, points = pts) |> expect_s3_class("ggplot")
  plot_sd_bounds_pomp(cv, points = pts, reference = "parity") |> expect_s3_class("ggplot")
  plot_sd_bounds_pomp(cv, points = pts, reference = "sharp") |> expect_s3_class("ggplot")
  um <- umbrella_data(n = 12, l = 1, u = 7, digits = 2)
  expect_s3_class(plot_umbrella(um, curve = sd_bounds_curve(l = 1, u = 7, n = 12,
                                                            granularity = "quasiinteger")), "ggplot")
})

# Regression: plot_sd_bounds() used to infer the mean-grid spacing from
# median(diff(curve$mean)). sd_bounds_curve()'s grid is deliberately not
# uniform - it adds every kink of the 1/(n * n_items) lattice plus a pair of
# neighbours 1e-9 away - so once the kinks outnumber the uniform grid (about
# n * (u - l) > 333 at the default `by`) the median fell BELOW the plain grid
# spacing. band_polygon() then read every ordinary interval as a gap and the
# feasible region was drawn as thousands of one-column slivers instead of one
# band. sd_bounds_curve() now records the spacing it used and plot_sd_bounds()
# reads it.

# the white knock-out rings plot_sd_bounds() builds under shade = "outside"
n_rings <- function(p) {
  for (ly in p$layers) {
    if (
      inherits(ly$geom, "GeomPolygon") &&
        is.data.frame(ly$data) &&
        "ring" %in% names(ly$data)
    ) {
      return(length(unique(ly$data$ring)))
    }
  }
  NA_integer_
}

test_that("sd_bounds_curve() records the grid spacing it used", {
  attr(sd_bounds_curve(l = 1, u = 7, n = 30), "step") |> expect_equal(6 / 1000)
  attr(sd_bounds_curve(l = 1, u = 7, n = 30, by = 0.01), "step") |> expect_equal(0.01)
})

test_that("a hole-free band is one ring, however dense the kink lattice", {
  # n = 30 was already fine; 60 and 200 are the cases the median guess broke,
  # at 921 and 2001 rings respectively
  for (n in c(30, 60, 200)) {
    cv <- sd_bounds_curve(l = 1, u = 7, n = n, granularity = "quasiinteger")
    expect_equal(n_rings(plot_sd_bounds(cv)), 1L,
                 info = sprintf("l = 1, u = 7, n = %d", n))
  }
  # and on a wider scale, where the threshold is crossed sooner
  cv <- sd_bounds_curve(l = 0, u = 10, n = 80, granularity = "quasiinteger")
  n_rings(plot_sd_bounds(cv)) |> expect_equal(1L)
})

test_that("genuine gaps in the band still split it into separate rings", {
  # strictly integer data is feasible only at GRIM-consistent means, so the
  # band is a comb: merging those into one ring would shade the impossible
  # means between the teeth as feasible, the error band_polygon() exists to
  # prevent. The fix must not over-merge.
  cv <- sd_bounds_curve(l = 1, u = 7, n = 30, granularity = "integer")
  expect_gt(n_rings(plot_sd_bounds(cv)), 50L)
})

test_that("a hand-built uniform curve still works without the attribute", {
  cv <- sd_bounds_curve(l = 1, u = 7, n = 60, granularity = "quasiinteger")
  flat <- tibble::tibble(mean = seq(1, 7, by = 0.006))
  bd <- do.call(rbind, lapply(flat$mean, function(m)
    sd_bounds(l = 1, u = 7, n = 60, mean = m, granularity = "quasiinteger")))
  flat$min_sd <- bd$min_sd
  flat$max_sd <- bd$max_sd
  flat$feasible <- bd$feasible
  expect_null(attr(flat, "step"))
  n_rings(plot_sd_bounds(flat)) |> expect_equal(1L)
})
