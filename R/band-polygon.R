# Geometry helpers shared by the plotting layer: feasible regions as closed
# rings, and the count-parity correction factor as a quantity in its own right.

#' A feasible region as closed rings
#'
#' Converts a `(mean, lo, hi)` band into one or more closed polygons, so that a
#' region can be drawn as a ring rather than as the area between two curves.
#'
#' Drawing a band with `ggplot2::geom_ribbon()` assumes the band exists at every
#' mean. Some constraint sets break that assumption: with a reported Cronbach's
#' alpha there are stretches near each scale limit where no composite satisfies
#' the constraints at all, and `sd_region_data()` returns `NA` there. A ribbon
#' simply draws nothing across such a gap. If the plot shades the *feasible*
#' region the gap is merely invisible; if it shades the *infeasible* region — as
#' [plot_sd_bounds()] and [plot_sd_region()] do by default — the gap becomes an
#' unshaded bar, which asserts that any SD whatsoever is possible at those
#' means. That is the one thing such a figure must not say.
#'
#' Rings avoid this by construction: shade the whole panel, knock the rings out,
#' and stroke their outlines. Anywhere the band is undefined stays shaded
#' because nothing was knocked out of it, rather than because the code
#' remembered to special-case it. A ring also closes itself, so a region ending
#' mid-scale gets a vertical edge instead of two curves trailing off.
#'
#' Runs are split wherever consecutive means are more than `1.5 * by` apart,
#' which distinguishes a genuine gap from floating-point jitter in `seq()`. A
#' run of a single row becomes a zero-area ring, which draws nothing visible;
#' these are left in rather than dropped, since whether an isolated feasible
#' point should be shown is the caller's decision.
#'
#' @param d A data frame with `mean`, `lo` and `hi`, e.g. from
#'   [sd_region_data()] under a band rule.
#' @param by The mean-grid spacing of `d`, used to tell a gap from jitter.
#' @return A long-format data frame with `mean`, `y` and a `ring` id, ready for
#'   `geom_polygon(aes(group = ring))`, or `NULL` if no row is finite.
#' @seealso [sd_region_data()], [plot_sd_region()].
#' @examples
#' # a reported alpha leaves stretches near each limit where no composite exists
#' d <- sd_region_data(0, 3, 7, rule = "alpha", n_items = 2, alpha = 0.70,
#'                     by = 0.001)
#' sum(is.na(d$lo))
#'
#' # the feasible region is therefore three rings, not one band
#' rings <- band_polygon(d, by = 0.001)
#' length(unique(rings$ring))
#' @export
band_polygon <- function(d, by) {
  d <- d[is.finite(d$lo) & is.finite(d$hi), , drop = FALSE]
  if (!nrow(d)) {
    return(NULL)
  }
  d <- d[order(d$mean), , drop = FALSE]
  run <- c(0, cumsum(diff(d$mean) > by * 1.5))
  out <- do.call(
    rbind,
    lapply(unique(run), function(k) {
      g <- d[run == k, , drop = FALSE]
      data.frame(
        mean = c(g$mean, rev(g$mean), g$mean[1]),
        y = c(g$hi, rev(g$lo), g$hi[1]),
        ring = k
      )
    })
  )
  rownames(out) <- NULL
  out
}


#' The outline of an umbrella grid
#'
#' The envelope of a lattice of reportable `(mean, sd)` tuples, as closed rings:
#' at each mean the lowest and highest consistent SD, handed to
#' [band_polygon()]. This is what `plot_umbrella(style = "contour")` draws, and
#' it is returned as data so the outline can be drawn elsewhere, or measured.
#'
#' `by` is the spacing that tells a true gap in the umbrella from its ordinary
#' striping, and it defaults to the median gap between the means that *have* a
#' consistent tuple rather than to the reporting grid's step. Under `Z =
#' "integer"` only the means an integer sum can round to survive, so on the
#' reporting grid almost every neighbouring pair of means is a gap; splitting
#' there would make each stripe its own one-column ring, and a one-column ring
#' has no width to draw. The median of the realised spacing is the striping's
#' own period, so the stripes join into the envelope and only a true stretch of
#' infeasible means — wider than `1.5 * by` — still splits the region.
#'
#' The contour therefore claims less than the tuples do: it says which `(mean,
#' sd)` pairs are near the region the tests admit, not which ones are in it.

#' @param umbrella Output of [umbrella_data()], or an already-filtered lattice
#'   from `sd_region_data(rule = "integer")` — anything with `mean` and `sd`,
#'   filtered by `consistent` when that column is present.
#' @param by Mean spacing that separates a true gap from the striping;
#'   defaults to the median spacing between the means that have a consistent
#'   tuple. Give it explicitly to keep a narrow true gap from being bridged.
#' @return A long-format data frame with `mean`, `y` and a `ring` id, as
#'   [band_polygon()] returns, or `NULL` if no tuple is consistent.
#' @seealso [plot_umbrella()] to draw it, [band_polygon()] for the ring
#'   construction, [umbrella_data()] for the lattice.
#' @examples
#' grid <- umbrella_data(n = 12, l = 1, u = 3, digits = 1)
#' rings <- umbrella_contour(grid)
#'
#' # one region here, and it spans the scale
#' length(unique(rings$ring))
#' range(rings$mean)
#' @export
umbrella_contour <- function(umbrella, by = NULL) {
  pts <- if ("consistent" %in% names(umbrella)) {
    umbrella[!is.na(umbrella$consistent) & umbrella$consistent, , drop = FALSE]
  } else {
    umbrella
  }

  ms <- sort(unique(pts$mean))

  if (!length(ms)) {
    return(NULL)
  }

  i <- match(pts$mean, ms)
  lo <- rep(NA_real_, length(ms))
  hi <- lo

  agg_lo <- tapply(pts$sd, i, min)
  agg_hi <- tapply(pts$sd, i, max)

  lo[as.integer(names(agg_lo))] <- agg_lo
  hi[as.integer(names(agg_hi))] <- agg_hi

  if (is.null(by)) {
    by <- if (length(ms) > 1) stats::median(diff(ms)) else 1
  }

  band_polygon(
    tibble::tibble(mean = ms, lo = lo, hi = hi),
    by = by
  )
}


#' The count-parity correction to the mean-conditional ceiling
#'
#' The factor relating Muilwijk's mean-conditional ceiling to the sharp one.
#'
#' Muilwijk's bound allows the counts at each scale limit to be fractional. They
#' cannot be: with `R = u - l`, a sample of `n` observations attaining the
#' ceiling puts `n_l = n(u - m)/R` at the lower limit and `n_u = n(m - l)/R` at
#' the upper, and those must be whole numbers. Rounding them costs
#'
#' \deqn{\delta = 1 - n \, \mathrm{frac}(n_l)\, \mathrm{frac}(n_u) / (n_l n_u)}
#'
#' so that the sharp maximum is
#' \eqn{s_{\max}^2 = \frac{n}{n-1}(u - m)(m - l)\,\delta}. [sd_max_muilwijk()]
#' computes the uncorrected arch and [sd_max_structure_s()] the corrected
#' maximum; this returns the correction between them.
#'
#' `delta` is interpretable on its own: it is the fraction of Muilwijk's ceiling
#' that is actually reachable. It equals 1 exactly when both counts are whole,
#' and falls toward 0 as the mean approaches a scale limit, where the ceiling is
#' least attainable.
#'
#' Computed through the corrected sum of squares rather than as a ratio, because
#' `n_l * n_u` vanishes at the limits and the ratio form is then `0/0`.
#' Fractional parts are snapped with `tol` so that counts whole up to floating
#' point give exactly 1.
#'
#' @param mean Numeric vector, mean(s) in `[lower, upper]`.
#' @param n Integer scalar, sample size, `n >= 2`.
#' @param lower,upper Numeric scalars, the effective limits.
#' @param tol Tolerance for snapping a near-whole count to whole.
#' @return Numeric vector in `[0, 1]`, the fraction of the Muilwijk ceiling's
#'   squared value that the sharp ceiling attains. 1 where the Muilwijk bound
#'   is itself 0.
#' @seealso [sd_max_muilwijk()] for the uncorrected arch,
#'   [sd_max_structure_s()] for the corrected maximum.
#' @examples
#' # whole counts at the midpoint of a 1-5 scale at n = 7: delta = 6/7
#' sd_delta(mean = 3, n = 7, lower = 1, upper = 5)
#'
#' # close to a limit the ceiling is barely attainable
#' sd_delta(mean = 1.07, n = 7, lower = 1, upper = 5)
#'
#' # it is exactly the gap between the two ceilings strait already exports
#' m <- c(2, 3, 4)
#' sd_max_structure_s(m, 7, 1, 5)^2
#' sd_max_muilwijk(m, 7, 1, 5)^2 * sd_delta(m, 7, 1, 5)
#' @export
sd_delta <- function(mean, n, lower, upper, tol = 1e-9) {
  R <- upper - lower
  nl <- n * (upper - mean) / R
  nu <- n * (mean - lower) / R
  fr <- function(x) {
    f <- x - floor(x)
    ifelse(f > 1 - tol, 0, f)
  }
  muil <- n * (upper - mean) * (mean - lower) # SS of the uncorrected bound
  ss <- muil - fr(nl) * fr(nu) * R^2 # SS of the sharp bound
  ifelse(muil <= 0, 1, pmax(0, ss) / muil)
}
