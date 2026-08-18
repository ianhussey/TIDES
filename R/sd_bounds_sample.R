# Constructive attaining samples for the SD bounds.

#' Construct a sample attaining an SD bound
#'
#' Returns a length-`n` sample that attains the maximum (`which = "max"`) or
#' minimum (`which = "min"`) sample SD of [sd_bounds()] for the given scale
#' limits, sample size and mean, built directly from the closed-form attaining
#' configuration (no search). This is the constructive companion to the
#' closed-form [sd_bounds()]: the returned sample has exactly `mean` as its mean
#' and its SD equals the corresponding bound.
#'
#' The maximum places `n_l` observations at `l` and `n_u` at `u` with one
#' remainder observation on the Structure S curve; the minimum clusters `n - 1`
#' observations on the two integers bracketing the mean with one free remainder.
#' Under the quasi-integer default the single remainder may be non-integer (the
#' `Z_{n-1}` relaxation used throughout the package); pass `Z = "integer"` to
#' require a GRIM-consistent mean, so the whole sample is integer.
#'
#' @param l,u Numeric scalars, the scale limits.
#' @param n Integer scalar, sample size (`n >= 2`).
#' @param mean Numeric scalar, the target mean in `[l, u]`.
#' @param which `"max"` (default) or `"min"`.
#' @param Z `"quasiinteger"` (default, one observation may be non-integer) or
#'   `"integer"` (requires a GRIM-consistent mean).
#' @return A numeric vector of length `n` whose mean is `mean` and whose sample
#'   SD equals the requested bound.
#' @examples
#' x <- sd_bounds_sample(l = 1, u = 7, n = 30, mean = 2.9667, which = "max")
#' c(mean = mean(x), sd = sd(x))          # 2.9667 and the max-SD bound
#' table(sd_bounds_sample(l = 1, u = 7, n = 30, mean = 89/30, which = "max"))
#' @export
sd_bounds_sample <- function(
  l,
  u,
  n,
  mean,
  which = c("max", "min"),
  Z = c("quasiinteger", "integer")
) {
  which <- match.arg(which)
  Z <- match.arg(Z)
  if (is.null(n) || n < 2) {
    stop("n must be >= 2")
  }
  if (mean < l - 1e-9 || mean > u + 1e-9) {
    stop("mean must lie in [l, u]")
  }
  if (Z == "integer" && abs(n * mean - round(n * mean)) > 1e-9) {
    stop("Z = 'integer' requires a GRIM-consistent mean (n * mean an integer)")
  }

  if (which == "max") {
    R <- u - l
    S_excess <- n * (mean - l)
    n_u <- pmin(floor(S_excess / R + 1e-9), n - 1)
    x_r <- l + (S_excess - n_u * R)
    n_l <- n - n_u - 1
    c(rep(l, n_l), x_r, rep(u, n_u))
  } else {
    v <- floor(mean)
    n_hi <- pmin(floor(n * (mean - v) + 1e-9), n - 1) # observations at v + 1
    x_r <- n * mean - n_hi * (v + 1) - (n - n_hi - 1) * v # the free remainder
    c(rep(v, n - n_hi - 1), x_r, rep(v + 1, n_hi))
  }
}
