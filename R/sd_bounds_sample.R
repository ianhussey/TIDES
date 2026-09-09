# Constructive attaining samples for the SD bounds.

#' Construct a sample attaining an SD bound
#'
#' Returns a length-`n` sample that attains the maximum (`which = "max"`) or
#' minimum (`which = "min"`) sample SD of [sd_bounds()] for the given scale
#' limits, sample size and mean, built directly from the closed-form attaining
#' configuration (no search). The returned sample has exactly `mean` as its
#' mean and its SD equals the corresponding bound.
#'
#' The maximum places `n_lower` observations at `l` and `n_upper` at `u` with
#' one remainder observation on the Structure S curve; the minimum clusters
#' `n - 1` observations on the two integers bracketing the mean with one free
#' remainder. Under the quasi-integer default the single remainder may be
#' non-integer; pass `granularity = "integer"` to require a GRIM-consistent
#' mean, so the whole sample is integer.
#'
#' @param l,u Numeric scalars, the scale limits.
#' @param n Integer scalar, sample size (`n >= 2`).
#' @param mean Numeric scalar, the target mean in `[l, u]`.
#' @param which `"max"` (default) or `"min"`.
#' @param granularity `"quasiinteger"` (default, one observation may be
#'   non-integer) or `"integer"` (requires a GRIM-consistent mean).
#' @return A numeric vector of length `n` whose mean is `mean` and whose sample
#'   SD equals the requested bound.
#' @examples
#' x <- sd_bounds_sample(l = 1, u = 7, n = 30, mean = 2.9667, which = "max")
#' c(mean = mean(x), sd = sd(x))
#'
#' table(sd_bounds_sample(l = 1, u = 7, n = 30, mean = 89 / 30, which = "max"))
#' @export
sd_bounds_sample <- function(
  l,
  u,
  n,
  mean,
  which = c("max", "min"),
  granularity = c("quasiinteger", "integer")
) {
  which <- rlang::arg_match(which)
  granularity <- rlang::arg_match(granularity)

  if (is.null(n) || n < 2) {
    cli::cli_abort("{.arg n} must be >= 2.")
  }
  if (mean < l - TOLERANCE || mean > u + TOLERANCE) {
    cli::cli_abort("{.arg mean} must lie in {.code [l, u]}.")
  }
  if (granularity == "integer" && !is_grim_consistent(mean, n)) {
    cli::cli_abort(
      '{.code granularity = "integer"} requires a GRIM-consistent mean \\
       ({.code n * mean} an integer).'
    )
  }

  if (which == "max") {
    span <- u - l
    sum_excess <- n * (mean - l)
    n_upper <- pmin(floor(sum_excess / span + TOLERANCE), n - 1)
    remainder <- l + (sum_excess - n_upper * span)
    n_lower <- n - n_upper - 1

    c(rep(l, n_lower), remainder, rep(u, n_upper))
  } else {
    base_value <- floor(mean)
    # observations at base_value + 1
    n_above <- pmin(floor(n * (mean - base_value) + TOLERANCE), n - 1)
    remainder <- n *
      mean -
      n_above * (base_value + 1) -
      (n - n_above - 1) * base_value

    c(rep(base_value, n - n_above - 1), remainder, rep(base_value + 1, n_above))
  }
}
