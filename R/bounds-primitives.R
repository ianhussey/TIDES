# bounds-primitives.R
#
# The bounds of the sample standard deviation under nested constraints, as a
# family of small single-purpose functions governed by the sd_bounds()
# dispatcher (see sd_bounds.R). The mathematics is the closed-form (not
# constructive) derivation developed in the STRAIT article; each primitive here
# implements one formula.
#
# ARCHITECTURE
#   Layer 0  utilities: bessel_factor, frac, is_grim_consistent, sd_from_ss
#   Layer 1  bound primitives (one formula each):
#            ceilings:  sd_max_span, sd_max_span_n, sd_max_structure_s
#            floors:    sd_min_integer, sd_min_quasi_integer   (range-free)
#                       sd_min_two_pin, sd_min_one_pin         (attained)
#            feasibility: feasible_mean_band
#   Layer 2  alpha composite layer: v_max_alpha, sd_bounds_alpha
#   Layer 3  rounding layer: unround_interval
#
# CONVENTIONS
#   - All SDs are sample SDs (Bessel, n - 1 divisor).
#   - l, u are LOGICAL limits (walls: membership only); a, b are OBSERVED
#     extremes (attained: at least one observation equals each). a supersedes l
#     and b supersedes u, because a wall beyond an attained extreme is vacuous.
#   - granularity is the grid constraint: "continuous" (none), "integer" (all
#     observations on the grid), "quasiinteger" (all but one; GRIM-free).
#   - Sharp mean-conditional ceilings are Structure S in every case: it equals
#     the delta-corrected Muilwijk bound (verified to 1e-13), is valid for
#     continuous data (endpoint counts must be whole regardless), and is
#     unchanged by attainment inside the feasible mean band.

# ---- Layer 0: utilities ------------------------------------------------------

#' Bessel factor
#'
#' The factor `sqrt(n / (n - 1))` converting a population SD of an attaining
#' configuration into the sample SD.
#'
#' @param n Integer scalar or vector, sample size(s), `n >= 2`.
#' @return Numeric, `sqrt(n / (n - 1))`.
#' @keywords internal
bessel_factor <- function(n) sqrt(n / (n - 1))

#' Fractional part, robust to floating-point dust near an integer
#'
#' `x - floor(x)`, but a value within `tol` of an integer (on either side) is
#' snapped to `0`. Scaling a rounded reported mean onto the integer grid can
#' land a hair away from an integer, where a naive fractional part produces a
#' spurious granularity floor `sqrt(d(1 - d))`. The snap is safe: the smallest
#' legitimate on-grid fractional part is `1/n`, orders of magnitude larger
#' than `tol`.
#'
#' @param x Numeric vector.
#' @param tol Numeric tolerance for snapping to the nearest integer.
#' @return `x - floor(x)`, in `[0, 1)`, dust near either integer snapped to 0.
#' @keywords internal
frac <- function(x, tol = TOLERANCE) {
  fractional <- x - floor(x)

  dplyr::if_else(fractional < tol | fractional > 1 - tol, 0, fractional)
}

#' Is a mean GRIM-consistent?
#'
#' A mean is GRIM-consistent at sample size `n` when `n * mean` is an integer,
#' i.e. some sample of `n` integers has exactly that mean.
#'
#' @param mean Numeric scalar, the (exact) mean.
#' @param n Integer scalar, sample size.
#' @param tol Numeric tolerance on `n * mean`'s distance from an integer.
#' @return Logical scalar.
#' @keywords internal
is_grim_consistent <- function(mean, n, tol = TOLERANCE) {
  dplyr::near(n * mean, round(n * mean), tol = tol)
}

#' Sample SD from a sum of squared deviations
#'
#' @param sum_squares Numeric vector, sum of squared deviations about the mean.
#' @param n Integer scalar, sample size.
#' @return `sqrt(pmax(0, sum_squares) / (n - 1))`.
#' @keywords internal
sd_from_ss <- function(sum_squares, n) sqrt(pmax(0, sum_squares) / (n - 1))

# ---- Layer 1: ceilings -------------------------------------------------------

#' Maximum SD given two limits only (no n)
#'
#' The n-free maximum over all sample sizes: attained at `n = 2` with one
#' observation at each limit (Popoviciu's configuration). Valid whether the
#' limits are walls or attained extremes.
#'
#' @param lower,upper Numeric scalars, the effective lower and upper limits.
#' @return `(upper - lower) / sqrt(2)`.
#' @examples
#' sd_max_span(lower = 1, upper = 7)
#' @export
sd_max_span <- function(lower, upper) (upper - lower) / sqrt(2)

#' Maximum SD given two limits and n (parity ceiling)
#'
#' The mean-agnostic maximum for a given sample size: half the observations at
#' each limit, with the odd-n parity correction (Popoviciu 1935 sec 4, Petocz
#' 2005). Unchanged by attainment, since the attaining configuration occupies
#' both limits.
#'
#' @param lower,upper Numeric scalars, the effective limits.
#' @param n Integer scalar, sample size, `n >= 2`.
#' @return Numeric scalar, the maximum sample SD.
#' @examples
#' # an even n splits exactly; an odd n pays the parity correction
#' sd_max_span_n(lower = 1, upper = 7, n = 30)
#' sd_max_span_n(lower = 1, upper = 7, n = 31)
#' @export
sd_max_span_n <- function(lower, upper, n) {
  span <- upper - lower

  if (n %% 2 == 0) {
    (span / 2) * sqrt(n / (n - 1))
  } else {
    (span / 2) * sqrt((n + 1) / n)
  }
}

#' Smooth mean-conditional maximum SD (Muilwijk / Bhatia-Davis)
#'
#' The mean-conditional ceiling
#' `sqrt(n/(n-1)) * sqrt((upper - mean)(mean - lower))`, due to Muilwijk (1966)
#' and rediscovered by Bhatia and Davis (2000). It is the smooth arch obtained
#' by allowing fractional counts at the two limits, so it is *not* sharp: the
#' sharp ceiling applies the count-parity correction and is
#' [sd_max_structure_s()], which never exceeds this.
#'
#' @param mean Numeric vector, mean(s) in `[lower, upper]`.
#' @param n Integer scalar, sample size, `n >= 2`.
#' @param lower,upper Numeric scalars, the limits.
#' @return Numeric vector, the Muilwijk / Bhatia-Davis maximum sample SD.
#' @examples
#' sd_max_muilwijk(mean = c(2, 4, 6), n = 30, lower = 1, upper = 7)
#'
#' # it is not sharp: the parity-corrected ceiling never exceeds it
#' sd_max_structure_s(mean = c(2, 4, 6), n = 30, lower = 1, upper = 7)
#' @export
sd_max_muilwijk <- function(mean, n, lower, upper) {
  bessel_factor(n) * sqrt(pmax(0, (upper - mean) * (mean - lower)))
}

#' Sharp mean-conditional maximum SD (Structure S)
#'
#' The sharp maximum given limits, n, and the mean: `n_upper` observations at
#' the upper limit, `n_lower` at the lower, at most one interior remainder
#' (Mestdagh et al. 2018 "Structure S"). Equal to the count-parity (delta)
#' corrected Muilwijk bound, so it is the sharp ceiling for continuous data too,
#' and it is unchanged by attainment of the limits inside the feasible mean
#' band.
#'
#' @param mean Numeric vector, mean(s) in `[lower, upper]`.
#' @param n Integer scalar, sample size, `n >= 2`.
#' @param lower,upper Numeric scalars, the effective limits.
#' @return Numeric vector, the maximum sample SD at each mean.
#' @examples
#' sd_max_structure_s(mean = c(2, 4, 6), n = 30, lower = 1, upper = 7)
#'
#' # a mean at a limit forces every observation there, so the SD is 0
#' sd_max_structure_s(mean = 1, n = 30, lower = 1, upper = 7)
#' @export
sd_max_structure_s <- function(mean, n, lower, upper) {
  span <- upper - lower
  sum_excess <- n * (mean - lower)
  n_upper <- pmin(floor(sum_excess / span + TOLERANCE), n - 1)
  remainder <- lower + (sum_excess - n_upper * span)
  n_lower <- n - n_upper - 1
  sum_squares <- n_lower *
    lower^2 +
    remainder^2 +
    n_upper * upper^2 -
    n * mean^2

  sd_from_ss(sum_squares, n)
}

# ---- Layer 1: range-free floors (granularity only) ---------------------------

#' Minimum SD for strictly integer data (Bernoulli floor)
#'
#' All observations on the two integers adjacent to the mean (Pesant & Regin
#' 2005; Bernoulli form). Range-free. Defined only at GRIM-consistent means;
#' the caller is responsible for that check (see [is_grim_consistent()]),
#' otherwise use [sd_min_quasi_integer()].
#'
#' @param mean Numeric vector, GRIM-consistent mean(s).
#' @param n Integer scalar, sample size, `n >= 2`.
#' @return Numeric vector, the minimum sample SD at each mean.
#' @examples
#' # zero at a whole-number mean, largest half way between two integers
#' sd_min_integer(mean = 3, n = 30)
#' sd_min_integer(mean = 3.5, n = 30)
#' @export
sd_min_integer <- function(mean, n) {
  fractional <- frac(mean)

  bessel_factor(n) * sqrt(fractional * (1 - fractional))
}

#' Minimum SD for quasi-integer data (GRIM-free floor)
#'
#' All but one observation an integer. Defined at every real mean; coincides
#' with [sd_min_integer()] at GRIM-consistent means and dips between them, so it
#' is a valid floor under either the strict or the relaxed hypothesis.
#' Range-free.
#'
#' @param mean Numeric vector, mean(s).
#' @param n Integer scalar, sample size, `n >= 2`.
#' @return Numeric vector, the minimum sample SD at each mean.
#' @examples
#' # at a GRIM-consistent mean it agrees with the strict-integer floor
#' sd_min_quasi_integer(mean = 3.5, n = 30)
#' sd_min_integer(mean = 3.5, n = 30)
#'
#' # between GRIM means it dips below, so it is defined at every mean
#' sd_min_quasi_integer(mean = 3.51, n = 30)
#' @export
sd_min_quasi_integer <- function(mean, n) {
  fractional <- frac(mean)
  fractional_sum <- frac(n * mean)

  sd_from_ss(
    n * fractional * (1 - fractional) - fractional_sum * (1 - fractional_sum),
    n
  )
}

# ---- Layer 1: attained-extremes floors ---------------------------------------

#' Minimum SD with both observed extremes attained (two-pin floor)
#'
#' One observation pinned at each observed extreme. Without a mean this is the
#' Nagy (1918) / Thomson (1955) floor `W / sqrt(2 (n - 1))`. With a mean it is
#' the pinned-plus-interior decomposition of the attained-extremes analysis,
#' plus, for integer data, the Pesant-Regin clustering of the `n - 2` interior
#' observations. The mean must lie in the feasible band (see
#' [feasible_mean_band()]); this function does not check it.
#'
#' @param a,b Numeric scalars, observed minimum and maximum (attained).
#' @param n Integer scalar, sample size, `n >= 2`.
#' @param mean Numeric scalar or NULL. If NULL, the mean-agnostic floor.
#' @param granularity One of "continuous", "integer", "quasiinteger".
#' @return Numeric scalar, the minimum sample SD.
#' @examples
#' # mean-agnostic floor: one observation pinned at each observed extreme
#' sd_min_two_pin(a = 1, b = 7, n = 30)
#'
#' # knowing the mean sharpens it, and granularity sharpens it again
#' sd_min_two_pin(a = 1, b = 7, n = 30, mean = 3)
#' sd_min_two_pin(a = 1, b = 7, n = 30, mean = 3, granularity = "integer")
#' @export
sd_min_two_pin <- function(a, b, n, mean = NULL, granularity = "continuous") {
  span <- b - a
  if (is.null(mean)) {
    return(span / sqrt(2 * (n - 1)))
  }
  if (n == 2) {
    # the feasible band is then the single midpoint mean
    return(span / sqrt(2))
  }

  n_interior <- n - 2
  dist_lower <- mean - a
  dist_upper <- b - mean
  sum_squares_continuous <- dist_lower^2 +
    dist_upper^2 +
    (dist_lower - dist_upper)^2 / n_interior
  sum_squares <- sum_squares_continuous

  if (granularity %in% c("integer", "quasiinteger")) {
    sum_interior <- n * mean - a - b
    fractional <- frac(sum_interior / n_interior)
    sum_squares <- sum_squares_continuous +
      n_interior * fractional * (1 - fractional)
    if (granularity == "quasiinteger") {
      sum_squares <- sum_squares -
        frac(sum_interior) * (1 - frac(sum_interior))
    }
    sum_squares <- max(sum_squares, sum_squares_continuous)
  }

  sd_from_ss(sum_squares, n)
}

#' Minimum SD with one observed extreme attained (one-pin floor)
#'
#' One observation pinned at a single attained extreme, the other side at most
#' walled. NOT the two-pin formula with a pin deleted: the interior count is
#' `n - 1`, not `n - 2`. Continuous floor `SS = n q^2 / (n - 1)` (strictly
#' sharper than the Laguerre-Samuelson corollary), plus the same granularity
#' terms on the `n - 1` free observations. Requires a mean; without one the
#' floor is zero, since all observations may sit at the pin.
#'
#' @param pin Numeric scalar, the attained extreme.
#' @param n Integer scalar, sample size, `n >= 2`.
#' @param mean Numeric scalar.
#' @param granularity One of "continuous", "integer", "quasiinteger".
#' @param side "max" if the pin is the observed maximum, "min" if the minimum.
#' @return Numeric scalar, the minimum sample SD.
#' @examples
#' # an observed maximum of 7 alongside a mean of 3 forces some spread
#' sd_min_one_pin(pin = 7, n = 30, mean = 3)
#'
#' # the mirror case: an observed minimum of 1 alongside a mean of 5
#' sd_min_one_pin(pin = 1, n = 30, mean = 5, side = "min")
#' @export
sd_min_one_pin <- function(
  pin,
  n,
  mean,
  granularity = "continuous",
  side = c("max", "min")
) {
  side <- rlang::arg_match(side)
  if (side == "min") {
    # reflect the observed-minimum case onto the observed-maximum formulas
    return(sd_min_one_pin(-pin, n, -mean, granularity, side = "max"))
  }

  n_free <- n - 1
  dist_to_pin <- pin - mean
  sum_squares_continuous <- n * dist_to_pin^2 / (n - 1)
  sum_squares <- sum_squares_continuous

  if (granularity %in% c("integer", "quasiinteger")) {
    sum_free <- n * mean - pin
    fractional <- frac(sum_free / n_free)
    sum_squares <- sum_squares_continuous +
      n_free * fractional * (1 - fractional)
    if (granularity == "quasiinteger") {
      sum_squares <- sum_squares - frac(sum_free) * (1 - frac(sum_free))
    }
    sum_squares <- max(sum_squares, sum_squares_continuous)
  }

  sd_from_ss(sum_squares, n)
}

# ---- Layer 1: feasibility ----------------------------------------------------

#' Feasible mean interval given the side constraints
#'
#' Walls admit any mean in `[lower, upper]`. Each attained extreme pins one
#' observation and narrows the band: both attained gives `[a + W/n, b - W/n]`;
#' an attained maximum with (at most) a walled minimum gives
#' `[lower + (b - lower)/n, b]`, and the mirror for an attained minimum.
#' `-Inf`/`Inf` when a side is absent.
#'
#' @param lower,upper Numeric scalars or NULL, the effective limits.
#' @param lower_attained,upper_attained Logical, is each limit attained?
#' @param n Integer scalar or NULL.
#' @return Numeric length-2 vector, the closed feasible mean interval.
#' @examples
#' # walls alone admit any mean on the scale
#' feasible_mean_band(lower = 1, upper = 7, n = 30)
#'
#' # an attained maximum pins one observation at 7, lifting the lowest mean
#' feasible_mean_band(lower = 1, upper = 7, upper_attained = TRUE, n = 30)
#'
#' # both extremes attained narrows the band at each end
#' feasible_mean_band(
#'   lower = 1, upper = 7, lower_attained = TRUE,
#'   upper_attained = TRUE, n = 30
#' )
#' @export
feasible_mean_band <- function(
  lower = NULL,
  upper = NULL,
  lower_attained = FALSE,
  upper_attained = FALSE,
  n = NULL
) {
  band_lo <- if (is.null(lower)) -Inf else lower
  band_hi <- if (is.null(upper)) Inf else upper

  if (!is.null(n)) {
    if (lower_attained && upper_attained) {
      span <- upper - lower
      return(c(lower + span / n, upper - span / n))
    }
    if (upper_attained && is.finite(band_lo)) {
      band_lo <- max(band_lo, lower + (upper - lower) / n)
    }
    if (lower_attained && is.finite(band_hi)) {
      band_hi <- min(band_hi, upper - (upper - lower) / n)
    }
  }

  c(band_lo, band_hi)
}

# ---- Layer 2: alpha composite layer ------------------------------------------

#' Allocation maximum of the item-variance sum (Theorem H3)
#'
#' The quasi-integer maximum of the summed item variances of a k-item integer
#' battery at a given sum-score mean (the STRAIT article, Appendix H).
#'
#' @param mean_sum Numeric vector, sum-score mean(s).
#' @param k Integer scalar, number of items.
#' @param n Integer scalar, sample size.
#' @param item_l,item_u Numeric scalars, the per-item limits.
#' @return Numeric vector, the maximum of the population item-variance sum.
#' @examples
#' # a 3-item 1-5 battery (sum score 3-15) at a sum-score mean of 9
#' v_max_alpha(mean_sum = 9, k = 3, n = 50, item_l = 1, item_u = 5)
#' @keywords internal
#' @export
v_max_alpha <- function(mean_sum, k, n, item_l, item_u) {
  step <- (item_u - item_l) / n
  total <- (mean_sum - k * item_l) / step
  per_item <- total / k
  floor_count <- pmin(floor(per_item + TOLERANCE_TIGHT), n - 1)
  theta <- per_item - floor_count
  phi <- total - floor(total + TOLERANCE_TIGHT)
  phi <- dplyr::if_else(phi > 1 - TOLERANCE, 0, phi)
  allocation <- function(count) step^2 * count * (n - count)

  pmax(
    0,
    k *
      ((1 - theta) *
        allocation(floor_count) +
        theta * allocation(floor_count + 1)) -
      step^2 * (n - 1) * phi * (1 - phi)
  )
}

# Internal: all non-negative integer count vectors (c_0, ..., c_{slots-1}) that
# sum to `total`. Rows enumerate compositions; used for the exact Gini envelope.
count_vectors <- function(slots, total) {
  if (slots == 1L) {
    return(matrix(total, ncol = 1))
  }

  # rbind() rather than a tidyverse binder: these are matrices, and building
  # them as data frames would cost more than the enumeration itself.
  do.call(
    rbind,
    purrr::map(
      0:total,
      \(count) cbind(count, count_vectors(slots - 1L, total - count))
    )
  )
}

#' Sharp alpha-conditional composite floor via the Gini mean difference
#'
#' The exact alpha-conditional minimum SD for a *strictly integer* composite
#' (STRAIT article, Theorem H5 and Corollaries H5a-H5b): the least sample SD
#' among integer sum-score profiles \eqn{S} with the given mean, integer values
#' in \eqn{[l, u]}, whose design-factor cap
#' \eqn{m_{\max}(S) = 2 n\, SS_S / \sum_{s,t}|S_s - S_t|} is at least the
#' reported design factor \eqn{m = 1/(1 - c\alpha)}. Unlike the amplified floor
#' it is strictly positive whenever \eqn{m > 1}, including at whole-number
#' sum-score means (Corollary H5b).
#'
#' Returns `NULL` (caller falls back to the proven amplified floor) when no
#' integer composite exists (\eqn{n \cdot mean} not an integer) or the
#' enumeration `choose(n + W, W)`, \eqn{W = u - l}, would exceed `max_profiles`.
#'
#' @param l,u Integer scalars, sum-score limits.
#' @param n Integer scalar, sample size.
#' @param mean Numeric scalar, sum-score mean (with \eqn{n \cdot mean} an
#'   integer).
#' @param m Numeric scalar, the reported design factor
#'   \eqn{1/(1 - c\alpha)}, \eqn{\ge 1}.
#' @param max_profiles Numeric, enumeration budget.
#' @return Numeric scalar, the exact floor in sum-score SD units, or `NULL`.
#' @keywords internal
sd_min_alpha_gini <- function(
  l,
  u,
  n,
  mean,
  m,
  max_profiles = MAX_GINI_PROFILES
) {
  sum_total <- round(n * mean)
  if (!dplyr::near(n * mean, sum_total, tol = TOLERANCE)) {
    return(NULL) # no integer composite
  }

  target <- sum_total - n * l # required sum of (value - l)
  envelope <- alpha_gini_envelope(l, u, n, m, max_profiles)
  if (is.null(envelope)) {
    return(NULL) # over budget, fall back
  }

  index <- match(target, envelope$target)
  if (is.na(index)) {
    return(NULL)
  }

  sqrt(envelope$sum_squares[index] / (n - 1))
}

# Internal: the whole Gini envelope in ONE enumeration pass, memoized.
# Grouping by sum yields every mean's floor from a single pass, which is what
# makes mean-sweeps (figures, umbrella grids) affordable; filtering a full
# enumeration down to one sum would re-enumerate the same 10^5-10^6 profiles
# once per mean. Cached on (l, u, n, m, budget) because `m` decides which
# profiles can support the reported alpha.
gini_envelope_cache <- new.env(parent = emptyenv())

alpha_gini_envelope <- function(
  l,
  u,
  n,
  m,
  max_profiles = MAX_GINI_PROFILES
) {
  key <- paste(l, u, n, signif(m, 12), max_profiles, sep = "|")
  hit <- gini_envelope_cache[[key]]
  if (!is.null(hit)) {
    return(if (identical(hit, NA)) NULL else hit)
  }

  spread <- as.integer(round(u - l))
  over_budget <- spread < 1L ||
    !is.finite(suppressWarnings(choose(n + spread, spread))) ||
    choose(n + spread, spread) > max_profiles
  if (over_budget) {
    assign(key, NA, envir = gini_envelope_cache)
    return(NULL)
  }

  values <- 0:spread # deviations above l
  counts <- count_vectors(spread + 1L, n)
  totals <- as.vector(counts %*% values)
  # SS_S, shift-invariant
  sum_squares <- as.vector(counts %*% (values^2)) - totals^2 / n

  # Gini double-sum over unordered value pairs:
  # sum_{a<b} c_a c_b (v_b - v_a)
  gini_sum <- numeric(nrow(counts))
  for (lower_current in seq_len(spread)) {
    for (upper_current in (lower_current + 1L):(spread + 1L)) {
      gini_sum <- gini_sum +
        counts[, lower_current] *
          counts[, upper_current] *
          (values[upper_current] - values[lower_current])
    }
  }

  # m_max(S) = n SS / (n V_min)
  design_factor_max <- dplyr::if_else(
    gini_sum > 0,
    n * sum_squares / gini_sum,
    0
  )
  # non-constant, and able to support the reported alpha
  supports_alpha <- design_factor_max >= m - TOLERANCE &
    sum_squares > TOLERANCE_TIGHT

  out <- if (!any(supports_alpha)) {
    list(target = integer(0), sum_squares = numeric(0))
  } else {
    best <- tapply(sum_squares[supports_alpha], totals[supports_alpha], min)
    list(target = as.integer(names(best)), sum_squares = as.numeric(best))
  }
  assign(key, out, envir = gini_envelope_cache)

  out
}

#' SD bounds for a k-item composite with reported Cronbach's alpha
#'
#' Bounds of the sum score's sample SD given the composite limits, n, sum-score
#' mean, granularity, and a reported alpha (the STRAIT article, Appendix H).
#' Ceiling: the mean-conditional ceiling divided by `sqrt(k - (k - 1) alpha)`;
#' under granularity, additionally the sharper quasi-integer `V_max` form.
#' Floor under granularity: the alpha-amplified quasi-integer floor, sharpened
#' for strictly integer composites to the exact Gini envelope of
#' [sd_min_alpha_gini()] when that enumeration is affordable. The floor can
#' exceed the ceiling for small k (a genuinely infeasible constraint set); this
#' is reported rather than clipped.
#'
#' @param l,u Numeric scalars, limits of the SUM score (composite units).
#' @param n Integer scalar, sample size.
#' @param mean Numeric scalar, sum-score mean.
#' @param granularity One of "continuous", "integer", "quasiinteger"
#'   (item-level grid).
#' @param alpha Numeric scalar, reported Cronbach's alpha, `alpha < 1`.
#' @param k_items Integer scalar, number of items, `k_items >= 1`.
#' @return List with `min_sd`, `max_sd`, `feasible`, `min_rule`, `note`.
#' @examples
#' # a 3-item 1-5 battery, so the sum score runs 3-15
#' sd_bounds_alpha(
#'   l = 3, u = 15, n = 50, mean = 9.4, granularity = "continuous",
#'   alpha = 0.8, k_items = 3
#' )
#'
#' # integer items lift the floor off zero and tighten the ceiling
#' sd_bounds_alpha(
#'   l = 3, u = 15, n = 50, mean = 9.4, granularity = "integer",
#'   alpha = 0.8, k_items = 3
#' )
#' @export
sd_bounds_alpha <- function(l, u, n, mean, granularity, alpha, k_items) {
  alpha_weight <- (k_items - 1) / k_items
  # k - (k - 1) alpha = k * design_divisor
  design_divisor <- 1 - alpha_weight * alpha
  if (design_divisor <= TOLERANCE_TIGHT) {
    cli::cli_abort(
      "{.arg alpha} is too high for {.arg k_items} = {k_items}: \\
       {.code 1 - alpha * (k_items - 1) / k_items} must be positive."
    )
  }

  ceiling_free <- sd_max_structure_s(mean, n, l, u) # alpha-free sharp ceiling
  ceiling_smooth <- bessel_factor(n) *
    sqrt(pmax(0, (u - mean) * (mean - l)) / (k_items * design_divisor))
  ceiling_sd <- min(ceiling_free, ceiling_smooth)
  floor_sd <- 0
  floor_rule <- "s >= 0"
  note <- NA_character_

  if (granularity %in% c("integer", "quasiinteger")) {
    ceiling_vmax <- sqrt(
      (n / (n - 1)) *
        v_max_alpha(mean, k_items, n, l / k_items, u / k_items) /
        design_divisor
    )
    ceiling_sd <- min(ceiling_sd, ceiling_vmax)
    # proven amplified floor
    floor_sd <- sd_min_quasi_integer(mean, n) / sqrt(design_divisor)
    floor_rule <- "alpha-amplified quasi-integer floor"

    # sharpen with the exact Gini envelope for strictly integer composites; it
    # is >= the amplified floor and strictly positive at whole-number means.
    if (granularity == "integer" && alpha > 0) {
      exact_floor <- sd_min_alpha_gini(
        l,
        u,
        n,
        mean,
        m = 1 / design_divisor
      )
      if (is.null(exact_floor)) {
        note <- paste(
          "exact alpha floor (Theorem H5) not evaluated (n*mean non-integer",
          "or composite window over budget); proven amplified floor used,",
          "conservative near whole-number means"
        )
      } else if (exact_floor > floor_sd + TOLERANCE_TIGHT) {
        floor_sd <- exact_floor
        floor_rule <- "alpha-conditional Gini envelope floor (Theorem H5)"
      }
    }
  }

  list(
    min_sd = floor_sd,
    max_sd = ceiling_sd,
    feasible = floor_sd <= ceiling_sd + TOLERANCE,
    min_rule = floor_rule,
    note = note
  )
}

# ---- Layer 3: rounding / truncation of reported inputs -----------------------

#' Reconstruct the interval of exact values behind a rounded/truncated report
#'
#' Given a reported value and the rounding procedure that produced it, returns
#' the interval of exact values consistent with the report. With
#' `unit = 10^-digits` and `h = unit / 2`:
#'
#' * `"up_or_down"` (default): `[x - h, x + h]`, both endpoints included — the
#'   agnostic choice when the direction at ties is unknown.
#' * `"up"` (round half up): `[x - h, x + h)`.
#' * `"down"` (round half down): `(x - h, x + h]`.
#' * `"even"` (banker's): `[x - h, x + h]`, both included.
#' * `"ceiling"`: `(x - unit, x]`; `"floor"`: `[x, x + unit)`.
#' * `"trunc"` (toward zero): `[x, x + unit)` for `x >= 0`, `(x - unit, x]` for
#'   `x < 0`. `"anti_trunc"` (away from zero) is its mirror.
#'
#' The interval itself comes from [roundwork::unround()], so this package and
#' the rest of the error-detection ecosystem unround a report the same way.
#' This function is the adapter: it names the endpoints `lo`/`hi`/`lo_incl`/
#' `hi_incl`, returns a plain list, and carries `digits` back out.
#'
#' A numeric `x` needs `digits` supplied explicitly, since trailing zeros
#' cannot be recovered from it. A character `x` is accepted as a convenience,
#' with digits inferred via [scrutiny::decimal_places()].
#'
#' @param x Numeric scalar, the reported value (character accepted for digit
#'   inference only).
#' @param digits Integer, reported decimal places. Required when `x` is
#'   numeric.
#' @param rounding One of the options above.
#' @return List: `lo`, `hi` (numeric), `lo_incl`, `hi_incl` (logical), `digits`
#'   (integer used).
#' @examples
#' # the exact values that could have been reported as 2.97
#' unround_interval(2.97, digits = 2)
#'
#' # truncation rather than rounding shifts the interval
#' unround_interval(2.97, digits = 2, rounding = "trunc")
#'
#' # digits are inferred from a string, so trailing zeros are honoured
#' unround_interval("2.90")
#' @export
unround_interval <- function(
  x,
  digits = NULL,
  rounding = c(
    "up_or_down",
    "up",
    "down",
    "even",
    "ceiling",
    "floor",
    "trunc",
    "anti_trunc"
  )
) {
  rounding <- rlang::arg_match(rounding)
  if (is.null(digits)) {
    if (!is.character(x)) {
      cli::cli_abort(c(
        "{.arg digits} must be supplied for a numeric {.arg x}.",
        i = "Trailing zeros are not recoverable from a numeric.",
        i = "Or pass the reported value as a string, to infer them."
      ))
    }
    digits <- scrutiny::decimal_places(x)
  }

  interval <- roundwork::unround(x, rounding = rounding, digits = digits)
  lo_incl <- interval$incl_lower
  hi_incl <- interval$incl_upper

  # roundwork is still taking over scrutiny's rounding layer, and its unround()
  # does not yet agree with its own rounding functions in two rules. Both are
  # patched here rather than worked around at each call site; drop this block
  # once roundwork settles, and the verdicts will not move.
  #
  #  * "even" returns NA for both flags, since whether a tie rounds to x
  #    depends on the parity of x's last digit. Reading NA as included is the
  #    conservative direction -- it admits more exact values behind the report,
  #    never fewer -- and it is what scrutiny::unround() returns today.
  #  * "anti_trunc" returns them the wrong way round relative to
  #    roundwork::round_anti_trunc(), which maps (2.96, 2.97] to 2.97: the
  #    lower endpoint is excluded and the upper included, not the mirror.
  if (rounding == "even") {
    lo_incl <- TRUE
    hi_incl <- TRUE
  } else if (rounding == "anti_trunc") {
    swap <- lo_incl
    lo_incl <- hi_incl
    hi_incl <- swap
  }

  list(
    lo = interval$lower,
    hi = interval$upper,
    lo_incl = lo_incl,
    hi_incl = hi_incl,
    digits = digits
  )
}
