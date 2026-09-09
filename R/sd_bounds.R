# sd_bounds() dispatcher and its exact-mean core, plus the POMP transforms.
# Part of the strait package; see bounds-primitives.R for the bound formulas.

# ---- Layer 4: dispatcher -----------------------------------------------------

# Internal: bounds for one EXACT mean (or no mean), with sides already
# resolved. Returns a list; feasibility gates for band/GRIM are applied here.
sd_bounds_core <- function(
  lower,
  upper,
  lower_attained,
  upper_attained,
  n,
  mean,
  granularity,
  alpha,
  k_items
) {
  out <- list(
    min_sd = 0,
    max_sd = Inf,
    feasible = TRUE,
    min_rule = "s >= 0",
    max_rule = "unbounded",
    note = NA_character_
  )
  fail <- function(note) {
    list(
      min_sd = NA_real_,
      max_sd = NA_real_,
      feasible = FALSE,
      min_rule = NA_character_,
      max_rule = NA_character_,
      note = note
    )
  }

  # feasibility: mean inside the band implied by walls/pins
  if (!is.null(mean) && (!is.null(lower) || !is.null(upper))) {
    band <- feasible_mean_band(
      lower,
      upper,
      lower_attained,
      upper_attained,
      n
    )
    if (mean < band[1] - TOLERANCE || mean > band[2] + TOLERANCE) {
      return(fail(
        if (lower_attained || upper_attained) {
          sprintf(
            "mean outside the feasible band [%.6g, %.6g] implied by the attained extreme(s)",
            band[1],
            band[2]
          )
        } else {
          "mean outside [l, u]"
        }
      ))
    }
  }

  # feasibility: GRIM under strict integer
  if (
    identical(granularity, "integer") &&
      !is.null(mean) &&
      !is_grim_consistent(mean, n)
  ) {
    return(fail(
      "mean is GRIM-inconsistent: no strictly integer sample has this mean"
    ))
  }

  # alpha branch (walls only; enforced by the dispatcher)
  if (!is.null(alpha)) {
    bounds_alpha <- sd_bounds_alpha(
      lower,
      upper,
      n,
      mean,
      granularity,
      alpha,
      k_items
    )
    if (!bounds_alpha$feasible) {
      return(fail(
        "alpha floor exceeds alpha ceiling: no sample satisfies all constraints"
      ))
    }
    out$min_sd <- bounds_alpha$min_sd
    out$max_sd <- bounds_alpha$max_sd
    out$min_rule <- if (bounds_alpha$min_sd > 0) {
      bounds_alpha$min_rule
    } else {
      "s >= 0"
    }
    out$max_rule <- "alpha ceiling"
    if (!is.na(bounds_alpha$note)) {
      out$note <- bounds_alpha$note
    }
    return(out)
  }

  # ceiling: min over applicable ceilings
  if (!is.null(lower) && !is.null(upper)) {
    if (is.null(n)) {
      out$max_sd <- sd_max_span(lower, upper)
      out$max_rule <- "span/sqrt(2) (n = 2 maximum)"
    } else if (is.null(mean)) {
      out$max_sd <- sd_max_span_n(lower, upper, n)
      out$max_rule <- "parity ceiling"
    } else {
      out$max_sd <- sd_max_structure_s(mean, n, lower, upper)
      out$max_rule <- "Structure S (sharp mean-conditional ceiling)"
    }
  }

  # floor: max over applicable floors
  bump <- function(candidate, rule) {
    if (candidate > out$min_sd + TOLERANCE_TIGHT) {
      out$min_sd <<- candidate
      out$min_rule <<- rule
    }
  }
  if (granularity %in% c("integer", "quasiinteger") && !is.null(mean)) {
    floor_sd <- if (granularity == "integer") {
      sd_min_integer(mean, n)
    } else {
      sd_min_quasi_integer(mean, n)
    }
    bump(floor_sd, sprintf("%s floor (range-free)", granularity))
  }
  if (!is.null(n)) {
    if (lower_attained && upper_attained) {
      bump(
        sd_min_two_pin(lower, upper, n, mean, granularity),
        if (is.null(mean)) {
          "two-pin floor W/sqrt(2(n-1))"
        } else {
          "two-pin attained floor"
        }
      )
    } else if (upper_attained && !is.null(mean)) {
      bump(
        sd_min_one_pin(upper, n, mean, granularity, side = "max"),
        "one-pin attained floor (max)"
      )
    } else if (lower_attained && !is.null(mean)) {
      bump(
        sd_min_one_pin(lower, n, mean, granularity, side = "min"),
        "one-pin attained floor (min)"
      )
    }
  }

  if (out$min_sd > out$max_sd + TOLERANCE) {
    return(fail("floor exceeds ceiling: no sample satisfies all constraints"))
  }

  out
}

# Internal: candidate exact means inside [m_lo, m_hi] for the envelope -- a
# dense grid plus every analytic breakpoint of the piecewise bounds (multiples
# of 1/n, where both floors and Structure S can kink).
candidate_means <- function(
  mean_lo,
  mean_hi,
  n,
  grid_n = ENVELOPE_GRID_POINTS
) {
  candidates <- seq(mean_lo, mean_hi, length.out = grid_n)

  if (!is.null(n)) {
    sums <- ceiling(n * mean_lo - TOLERANCE):floor(n * mean_hi + TOLERANCE)
    if (length(sums) && length(sums) <= 5000) {
      candidates <- c(
        candidates,
        sums / n,
        pmax(mean_lo, pmin(mean_hi, sums / n + TOLERANCE_TIGHT)),
        pmax(mean_lo, pmin(mean_hi, sums / n - TOLERANCE_TIGHT))
      )
    }
  }

  sort(unique(pmax(mean_lo, pmin(mean_hi, candidates))))
}

#' Bounds of the sample SD under a chosen set of constraints
#'
#' Computes the smallest and largest sample standard deviations consistent with
#' whichever constraints are supplied, following the nested framework and the
#' attained-extremes extensions of the STRAIT article. All constraints default
#' to `NULL`; supplying more of them can only narrow the bounds. With none
#' supplied the only bound is `0 <= s < Inf`.
#'
#' Constraint semantics:
#' * `l`, `u` — logical scale limits (walls): observations lie in `[l, u]` but
#'   need not touch either limit.
#' * `a`, `b` — observed extremes (attained): at least one observation equals
#'   each. `a` supersedes `l` and `b` supersedes `u`. Attainment leaves
#'   ceilings unchanged but creates nonzero floors.
#' * `mean`, `digits_mean` — the reported mean with its number of reported
#'   decimal places, mirroring scrutiny's numeric-plus-digits API. With
#'   `rounding = NULL` the mean is treated as exact and `digits_mean` is
#'   ignored.
#' * `rounding` — `NULL` (mean exact) or a rule from [unround_interval()]. The
#'   returned bounds are then the envelope over all exact means consistent with
#'   the report, intersected with the feasible mean band.
#' * `sd`, `digits_sd` — the reported SD never changes the bounds; it enables
#'   `sd_in_bounds` and, under `granularity = "integer"` with `rounding`, the
#'   GRIMMER verdict.
#' * `granularity` — `"continuous"` (none; the default), `"integer"` (all
#'   observations on the grid), or `"quasiinteger"` (all but one; defined at
#'   every mean, the GRIM-free relaxation).
#' * `scoring`, `n_items` — `"singleitem"` (default, requires `n_items = 1`);
#'   `"sumscored"` (observations are sums of `n_items` integer items, still
#'   integer-gridded); `"meanscored"` (observations sit on a `1/n_items` grid).
#' * `alpha` — reported Cronbach's alpha; requires `scoring = "sumscored"` or
#'   `"meanscored"` plus `l`, `u`, `mean`, `n`. Not combinable with `a`/`b`.
#'
#' GRIM and GRIMMER verdicts are deferred to [scrutiny::grim()] and
#' [scrutiny::grimmer()] rather than reimplemented. The bounds themselves come
#' from this package's own enumeration of GRIM-consistent means; a disagreement
#' between the two is surfaced in `note` rather than resolved, since scrutiny
#' has documented floating-point boundary cases.
#'
#' @param l,u Numeric scalars or NULL, logical limits.
#' @param a,b Numeric scalars or NULL, observed (attained) extremes.
#' @param n Integer scalar or NULL, sample size (`n >= 2` when given).
#' @param mean Numeric scalar or NULL, the reported sample mean.
#' @param digits_mean Integer scalar or NULL, decimal places of the reported
#'   mean (required with `rounding`).
#' @param sd Numeric scalar or NULL, the reported sample SD.
#' @param digits_sd Integer scalar or NULL, decimal places of the reported SD.
#' @param rounding NULL (mean is exact) or a rounding rule; see
#'   [unround_interval()].
#' @param granularity "continuous" (default), "integer", or "quasiinteger".
#' @param scoring "singleitem" (default), "sumscored", or "meanscored".
#' @param n_items Positive whole number of response items (default 1).
#' @param alpha Numeric scalar or NULL, reported Cronbach's alpha.
#'
#' @return A one-row tibble: `min_sd`, `max_sd` (`Inf` if unbounded),
#'   `feasible` (FALSE when the constraint set admits no sample, with
#'   `min_sd`/`max_sd` `NA`), `min_rule`, `max_rule`, `grim`, `grimmer`
#'   (scrutiny's verdicts; `NA` when not applicable), `sd_in_bounds`, `note`.
#'
#' @examples
#' sd_bounds() # 0 <= s < Inf
#' sd_bounds(l = 1, u = 5, n = 9, mean = 2) # sharp mean-conditional
#' sd_bounds(a = 2, b = 4, n = 9, mean = 2.44, granularity = "integer")
#'
#' # a mean rounded to 2 dp: envelope over [2.965, 2.975]
#' sd_bounds(
#'   l = 1, u = 7, n = 30, mean = 2.97, digits_mean = 2,
#'   rounding = "up_or_down", granularity = "integer"
#' )
#'
#' # sum-scored 3-item composite (l, u, mean in sum units) with alpha
#' sd_bounds(
#'   l = 3, u = 15, n = 20, mean = 9, granularity = "integer",
#'   scoring = "sumscored", n_items = 3, alpha = 0.8
#' )
#' @export
sd_bounds <- function(
  l = NULL,
  u = NULL,
  a = NULL,
  b = NULL,
  n = NULL,
  mean = NULL,
  digits_mean = NULL,
  sd = NULL,
  digits_sd = NULL,
  rounding = NULL,
  granularity = c("continuous", "integer", "quasiinteger"),
  scoring = c("singleitem", "sumscored", "meanscored"),
  n_items = 1,
  alpha = NULL
) {
  as_row <- function(bounds) {
    tibble::tibble(
      min_sd = bounds$min_sd,
      max_sd = bounds$max_sd,
      feasible = bounds$feasible,
      min_rule = bounds$min_rule,
      max_rule = bounds$max_rule,
      grim = bounds$grim,
      grimmer = bounds$grimmer,
      sd_in_bounds = bounds$sd_in_bounds,
      note = bounds$note
    )
  }

  # -- validate ----------------------------------------------------------------
  granularity <- rlang::arg_match(granularity)
  scoring <- rlang::arg_match(scoring)

  if (is.character(mean) || is.character(sd)) {
    cli::cli_abort(c(
      "{.arg mean} and {.arg sd} must be numeric.",
      i = "Parse reported strings upstream, e.g. with \\
           {.fn scrutiny::decimal_places} and {.fn as.numeric}."
    ))
  }
  if (
    is.null(n_items) ||
      length(n_items) != 1 ||
      n_items < 1 ||
      !dplyr::near(n_items, round(n_items), tol = TOLERANCE)
  ) {
    cli::cli_abort("{.arg n_items} must be a positive whole number.")
  }
  n_items <- as.integer(round(n_items))

  if (scoring == "singleitem") {
    if (n_items != 1L) {
      cli::cli_abort(
        '{.code scoring = "singleitem"} requires {.code n_items = 1}.'
      )
    }
    if (!is.null(alpha)) {
      cli::cli_abort(c(
        '{.code scoring = "singleitem"} cannot take {.arg alpha}.',
        i = "A single item has no internal consistency.",
        i = 'Use {.code scoring = "sumscored"} or {.code "meanscored"}.'
      ))
    }
  }
  if (!is.null(n) && n < 2) {
    cli::cli_abort("{.arg n} must be >= 2 for a sample SD.")
  }
  if (!is.null(mean) && is.null(n)) {
    cli::cli_abort("Mean-conditional bounds require {.arg n}.")
  }
  if (!is.null(rounding) && is.null(mean)) {
    cli::cli_abort("{.arg rounding} requires a {.arg mean}.")
  }
  if (!is.null(rounding) && is.null(digits_mean)) {
    cli::cli_abort("{.arg rounding} requires {.arg digits_mean}.")
  }
  if (!is.null(sd) && !is.null(rounding) && is.null(digits_sd)) {
    cli::cli_abort(
      "A reported {.arg sd} with {.arg rounding} requires {.arg digits_sd}."
    )
  }
  if (!is.null(l) && !is.null(u) && u <= l) {
    cli::cli_abort("Need {.code u > l}.")
  }
  if (!is.null(a) && !is.null(b) && b < a) {
    cli::cli_abort("Need {.code b >= a}.")
  }
  if (!is.null(a) && !is.null(l) && a < l) {
    cli::cli_abort(
      "Observed minimum {.arg a} cannot lie below the wall {.arg l}."
    )
  }
  if (!is.null(b) && !is.null(u) && b > u) {
    cli::cli_abort(
      "Observed maximum {.arg b} cannot lie above the wall {.arg u}."
    )
  }
  # granularity limits live on the reported scale; the 1/n_items mean-score grid
  # maps to integers under w = n_items * x, so reported limits must be integers.
  if (granularity %in% c("integer", "quasiinteger")) {
    limits <- c(l, u, a, b)
    if (length(limits) && !all(dplyr::near(limits, round(limits)))) {
      cli::cli_abort(
        "Integer and quasi-integer constraints require integer-valued limits."
      )
    }
  }

  # -- granularity multiplier and alpha item count -----------------------------
  # `multiplier` rescales the reported (mean-score) scale to the integer sum
  # scale the granularity and GRIM machinery operate on: w = multiplier * x.
  # Single-item and sum-scored observations are already integer-gridded, so it
  # is 1; mean-scored data sit on a 1/n_items grid. After that rescaling a mean
  # score becomes an n_items-item sum, so the same value serves as both the
  # granularity divisor and the composite item count.
  multiplier <- if (scoring == "meanscored") n_items else 1L

  if (!is.null(alpha)) {
    if (alpha >= 1) {
      cli::cli_abort("{.arg alpha} must be < 1.")
    }
    if (n_items == 1L) {
      alpha <- NULL # alpha is inert at one item
    } else if (is.null(l) || is.null(u) || is.null(n) || is.null(mean)) {
      cli::cli_abort(
        "Alpha bounds require {.arg l}, {.arg u}, {.arg n}, and {.arg mean}."
      )
    } else if (!is.null(a) || !is.null(b)) {
      cli::cli_abort(
        "{.arg alpha} with attained extremes is not supported (open problem)."
      )
    }
  }

  # -- resolve sides: attained supersedes wall ---------------------------------
  lower <- if (!is.null(a)) a else l
  upper <- if (!is.null(b)) b else u
  lower_attained <- !is.null(a)
  upper_attained <- !is.null(b)

  # scaled (w = multiplier * x) copies for the integer-grid core; SD outputs
  # divide by the multiplier again.
  scale_up <- function(x) if (is.null(x)) NULL else x * multiplier
  lower_scaled <- scale_up(lower)
  upper_scaled <- scale_up(upper)

  bounds_at <- function(mean_current) {
    bounds <- sd_bounds_core(
      lower_scaled,
      upper_scaled,
      lower_attained,
      upper_attained,
      n,
      if (is.null(mean_current)) NULL else mean_current * multiplier,
      granularity,
      alpha,
      n_items
    )
    if (!is.na(bounds$min_sd)) {
      bounds$min_sd <- bounds$min_sd / multiplier
    }
    if (!is.na(bounds$max_sd) && is.finite(bounds$max_sd)) {
      bounds$max_sd <- bounds$max_sd / multiplier
    }
    bounds
  }

  # -- scrutiny verdicts (deferred, never reimplemented) -----------------------
  # GRIM and GRIMMER apply when the mean is a rounded report and the data are
  # strictly integer. Verdicts are reported verbatim; any divergence from our
  # own enumeration is surfaced in `note`.
  grim_verdict <- NA
  grimmer_verdict <- NA
  if (granularity == "integer" && !is.null(rounding) && !is.null(mean)) {
    grim_verdict <- isTRUE(as.logical(unname(scrutiny::grim(
      x = mean,
      n = n,
      digits_x = digits_mean,
      items = n_items,
      rounding = rounding
    )))[1])
    if (!is.null(sd)) {
      grimmer_verdict <- isTRUE(as.logical(unname(scrutiny::grimmer(
        x = mean,
        sd = sd,
        n = n,
        digits_x = digits_mean,
        digits_sd = digits_sd,
        items = n_items,
        rounding = rounding
      )))[1])
    }
  }

  # does the reported sd (as an interval if rounded) overlap the bounds?
  sd_overlaps <- function(min_sd, max_sd) {
    if (is.null(sd) || is.na(min_sd)) {
      return(NA)
    }
    if (!is.null(rounding) && !is.null(digits_sd)) {
      interval <- unround_interval(sd, digits_sd, rounding)
      interval$hi >= min_sd - TOLERANCE && interval$lo <= max_sd + TOLERANCE
    } else {
      sd >= min_sd - TOLERANCE && sd <= max_sd + TOLERANCE
    }
  }

  finish <- function(bounds, extra_note = NULL) {
    bounds$grim <- grim_verdict
    bounds$grimmer <- grimmer_verdict
    bounds$sd_in_bounds <- sd_overlaps(bounds$min_sd, bounds$max_sd)
    if (!is.null(extra_note)) {
      bounds$note <- if (is.na(bounds$note)) {
        extra_note
      } else {
        paste(bounds$note, extra_note, sep = "; ")
      }
    }
    as_row(bounds)
  }

  # -- exact-mean path ---------------------------------------------------------
  if (is.null(rounding) || is.null(mean)) {
    return(finish(bounds_at(mean)))
  }

  # -- rounded/truncated mean: envelope over the rounding interval -------------
  interval <- unround_interval(mean, digits_mean, rounding)
  band <- feasible_mean_band(
    lower,
    upper,
    lower_attained,
    upper_attained,
    n
  )
  mean_lo <- max(interval$lo, band[1])
  mean_hi <- min(interval$hi, band[2])

  if (mean_lo > mean_hi + TOLERANCE_TIGHT) {
    return(finish(list(
      min_sd = NA_real_,
      max_sd = NA_real_,
      feasible = FALSE,
      min_rule = NA_character_,
      max_rule = NA_character_,
      note = sprintf(
        "no mean in the rounding interval [%.6g, %.6g] lies in the feasible band [%.6g, %.6g]",
        interval$lo,
        interval$hi,
        band[1],
        band[2]
      )
    )))
  }

  divergence <- NULL
  if (granularity == "integer") {
    # exact: enumerate the GRIM-consistent means in the interval. Under the
    # 1/n_items mean grid a mean is GRIM-consistent when n * multiplier * mean
    # is an integer, so enumerate over that product.
    grid_n <- n * multiplier
    sums <- ceiling(grid_n * mean_lo - TOLERANCE):floor(
      grid_n * mean_hi + TOLERANCE
    )
    candidates <- sums[
      sums >= grid_n * mean_lo - TOLERANCE &
        sums <= grid_n * mean_hi + TOLERANCE
    ] /
      grid_n
    if (!interval$lo_incl) {
      candidates <- candidates[
        !dplyr::near(candidates, interval$lo, tol = TOLERANCE)
      ]
    }
    if (!interval$hi_incl) {
      candidates <- candidates[
        !dplyr::near(candidates, interval$hi, tol = TOLERANCE)
      ]
    }
    # surface any disagreement with scrutiny's deferred GRIM verdict
    if (!is.na(grim_verdict) && grim_verdict != (length(candidates) > 0)) {
      divergence <- paste(
        "scrutiny::grim verdict disagrees with the package's GRIM-mean",
        "enumeration (a known scrutiny floating-point boundary case)"
      )
    }
    if (!length(candidates)) {
      return(finish(
        list(
          min_sd = NA_real_,
          max_sd = NA_real_,
          feasible = FALSE,
          min_rule = NA_character_,
          max_rule = NA_character_,
          note = "no GRIM-consistent mean lies in the rounding interval: the reported mean is not attainable by integer data at this n"
        ),
        extra_note = divergence
      ))
    }
  } else {
    candidates <- candidate_means(mean_lo, mean_hi, n * multiplier)
  }

  results <- candidates |>
    purrr::map(bounds_at) |>
    purrr::keep(\(bounds) bounds$feasible)

  if (!length(results)) {
    return(finish(
      list(
        min_sd = NA_real_,
        max_sd = NA_real_,
        feasible = FALSE,
        min_rule = NA_character_,
        max_rule = NA_character_,
        note = "no mean in the rounding interval yields a feasible constraint set"
      ),
      extra_note = divergence
    ))
  }

  mins <- purrr::map_dbl(results, \(bounds) bounds$min_sd)
  maxs <- purrr::map_dbl(results, \(bounds) bounds$max_sd)
  at_min <- which.min(mins)
  at_max <- which.max(maxs)

  finish(
    list(
      min_sd = mins[at_min],
      max_sd = maxs[at_max],
      feasible = TRUE,
      min_rule = sprintf(
        "%s (envelope over rounding interval)",
        results[[at_min]]$min_rule
      ),
      max_rule = sprintf(
        "%s (envelope over rounding interval)",
        results[[at_max]]$max_rule
      ),
      note = sprintf(
        "mean treated as %s-%s to %d dp: envelope over [%.6g, %.6g]%s",
        if (rounding %in% c("trunc", "anti_trunc", "ceiling", "floor")) {
          "truncated"
        } else {
          "rounded"
        },
        rounding,
        interval$digits,
        interval$lo,
        interval$hi,
        if (mean_lo > interval$lo || mean_hi < interval$hi) {
          sprintf(", clipped to feasible [%.6g, %.6g]", mean_lo, mean_hi)
        } else {
          ""
        }
      )
    ),
    extra_note = divergence
  )
}

# ---- Layer 4b: percent-of-maximum-possible (POMP) transforms -----------------

# Internal: POMP location and two POMP dispersion scores, on the reported scale.
#   pomp_mean       (mean - lower) / (upper - lower), in [0, 1]
#   pomp_sd_parity  s / s_max_parity, the mean-agnostic parity (Popoviciu) max.
#                   A LINEAR rescaling: floor 0, mean-independent, so points
#                   that differ only in s stay ordered and the umbrella geometry
#                   is undistorted.
#   pomp_sd_sharp   (s - min_sd) / (max_sd - min_sd), the position of s within
#                   the SHARP mean-conditional band actually returned for this
#                   constraint set. Non-linear and mean-dependent, but "1" means
#                   attained-at-this-mean and it pools scales onto one
#                   relative-dispersion axis. NA without a mean or on a
#                   degenerate band.
# lower/upper are the EFFECTIVE limits (a supersedes l, b supersedes u).
pomp_cols <- function(mean, sd, min_sd, max_sd, lower, upper, n, has_mean) {
  location <- if (
    !is.null(mean) && !is.null(lower) && !is.null(upper) && (upper - lower) > 0
  ) {
    (mean - lower) / (upper - lower)
  } else {
    NA_real_
  }

  parity_max <- if (!is.null(lower) && !is.null(upper) && !is.null(n)) {
    sd_max_span_n(lower, upper, n)
  } else {
    NA_real_
  }
  dispersion_parity <- if (
    !is.null(sd) && !is.na(parity_max) && parity_max > 0
  ) {
    sd / parity_max
  } else {
    NA_real_
  }

  dispersion_sharp <- if (
    !is.null(sd) &&
      has_mean &&
      !is.na(min_sd) &&
      !is.na(max_sd) &&
      is.finite(max_sd) &&
      (max_sd - min_sd) > TOLERANCE_TIGHT
  ) {
    (sd - min_sd) / (max_sd - min_sd)
  } else {
    NA_real_
  }

  list(
    pomp_mean = location,
    pomp_sd_parity = dispersion_parity,
    pomp_sd_sharp = dispersion_sharp
  )
}
