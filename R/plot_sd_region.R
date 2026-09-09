# One plotting entry point for the feasible-region panels: each `rule` is a
# constraint set from the nested framework, drawn either as a continuous band or
# (for strictly integer data) as the lattice of attainable reported tuples.

# The rule vocabulary, shared by sd_region_data() and plot_sd_region().
REGION_RULES <- c(
  "quasi",
  "range",
  "range_n",
  "mean",
  "muilwijk",
  "mean_naive_floor",
  "mestdagh",
  "pesant_regin",
  "alpha",
  "integer",
  "integer_alpha",
  "attainable",
  "attainable_alpha"
)

ROUNDING_RULES <- c(
  "half_up",
  "half_down",
  "native",
  "ceiling",
  "floor",
  "trunc",
  "anti_trunc"
)

SCORING_KINDS <- c("singleitem", "sumscored", "meanscored")

# Internal: the alpha-free quasi-integer sharp band. `multiplier` sends the
# reported grid to the integers under w = multiplier * (x - l), so it is
# n_items for mean-scored composites and 1 otherwise.
quasi_band <- function(mean, n, l, u, multiplier) {
  tibble::tibble(
    mean = mean,
    lo = sd_min_quasi_integer(multiplier * (mean - l), n) / multiplier,
    hi = sd_max_structure_s(
      multiplier * (mean - l),
      n,
      0,
      multiplier * (u - l)
    ) /
      multiplier
  )
}

# Internal: resolve a `scoring` choice into the quantities the region formulas
# need. `multiplier` is the granularity multiplier (above); `to_sum_scale` maps
# a reported mean to the sum-score mean the alpha bounds are stated in;
# `sd_divisor` converts a sum-score SD back to the reported SD's units.
scoring_geometry <- function(scoring, k, l, u) {
  switch(
    scoring,
    singleitem = {
      if (k != 1L) {
        cli::cli_abort(
          '{.code scoring = "singleitem"} requires {.code n_items = 1}.'
        )
      }
      list(
        multiplier = 1,
        item_l = l,
        item_u = u,
        to_sum_scale = function(mean) mean,
        sd_divisor = 1
      )
    },
    meanscored = list(
      multiplier = k,
      item_l = l,
      item_u = u,
      to_sum_scale = function(mean) k * mean,
      sd_divisor = k
    ),
    sumscored = {
      if (
        !dplyr::near(l / k, round(l / k), tol = TOLERANCE) ||
          !dplyr::near(u / k, round(u / k), tol = TOLERANCE)
      ) {
        cli::cli_abort(c(
          '{.code scoring = "sumscored"} needs {.arg l} and {.arg u} divisible \\
           by {.arg n_items}.',
          i = "They are the composite's limits, {.arg n_items} times the \\
               per-item limits."
        ))
      }
      list(
        multiplier = 1,
        item_l = l / k,
        item_u = u / k,
        to_sum_scale = function(mean) mean,
        sd_divisor = 1
      )
    }
  )
}

gcd2 <- function(a, b) {
  while (b) {
    carry <- b
    b <- a %% b
    a <- carry
  }

  a
}

# Internal: apply one reporting-rounding convention, using roundwork's
# implementations so that the vocabulary matches the rest of the
# error-detection ecosystem. "native" is base R's round(), halves to even.
round_reported <- function(x, digits, rounding) {
  switch(
    rounding,
    half_up = roundwork::round_up(x, digits),
    half_down = roundwork::round_down(x, digits),
    native = round(x, digits),
    ceiling = roundwork::round_ceiling(x, digits),
    floor = roundwork::round_floor(x, digits),
    trunc = roundwork::round_trunc(x, digits),
    anti_trunc = roundwork::round_anti_trunc(x, digits),
    cli::cli_abort("Unknown rounding rule: {.val {rounding}}.")
  )
}

# Internal: round a lattice of (mean, sd) pairs to reporting precision and
# collapse the duplicates that creates. Distinct attainable pairs routinely
# round into one reported cell, which is exactly why a rounded lattice can look
# solid where the exact one is full of holes.
round_lattice <- function(lattice, digits, rounding) {
  if (is.null(digits)) {
    return(lattice)
  }

  lattice |>
    dplyr::mutate(
      mean = round_reported(.data$mean, digits, rounding),
      sd = round_reported(.data$sd, digits, rounding),
      # some rounding rules return negative zero at 0, which compares equal but
      # prints and string-matches differently; normalise it away
      mean = .data$mean + 0,
      sd = .data$sd + 0
    ) |>
    dplyr::distinct() |>
    dplyr::arrange(.data$mean, .data$sd)
}

# Internal: the EXACT attainable (mean, SD) pairs for strictly integer data,
# with no reporting grid. A dynamic program over the n observations accumulates
# the reachable (sum, sum-of-squares) pairs on the shifted grid
# y = multiplier * (x - l) in 0..span, from which both statistics are
# recovered. Unlike the rounded "integer" lattice this shows the true interior
# holes, which rounding smears shut.
#
# Five things keep the state space small enough for this to stay fast in R:
#
#  (1) The second axis is NOT the sum of squares Q but R = span*S - Q =
#      sum y(span - y). Every term is non-negative and at most
#      floor(span^2/4), so R spans n*floor(span^2/4) rather than n*span^2:
#      about a quarter of the cells.
#  (2) R is then divided by the gcd of the achievable y(span - y). For odd
#      span every term is even, so that halves the axis again.
#  (3) Only the frontier reachable after t items is live at step t, so the
#      active window grows instead of the full grid being swept n times.
#  (4) The reachable set is symmetric under y -> span - y, which maps
#      S -> t*span - S and leaves R fixed. Only the lower half of the S axis is
#      computed; the upper half is a reversed copy.
#  (5) Cells are `raw` (one byte) rather than `logical` (four).
#
# Shifted whole-block ORs are used rather than scattering into which(...)
# indices: the block form is what lets (3) and (4) restrict the work.
attainable_lattice <- function(
  l,
  u,
  n,
  multiplier = 1,
  max_cells = MAX_LATTICE_CELLS
) {
  span <- multiplier * (u - l)
  if (!dplyr::near(span, round(span), tol = TOLERANCE)) {
    cli::cli_abort(
      "{.arg l}, {.arg u} and {.arg n_items} must put the scale limits on the \\
       integer grid."
    )
  }
  span <- as.integer(round(span))
  if (span < 1L) {
    cli::cli_abort("The scale limits must span at least one grid step.")
  }

  values <- 0:span
  contributions <- values * (span - values) # each item's contribution to R
  divisor <- Reduce(gcd2, contributions[contributions > 0])
  if (!length(divisor) || is.na(divisor) || divisor < 1) {
    divisor <- 1
  }
  axis_span <- max(contributions) / divisor # per-item span of the R axis
  sum_max <- n * span
  residual_max <- n * axis_span

  if ((sum_max + 1) * (residual_max + 1) > max_cells) {
    cli::cli_abort(c(
      "The exact lattice is too large to enumerate here \\
       ({format(sum_max + 1)} x {format(residual_max + 1)} cells).",
      i = 'Use {.code rule = "integer"} with a reporting precision instead.'
    ))
  }

  zero <- as.raw(0)
  reached <- matrix(zero, sum_max + 1L, residual_max + 1L)
  next_reached <- matrix(zero, sum_max + 1L, residual_max + 1L)
  reached[1L, 1L] <- as.raw(1)
  residual_steps <- contributions / divisor

  for (item_current in seq_len(n)) {
    sum_before <- (item_current - 1L) * span
    residual_before <- (item_current - 1L) * axis_span # frontier before it
    sum_after <- item_current * span
    residual_after <- item_current * axis_span # frontier after it
    half_rows <- floor(sum_after / 2) + 1L # rows computed; rest mirrored
    next_reached[1L:(sum_after + 1L), 1L:(residual_after + 1L)] <- zero

    for (index_current in seq_along(values)) {
      value_current <- values[index_current]
      residual_step <- residual_steps[index_current]
      row_hi <- min(sum_before + 1L + value_current, half_rows)
      if (row_hi < 1L + value_current) {
        next
      }
      next_reached[
        (1L + value_current):row_hi,
        (1L + residual_step):(residual_before + 1L + residual_step)
      ] <- next_reached[
        (1L + value_current):row_hi,
        (1L + residual_step):(residual_before + 1L + residual_step)
      ] |
        reached[
          1L:(row_hi - value_current),
          1L:(residual_before + 1L),
          drop = FALSE
        ]
    }

    if (half_rows < sum_after + 1L) {
      # y -> span - y symmetry
      next_reached[
        (half_rows + 1L):(sum_after + 1L),
        1L:(residual_after + 1L)
      ] <- next_reached[
        (sum_after + 1L - half_rows):1L,
        1L:(residual_after + 1L),
        drop = FALSE
      ]
    }

    swap <- reached
    reached <- next_reached
    next_reached <- swap
  }

  reached_index <- which(reached != zero)
  n_rows <- sum_max + 1L
  sums <- (reached_index - 1L) %% n_rows
  squares <- span * sums - ((reached_index - 1L) %/% n_rows) * divisor
  # the sum of squared deviations is translation-free
  sum_squares <- squares - sums^2 / n

  tibble::tibble(
    mean = l + sums / (n * multiplier),
    sd = sqrt(pmax(0, sum_squares) / (n - 1)) / multiplier
  ) |>
    dplyr::arrange(.data$mean, .data$sd)
}

#' Feasible-region data for one constraint set
#'
#' The `(mean, lo, hi)` band, or the lattice of attainable reported
#' `(mean, sd)` tuples, for a chosen constraint set. See [plot_sd_region()] for
#' the `rule` vocabulary; this is the data behind that plot.
#'
#' @inheritParams plot_sd_region
#' @return A tibble with `mean`, `lo`, `hi` and `sd`, and a `type` attribute of
#'   `"band"` or `"points"`. Band rules fill `lo` and `hi` and leave `sd` `NA`;
#'   the lattice rules (`"integer"`, `"integer_alpha"`, `"attainable"`,
#'   `"attainable_alpha"`) fill `sd` with the tuples passing every applicable
#'   test and leave `lo` and `hi` `NA`.
#' @examples
#' # a band rule fills lo/hi across the scale
#' band <- sd_region_data(l = 1, u = 5, n = 7, rule = "quasi")
#' attr(band, "type")
#'
#' # a lattice rule fills sd with the attainable tuples themselves
#' points <- sd_region_data(l = 1, u = 5, n = 7, rule = "attainable")
#' attr(points, "type")
#' @export
sd_region_data <- function(
  l,
  u,
  n,
  rule = REGION_RULES,
  scoring = NULL,
  n_items = 1,
  alpha = NULL,
  digits = 2,
  by = NULL,
  round_digits = NULL,
  rounding = ROUNDING_RULES
) {
  rule <- rlang::arg_match(rule)
  if (rule == "muilwijk") {
    rule <- "mean" # alias: named for its author
  }
  rounding <- rlang::arg_match(rounding)

  if (
    rule %in% c("alpha", "integer_alpha", "attainable_alpha") && is.null(alpha)
  ) {
    cli::cli_abort("{.code rule = {rule}} requires {.arg alpha}.")
  }

  # back-compatible default: a composite is in mean-score units unless told
  # otherwise, and one item is a single item.
  if (is.null(scoring)) {
    scoring <- if (n_items > 1) "meanscored" else "singleitem"
  }
  scoring <- rlang::arg_match(scoring, SCORING_KINDS)
  geometry <- scoring_geometry(scoring, n_items, l, u)

  as_points <- function(points) {
    out <- points |>
      round_lattice(round_digits, rounding) |>
      dplyr::transmute(
        mean = .data$mean,
        lo = NA_real_,
        hi = NA_real_,
        sd = .data$sd
      )
    attr(out, "type") <- "points"
    out
  }

  # lattice rules: only GRIM/GRIMMER-attainable reported tuples exist, so the
  # region is a set of points, not a band.
  if (rule %in% c("integer", "integer_alpha")) {
    points <- umbrella_data(
      n = n,
      l = l,
      u = u,
      digits = digits,
      granularity = "integer",
      scoring = scoring,
      n_items = n_items,
      alpha = if (rule == "integer_alpha") alpha else NULL
    ) |>
      dplyr::filter(.data$consistent) |>
      dplyr::select("mean", "sd")

    return(as_points(points))
  }

  # exact lattice rules: the true attainable tuples, with no reporting grid.
  if (rule %in% c("attainable", "attainable_alpha")) {
    points <- attainable_lattice(l, u, n, geometry$multiplier)

    if (rule == "attainable_alpha") {
      bounds <- unique(points$mean) |>
        purrr::map(function(mean_current) {
          at_mean <- sd_bounds(
            l = l,
            u = u,
            n = n,
            mean = mean_current,
            granularity = "integer",
            scoring = scoring,
            n_items = n_items,
            alpha = alpha
          )
          tibble::tibble(
            mean = mean_current,
            band_lo = at_mean$min_sd,
            band_hi = at_mean$max_sd,
            ok = isTRUE(at_mean$feasible) && !is.na(at_mean$min_sd)
          )
        }) |>
        purrr::list_rbind()

      points <- points |>
        dplyr::inner_join(bounds, by = "mean") |>
        dplyr::filter(
          .data$ok,
          .data$sd >= .data$band_lo - TOLERANCE,
          .data$sd <= .data$band_hi + TOLERANCE
        ) |>
        dplyr::select("mean", "sd") |>
        dplyr::arrange(.data$mean, .data$sd)
    }

    return(as_points(points))
  }

  if (is.null(by)) {
    by <- (u - l) / MEAN_GRID_DIVISIONS
  }
  means <- seq(l, u, by = by)
  quasi <- quasi_band(means, n, l, u, geometry$multiplier)
  # strict Bernoulli floor, defined off-grid too
  naive_floor <- sd_min_integer(geometry$multiplier * means, n) /
    geometry$multiplier

  out <- switch(
    rule,
    range = tibble::tibble(mean = means, lo = 0, hi = sd_max_span(l, u)),
    range_n = tibble::tibble(
      mean = means,
      lo = 0,
      hi = sd_max_span_n(l, u, n)
    ),
    mean = tibble::tibble(
      mean = means,
      lo = 0,
      hi = sd_max_muilwijk(means, n, l, u)
    ),
    mean_naive_floor = tibble::tibble(
      mean = means,
      lo = naive_floor,
      hi = sd_max_muilwijk(means, n, l, u)
    ),
    mestdagh = tibble::tibble(mean = means, lo = 0, hi = quasi$hi),
    pesant_regin = tibble::tibble(
      mean = means,
      lo = quasi$lo,
      hi = sd_max_span_n(l, u, n)
    ),
    quasi = tibble::tibble(mean = means, lo = quasi$lo, hi = quasi$hi),
    alpha = {
      if (n_items < 2) {
        cli::cli_abort(
          '{.code rule = "alpha"} needs {.code n_items >= 2} (alpha is inert \\
           at one item).'
        )
      }
      design_divisor <- 1 - ((n_items - 1) / n_items) * alpha
      if (design_divisor <= TOLERANCE_TIGHT) {
        cli::cli_abort("{.arg alpha} is too high for this {.arg n_items}.")
      }
      ceiling_alpha <- sqrt(
        (n / (n - 1)) *
          v_max_alpha(
            geometry$to_sum_scale(means),
            n_items,
            n,
            geometry$item_l,
            geometry$item_u
          ) /
          design_divisor
      ) /
        geometry$sd_divisor
      tibble::tibble(
        mean = means,
        # alpha-amplified quasi-integer floor; alpha can only tighten
        lo = quasi$lo / sqrt(design_divisor),
        hi = pmin(ceiling_alpha, quasi$hi)
      )
    }
  )

  # an empty band (floor above ceiling) is infeasible, not a negative region
  out <- out |>
    dplyr::mutate(
      infeasible = .data$lo > .data$hi + TOLERANCE_TIGHT,
      lo = dplyr::if_else(.data$infeasible, NA_real_, .data$lo),
      hi = dplyr::if_else(.data$infeasible, NA_real_, .data$hi),
      sd = NA_real_
    ) |>
    dplyr::select("mean", "lo", "hi", "sd")

  attr(out, "type") <- "band"

  out
}

#' Plot the feasible SD region for a chosen constraint set
#'
#' Draws the region of sample SDs a constraint set asserts to be possible, as a
#' function of the mean, with the sharp quasi-integer band optionally repeated
#' as a dashed reference. This reproduces the panels of the nested-constraints
#' figure in the STRAIT article from a single entry point.
#'
#' Rules, in the order they enter the framework (`n_items > 1` puts every rule
#' in mean-score units for a composite of that many integer items):
#'
#' \describe{
#'   \item{`"range"`}{`(u - l)/sqrt(2)`, no floor; ignores `n` and the mean
#'     (Popoviciu 1935, range only).}
#'   \item{`"range_n"`}{the parity ceiling, no floor (Popoviciu, as restored by
#'     Petocz 2005).}
#'   \item{`"mean"`, `"muilwijk"`}{Muilwijk's mean-conditional ceiling as
#'     originally stated, over a floor of zero (Muilwijk 1966; Bhatia-Davis
#'     2000). A valid ceiling at every mean, but attained only where the counts
#'     at each limit come out whole, so it is conservative elsewhere; see
#'     `"mestdagh"` and [sd_delta()]. The two names are aliases.}
#'   \item{`"mean_naive_floor"`}{that arch plus the naive Bernoulli floor
#'     (Fuenderich et al. 2025).}
#'   \item{`"pesant_regin"`}{integer minimum with a loose (parity) ceiling
#'     (Pesant and Regin 2005).}
#'   \item{`"mestdagh"`}{the sharp integer maximum, no floor (Mestdagh et al.
#'     2018): `"mean"` with the count-parity correction, hence never above it.}
#'   \item{`"quasi"`}{both bounds sharp and GRIM-free (this package's default).}
#'   \item{`"alpha"`}{additionally conditioning on a reported Cronbach's
#'     `alpha`; needs `n_items >= 2`.}
#'   \item{`"integer"`}{strictly integer data: the lattice of reported
#'     `(mean, sd)` tuples passing GRIM, GRIMMER and the bounds — exactly what
#'     [brimmer()] admits. Built via [umbrella_data()], so the pairs are
#'     rounded to `digits`, and drawn as points.}
#'   \item{`"integer_alpha"`}{that lattice, additionally inside the
#'     alpha-conditional bounds.}
#'   \item{`"attainable"`}{the EXACT attainable `(mean, sd)` tuples of strictly
#'     integer data, with no reporting grid. Unlike `"integer"` this shows the
#'     true interior holes, which rounding smears shut; it is enumerated by
#'     dynamic programming and errors if the lattice is too large.}
#'   \item{`"attainable_alpha"`}{those exact tuples, additionally inside the
#'     alpha-conditional bounds.}
#' }
#'
#' The two lattice rules are easy to confuse and the difference is not small.
#' GRIMMER is necessary but not sufficient for an integer sample to exist, so
#' the attainable set is a strict subset of the GRIMMER-consistent one: at
#' `l = 1, u = 5, n = 10, digits = 2` there are 491 GRIMMER-consistent reported
#' pairs against 447 attainable ones. Use [brimmest()] to certify a single
#' tuple rather than rounding a lattice to compare against.
#'
#' @param l,u Numeric scalars, the scale limits (mean-score units when
#'   `n_items > 1` and `scoring` is left at its default).
#' @param n Integer scalar, sample size.
#' @param rule Which constraint set to draw; see Details.
#' @param scoring One of `"singleitem"`, `"sumscored"`, `"meanscored"`, as in
#'   [sd_bounds()]. Defaults to `"meanscored"` when `n_items > 1` and
#'   `"singleitem"` otherwise.
#' @param n_items Integer, number of items in the composite (default 1).
#' @param alpha Reported Cronbach's alpha; required by the alpha rules.
#' @param digits Reported decimal places, used by the rounded lattice rules
#'   (default 2); ignored by the exact `"attainable"` rules.
#' @param round_digits Optionally round the returned `(mean, sd)` lattice to
#'   this many decimal places, collapsing pairs that round together. `NULL`
#'   (default) returns them unrounded. Band rules are unaffected: rounding a
#'   bound can turn it anti-conservative.
#' @param rounding How to round when `round_digits` is given. `"half_up"`
#'   (default) and the other named rules use roundwork's implementations;
#'   `"native"` is base R's `round()`, which rounds halves to even.
#' @param show_reference Draw the alpha-free sharp quasi-integer band as a
#'   dashed reference (default `TRUE`; skipped for `rule = "quasi"`, which is
#'   that band).
#' @param title Optional plot title.
#' @param by Optional mean-grid spacing (default `(u - l) / 1000`); ignored by
#'   the lattice rules, which step by `10^-digits`.
#' @param shade `"outside"` (default) shades the infeasible region and leaves
#'   the feasible one clear, matching [plot_sd_bounds()]; `"inside"` fills the
#'   band. The outside form is drawn with [band_polygon()], so a rule whose
#'   band is undefined at some means leaves those means shaded rather than
#'   blank.
#' @param expand Padding around the plotted region, as a proportion of the
#'   scale width `u - l` rather than a fixed number of SD units. Applied to
#'   both axes.
#' @param fill,line_colour,reference_colour,point_colour,point_size Appearance.
#' @return A ggplot object.
#' @examples
#' plot_sd_region(l = 1, u = 5, n = 7, rule = "mean") # Bhatia-Davis
#' plot_sd_region(l = 1, u = 5, n = 7, rule = "quasi") # sharp, GRIM-free
#' \donttest{
#' # digits = 1 keeps the reporting grid small: the lattice rules screen every
#' # cell with GRIMMER
#' plot_sd_region(l = 1, u = 5, n = 7, rule = "integer", digits = 1)
#' plot_sd_region(l = 1, u = 5, n = 7, rule = "alpha", n_items = 2, alpha = 0.7)
#' }
#' @export
plot_sd_region <- function(
  l,
  u,
  n,
  rule = REGION_RULES,
  scoring = NULL,
  n_items = 1,
  alpha = NULL,
  digits = 2,
  round_digits = NULL,
  rounding = "half_up",
  show_reference = TRUE,
  title = NULL,
  by = NULL,
  shade = c("outside", "inside"),
  expand = 0.03,
  fill = "grey85",
  line_colour = "black",
  reference_colour = "grey45",
  point_colour = COLOUR_LATTICE,
  point_size = 0.5
) {
  rule <- rlang::arg_match(rule)
  if (rule == "muilwijk") {
    rule <- "mean" # alias: named for its author
  }
  shade <- rlang::arg_match(shade)

  region <- sd_region_data(
    l = l,
    u = u,
    n = n,
    rule = rule,
    scoring = scoring,
    n_items = n_items,
    alpha = alpha,
    digits = digits,
    by = by,
    round_digits = round_digits,
    rounding = rounding
  )
  geometry <- scoring_geometry(
    if (is.null(scoring)) {
      if (n_items > 1) "meanscored" else "singleitem"
    } else {
      rlang::arg_match(scoring, SCORING_KINDS)
    },
    n_items,
    l,
    u
  )
  is_points <- identical(attr(region, "type"), "points")
  step <- if (is.null(by)) (u - l) / MEAN_GRID_DIVISIONS else by

  plot <- ggplot2::ggplot()

  # the region itself first, so the dashed reference stays visible on top of it
  if (!is_points) {
    plot <- plot +
      if (shade == "inside") {
        ggplot2::geom_ribbon(
          data = region,
          ggplot2::aes(
            x = .data$mean,
            ymin = .data$lo,
            ymax = .data$hi
          ),
          fill = fill,
          na.rm = TRUE
        )
      } else {
        # The alpha rules leave stretches near each limit where no composite
        # exists and `lo`/`hi` are NA; a ribbon draws nothing there, so
        # assembled shading would leave those means unshaded and imply every SD
        # is possible at them. Knocking rings out of a shaded panel keeps them
        # shaded. See band_polygon().
        knockout_layers(band_polygon(region, by = step))
      }
  }

  # dashed reference: the sharp alpha-free band this rule is judged against
  reference <- if (isTRUE(show_reference) && rule != "quasi") {
    quasi_band(seq(l, u, by = step), n, l, u, geometry$multiplier)
  }
  if (!is.null(reference)) {
    plot <- plot +
      ggplot2::geom_line(
        data = reference,
        ggplot2::aes(.data$mean, .data$hi),
        colour = reference_colour,
        linetype = "dashed",
        linewidth = 0.25,
        na.rm = TRUE
      ) +
      ggplot2::geom_line(
        data = reference,
        ggplot2::aes(.data$mean, .data$lo),
        colour = reference_colour,
        linetype = "dashed",
        linewidth = 0.25,
        na.rm = TRUE
      )
  }

  plot <- plot +
    if (is_points) {
      ggplot2::geom_point(
        data = region,
        ggplot2::aes(.data$mean, .data$sd),
        colour = point_colour,
        size = point_size,
        na.rm = TRUE
      )
    } else {
      list(
        ggplot2::geom_line(
          data = region,
          ggplot2::aes(.data$mean, .data$hi),
          colour = line_colour,
          linewidth = 0.4,
          na.rm = TRUE
        ),
        ggplot2::geom_line(
          data = region,
          ggplot2::aes(.data$mean, .data$lo),
          colour = line_colour,
          linewidth = 0.4,
          na.rm = TRUE
        )
      )
    }

  y_hi <- if (is_points) {
    max(region$sd, na.rm = TRUE)
  } else {
    max(c(region$hi, reference$hi), na.rm = TRUE)
  }

  plot +
    padded_coord(l, u, y_hi, expand) +
    ggplot2::labs(x = "Mean", y = "Sample standard deviation", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
