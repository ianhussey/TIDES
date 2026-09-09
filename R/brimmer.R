# Report-level consistency checking (single row and batch) and the data
# builders (bound curves, umbrella grid) that drive the plots.

# ---- Layer 5: report-level consistency check ---------------------------------

#' BRIMMER: check a reported (mean, SD, n) against the SD bounds
#'
#' Bounds-Related Inconsistency of Means and Errors Reported. The SD-side
#' bounds test, and the package's main entry point: it asks whether a reported
#' standard deviation is arithmetically attainable given the scale limits, the
#' sample size and the reported mean. [brim()] is the mean-side test, which
#' needs no SD.
#'
#' BRIMMER is nested on top of BRIM exactly as GRIMMER is nested on top of
#' GRIM. A report is consistent only if every applicable test passes, and each
#' names a distinct defect:
#'
#' * `in_scale_range` — the BRIM predicate: can the reported mean's rounding
#'   interval meet the feasible mean band at all? See [feasible_mean_band()].
#' * `bounds` — does the reported SD's own rounding interval overlap
#'   `[min_sd, max_sd]`?
#' * `feasibility` — a residual: the constraint set admits no sample for a
#'   reason none of the other tests already accounts for, such as an alpha
#'   floor exceeding the alpha ceiling.
#' * `grim`, `grimmer` — scrutiny's verdicts, deferred verbatim (only under
#'   `granularity = "integer"` with `rounding`). When GRIMMER is the only
#'   failing test the note says so, since that pattern matches scrutiny's
#'   documented false-flag family.
#'
#' Unlike [sd_bounds()], `rounding` here defaults to `"up_or_down"`, since a
#' report-checking context virtually always deals with rounded values. Pass
#' `rounding = NULL` to treat the inputs as exact.
#'
#' @inheritParams sd_bounds
#' @param n Integer scalar, sample size.
#' @param sd Numeric scalar or NULL, the reported SD — the thing being checked.
#'   NULL runs the mean-side tests only, which is what [brim()] does.
#' @param rounding Rounding rule for unrounding the reported values (default
#'   `"up_or_down"`); NULL treats inputs as exact. See [unround_interval()].
#'
#' @return A one-row tibble: `consistent` (logical), `failed_tests`
#'   (comma-separated names, `""` if none), the [sd_bounds()] columns `min_sd`,
#'   `max_sd`, `feasible`, the BRIM predicate `in_scale_range`, then `grim`,
#'   `grimmer`, `sd_in_bounds`, the POMP transforms `pomp_mean`,
#'   `pomp_sd_parity` and `pomp_sd_sharp`, and `note`.
#'
#' @seealso [brim()] for the mean-side test, [brimmest()] for the exact
#'   certificate, and [brimmer_map()] to apply this across a data frame.
#' @examples
#' # a perfectly ordinary report on a 1-7 scale
#' brimmer(
#'   l = 1, u = 7, n = 30, mean = 2.97, digits_mean = 2,
#'   sd = 2.83, digits_sd = 2, granularity = "integer"
#' )
#'
#' # SD below the quasi-integer floor: inconsistent via the bounds test
#' brimmer(
#'   l = 1, u = 7, n = 30, mean = 2.97, digits_mean = 2,
#'   sd = 0.10, digits_sd = 2, granularity = "quasiinteger"
#' )
#'
#' # GRIM-impossible mean: inconsistent via grim
#' brimmer(
#'   l = 1, u = 7, n = 30, mean = 3.51, digits_mean = 2,
#'   sd = 1.00, digits_sd = 2, granularity = "integer"
#' )
#' @export
brimmer <- function(
  l = NULL,
  u = NULL,
  a = NULL,
  b = NULL,
  n = NULL,
  mean = NULL,
  digits_mean = NULL,
  sd = NULL,
  digits_sd = NULL,
  rounding = "up_or_down",
  granularity = c("continuous", "integer", "quasiinteger"),
  scoring = c("singleitem", "sumscored", "meanscored"),
  n_items = 1,
  alpha = NULL
) {
  if (is.null(sd) && is.null(mean)) {
    cli::cli_abort(
      "{.fn brimmer} needs a reported {.arg sd}, a reported {.arg mean}, or both."
    )
  }
  if (is.null(n)) {
    cli::cli_abort("{.fn brimmer} requires {.arg n}.")
  }

  # a mean-free check is legitimate (e.g. sd vs the parity ceiling), but
  # rounding then applies to the sd only
  bounds <- sd_bounds(
    l = l,
    u = u,
    a = a,
    b = b,
    n = n,
    mean = mean,
    digits_mean = digits_mean,
    sd = sd,
    digits_sd = digits_sd,
    rounding = if (is.null(mean)) NULL else rounding,
    granularity = granularity,
    scoring = scoring,
    n_items = n_items,
    alpha = alpha
  )

  # mean-free path: sd_bounds() skipped unrounding; redo the sd overlap with
  # the sd's own interval if rounding is in force
  if (is.null(mean) && !is.null(rounding)) {
    if (is.null(digits_sd)) {
      cli::cli_abort(
        "{.arg rounding} requires {.arg digits_sd} for the reported {.arg sd}."
      )
    }
    interval <- unround_interval(sd, digits_sd, rounding)
    bounds$sd_in_bounds <- if (is.na(bounds$min_sd)) {
      NA
    } else {
      interval$hi >= bounds$min_sd - TOLERANCE &&
        interval$lo <= bounds$max_sd + TOLERANCE
    }
  }

  # POMP transforms on the reported scale (effective limits: a supersedes l,
  # b supersedes u). Descriptive, so they use the reported point mean/sd, not
  # the unrounding interval.
  effective_lower <- if (!is.null(a)) a else l
  effective_upper <- if (!is.null(b)) b else u

  # The BRIM predicate, computed independently of sd_bounds(). sd_bounds()
  # legitimately folds GRIM into `feasible` under strict granularity, but an
  # out-of-range mean and a granular-impossible mean are different defects and
  # must not be reported as one.
  in_scale_range <- NA
  if (
    !is.null(mean) &&
      (!is.null(l) || !is.null(u) || !is.null(a) || !is.null(b))
  ) {
    band <- feasible_mean_band(
      lower = effective_lower,
      upper = effective_upper,
      lower_attained = !is.null(a),
      upper_attained = !is.null(b),
      n = n
    )
    interval <- if (is.null(rounding)) {
      list(lo = mean, hi = mean)
    } else {
      unround_interval(mean, digits_mean, rounding)
    }
    in_scale_range <- interval$hi >= band[1] - TOLERANCE &&
      interval$lo <= band[2] + TOLERANCE
  }

  failed <- character(0)
  if (isFALSE(in_scale_range)) {
    failed <- c(failed, "in_scale_range")
  }
  if (isFALSE(bounds$sd_in_bounds)) {
    failed <- c(failed, "bounds")
  }
  if (isFALSE(bounds$grim)) {
    failed <- c(failed, "grim")
  }
  if (isFALSE(bounds$grimmer)) {
    failed <- c(failed, "grimmer")
  }
  # `feasibility` is the residual: no sample exists for a reason none of the
  # named tests above already accounts for.
  if (!isTRUE(bounds$feasible) && !length(failed)) {
    failed <- c(failed, "feasibility")
  }

  note <- bounds$note
  if (identical(failed, "grimmer")) {
    caveat <- paste(
      "only GRIMMER failed while feasibility and bounds pass; scrutiny's",
      "GRIMMER has documented false-flag cases - verify before flagging"
    )
    note <- if (is.na(note)) caveat else paste(note, caveat, sep = "; ")
  }

  pomp <- pomp_cols(
    mean,
    sd,
    bounds$min_sd,
    bounds$max_sd,
    effective_lower,
    effective_upper,
    n,
    has_mean = !is.null(mean)
  )

  dplyr::bind_cols(
    tibble::tibble(
      consistent = length(failed) == 0,
      failed_tests = paste(failed, collapse = ",")
    ),
    bounds[, c("min_sd", "max_sd", "feasible")],
    tibble::tibble(in_scale_range = in_scale_range),
    bounds[, c("grim", "grimmer", "sd_in_bounds")],
    tibble::tibble(
      pomp_mean = pomp$pomp_mean,
      pomp_sd_parity = pomp$pomp_sd_parity,
      pomp_sd_sharp = pomp$pomp_sd_sharp,
      note = note
    )
  )
}

#' BRIM: check a reported mean against the scale bounds
#'
#' Bounds-Related Inconsistency of Means. The mean-side test: is the reported
#' mean attainable at all, given the scale limits, the sample size and, when
#' supplied, the attained extremes? Where GRIM asks whether a mean is
#' attainable by strictly integer data at a given `n`, BRIM asks whether it is
#' attainable within the reporting range at all — a weaker condition that
#' applies to continuous data too, and one that tightens sharply once an
#' observed minimum or maximum is reported.
#'
#' A thin wrapper on [brimmer()] with no reported SD, so only `in_scale_range`
#' and (under `granularity = "integer"` with `rounding`) scrutiny's `grim`
#' verdict apply. The two name different defects: a mean can be out of range
#' but granular-attainable (7.50 on a 1-7 scale at `n = 30`), inside the range
#' but GRIM-impossible (3.51 at `n = 30`), or both (7.51).
#'
#' @inheritParams brimmer
#' @return A one-row tibble: `consistent`, `failed_tests`, `in_scale_range`,
#'   `grim`, the feasible mean band `band_lo` and `band_hi`, `pomp_mean`, and
#'   `note`.
#' @seealso [brimmer()] for the SD-side test, [brimmest()] for the exact
#'   certificate.
#' @examples
#' brim(l = 1, u = 7, n = 30, mean = 2.97, digits_mean = 2)
#'
#' # a mean above the scale maximum cannot be attained
#' brim(l = 1, u = 7, n = 30, mean = 7.5, digits_mean = 1)
#'
#' # attained extremes narrow the band, and can exclude a mean the bare
#' # scale limits would allow
#' brim(a = 1, b = 7, n = 30, mean = 1.10, digits_mean = 2)
#' @export
brim <- function(
  l = NULL,
  u = NULL,
  a = NULL,
  b = NULL,
  n = NULL,
  mean = NULL,
  digits_mean = NULL,
  rounding = "up_or_down",
  granularity = c("continuous", "integer", "quasiinteger"),
  scoring = c("singleitem", "sumscored", "meanscored"),
  n_items = 1
) {
  if (is.null(mean)) {
    cli::cli_abort(
      "{.fn brim} checks a reported mean: {.arg mean} is required."
    )
  }

  # sd = NULL, not NA: NULL is this package's "not supplied", whereas an NA
  # would propagate through the bounds arithmetic instead of switching the
  # SD-side tests off
  out <- brimmer(
    l = l,
    u = u,
    a = a,
    b = b,
    n = n,
    mean = mean,
    digits_mean = digits_mean,
    sd = NULL,
    digits_sd = NULL,
    rounding = rounding,
    granularity = granularity,
    scoring = scoring,
    n_items = n_items
  )

  # same resolution of sides and attainment that sd_bounds() applies, so the
  # reported band is the one feasibility was actually tested against
  band <- feasible_mean_band(
    lower = if (!is.null(a)) a else l,
    upper = if (!is.null(b)) b else u,
    lower_attained = !is.null(a),
    upper_attained = !is.null(b),
    n = n
  )

  dplyr::bind_cols(
    out[, c("consistent", "failed_tests", "in_scale_range", "grim")],
    tibble::tibble(band_lo = band[1], band_hi = band[2]),
    out[, c("pomp_mean", "note")]
  )
}

# ---- Layer 6: batch report checking ------------------------------------------

# Internal: resolve the recognised arguments of a batch call to per-row
# vectors. Shared by brimmer_map() and brimmest_map(), which take the same
# shape of input: a data frame whose columns supply per-row values, plus `...`
# constants broadcast to every row.
#
# `row_args` are the names that may come either way; `call_args` are names that
# apply to the whole call (a set of rounding rules, a cost budget) and so may
# only be constants -- passing one as a column is an error rather than a silent
# per-row reading. Returns only the arguments actually supplied, so the caller
# reads the names off the result rather than tracking presence itself.
resolve_row_args <- function(data, consts, row_args, call_args = character(0)) {
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.")
  }

  unknown <- setdiff(names(consts), c(row_args, call_args))
  if (length(unknown)) {
    cli::cli_abort(
      "{cli::qty(length(unknown))}Unknown constant argument{?s}: {.arg {unknown}}."
    )
  }

  clash <- intersect(call_args, names(data))
  if (length(clash)) {
    cli::cli_abort(c(
      "{cli::qty(length(clash))}{.arg {clash}} appl{?ies/y} to the whole \\
       call, not to one row.",
      i = "{cli::qty(length(clash))}Supply {?it/them} as a constant rather \\
           than a column."
    ))
  }

  cols <- row_args |>
    purrr::map(function(name_current) {
      in_column <- name_current %in% names(data)
      in_constant <- name_current %in%
        names(consts) &&
        !is.null(consts[[name_current]])
      if (in_column && in_constant) {
        cli::cli_abort(
          "{.arg {name_current}} is supplied as both a column and a constant."
        )
      }
      if (in_column) {
        data[[name_current]]
      } else if (in_constant) {
        rep(consts[[name_current]], length.out = nrow(data))
      } else {
        NULL
      }
    }) |>
    rlang::set_names(row_args)

  purrr::compact(cols)
}

# Internal: a grouping key over parallel vectors, for computing once per
# distinct input tuple instead of once per row. format() at a fixed nsmall
# keeps 3 and 3.0 in the same group without the ambiguity of pasting raw
# numerics, and "\r" cannot occur in a formatted number, so the join is
# unambiguous.
row_key <- function(values) {
  do.call(
    paste,
    c(purrr::map(values, format, nsmall = 6, trim = TRUE), sep = "\r")
  )
}

#' Check many reported (mean, SD, n) rows against the SD bounds
#'
#' Applies [brimmer()] to each row of a data frame, de-duplicating identical
#' constraint tuples so the rounding envelope is computed once per distinct
#' input and reused. Columns of `data` whose names match [brimmer()] arguments
#' are taken per row; any argument given in `...` is a constant broadcast to
#' every row. Supplying one name both ways is an error.
#'
#' Recognised names (column or constant): `l`, `u`, `a`, `b`, `n`, `mean`,
#' `digits_mean`, `sd`, `digits_sd`, `rounding`, `granularity`, `scoring`,
#' `n_items`, `alpha`.
#'
#' @param data A data frame, one reported statistic set per row.
#' @param ... Constant arguments applied to all rows.
#' @param include_inputs If TRUE (default), returns `data` column-bound to the
#'   results; if FALSE, only the result columns, in the same row order.
#' @return A tibble of the [brimmer()] columns, one row per input row
#'   (optionally with the inputs prepended).
#' @examples
#' reports <- tibble::tibble(
#'   mean = c(2.97, 3.51, 4.20),
#'   sd = c(2.83, 3.50, 0.90),
#'   n = c(30, 30, 30)
#' )
#'
#' # row 2 reports an SD above the ceiling for a 1-7 scale at n = 30
#' reports |>
#'   brimmer_map(l = 1, u = 7, digits_mean = 2, digits_sd = 2) |>
#'   dplyr::select("mean", "sd", "n", "consistent", "failed_tests")
#' @export
brimmer_map <- function(data, ..., include_inputs = TRUE) {
  arg_names <- c(
    "l",
    "u",
    "a",
    "b",
    "n",
    "mean",
    "digits_mean",
    "sd",
    "digits_sd",
    "rounding",
    "granularity",
    "scoring",
    "n_items",
    "alpha"
  )
  cols <- resolve_row_args(data, list(...), arg_names)
  present <- names(cols)
  if (!("sd" %in% present)) {
    cli::cli_abort(
      "A reported {.arg sd} is required, as a column of {.arg data} or a constant."
    )
  }

  # de-duplicate identical input tuples; compute once per unique tuple
  key <- row_key(cols)
  unique_rows <- which(!duplicated(key))
  back <- match(key, key[unique_rows])

  out <- unique_rows |>
    purrr::map(function(row_current) {
      args <- purrr::map(cols, \(column) column[row_current])
      rlang::exec(brimmer, !!!args)
    }) |>
    purrr::list_rbind()

  out <- out[back, , drop = FALSE]

  if (include_inputs) dplyr::bind_cols(data, out) else out
}

# ---- Layer 7: bound curves and the umbrella grid ------------------------------

#' SD-bound curves across the mean (hole-free under granularity)
#'
#' Traces the minimum and maximum sample SD as the mean sweeps `[l, u]`, on a
#' grid dense enough to resolve the piecewise floors and the Structure S
#' ceiling, whose kinks fall on the `1/(n * n_items)` lattice. Under a
#' granularity constraint the quasi-integer floor is defined at every mean, so
#' the curves have no gaps.
#'
#' @param l,u Numeric scalars, limits.
#' @param n Integer scalar, sample size.
#' @param granularity Granularity; `"quasiinteger"` (default) gives a hole-free
#'   floor.
#' @param scoring,n_items,alpha As in [sd_bounds()].
#' @param by Numeric or NULL; mean-grid spacing (default `(u - l) / 1000`).
#' @return A tibble: `mean`, `min_sd`, `max_sd`, `feasible`, `pomp_mean`,
#'   `parity_max`, `ceil_parity` (= `max_sd / parity_max`), `floor_parity`. The
#'   `"step"` attribute records the uniform spacing `by` actually used. The
#'   mean grid is deliberately NOT uniform — the kinks and their neighbourhoods
#'   are sampled far more densely than `by` — so the spacing cannot be
#'   recovered from the returned means, and a consumer that must tell a
#'   sampling gap from a genuine gap in the band (as [band_polygon()], and
#'   hence [plot_sd_bounds()], must) has to be told.
#' @examples
#' curve <- sd_bounds_curve(l = 1, u = 7, n = 30, by = 0.1)
#'
#' # the ceiling peaks near the scale midpoint
#' curve[which.max(curve$max_sd), c("mean", "min_sd", "max_sd")]
#' @export
sd_bounds_curve <- function(
  l,
  u,
  n,
  granularity = "quasiinteger",
  scoring = "singleitem",
  n_items = 1,
  alpha = NULL,
  by = NULL
) {
  multiplier <- if (scoring == "meanscored") n_items else 1L
  if (is.null(by)) {
    by <- (u - l) / MEAN_GRID_DIVISIONS
  }

  grid_n <- n * multiplier
  kinks <- unique(c(
    seq(ceiling(l * grid_n), floor(u * grid_n)) / grid_n,
    seq(ceiling(l * multiplier), floor(u * multiplier)) / multiplier
  ))
  means <- sort(unique(pmin(
    u,
    pmax(
      l,
      c(
        seq(l, u, by = by),
        kinks,
        kinks + TOLERANCE,
        kinks - TOLERANCE
      )
    )
  )))

  out <- means |>
    purrr::map(function(mean_current) {
      bounds <- sd_bounds(
        l = l,
        u = u,
        n = n,
        mean = mean_current,
        granularity = granularity,
        scoring = scoring,
        n_items = n_items,
        alpha = alpha
      )
      tibble::tibble(
        mean = mean_current,
        min_sd = bounds$min_sd,
        max_sd = bounds$max_sd,
        feasible = bounds$feasible
      )
    }) |>
    purrr::list_rbind()

  parity_max <- sd_max_span_n(l, u, n)
  out <- out |>
    dplyr::mutate(
      pomp_mean = (.data$mean - l) / (u - l),
      parity_max = parity_max,
      ceil_parity = .data$max_sd / parity_max,
      floor_parity = .data$min_sd / parity_max
    )

  attr(out, "step") <- by

  out
}

#' Build the GRIM x GRIMMER x bounds umbrella grid (full, with verdicts)
#'
#' For each reported mean on the `10^-digits` grid over `[l, u]`, computes the
#' SD-bounds envelope once, then evaluates every reported SD on the same grid
#' from 0 up to that mean's ceiling, tagging each `(mean, sd)` cell with the
#' bounds-overlap and (under `granularity = "integer"`) GRIMMER verdicts. This
#' returns the FULL grid — feasible means x candidate SDs — so the plotting
#' layer can render failures as well as the passing "umbrella".
#'
#' @param n Integer scalar, sample size.
#' @param l,u Numeric scalars, limits.
#' @param digits Integer, reported decimal places (grid step `10^-digits`).
#' @param granularity Granularity (default `"integer"`).
#' @param scoring,n_items As in [sd_bounds()].
#' @param alpha Optional reported Cronbach's alpha; when supplied the bounds
#'   each cell is tested against are the alpha-conditional ones.
#' @param rounding Rounding rule for mean and SD (default `"up_or_down"`).
#' @return A tibble: `mean`, `sd`, `min_sd`, `max_sd`, `in_bounds`, `grimmer`,
#'   `consistent`. `grimmer` is evaluated only for SDs inside the sharp bounds;
#'   GRIM-inconsistent means are pruned before any SD is tested, which is why
#'   they are absent from the grid.
#' @examples
#' # the full grid of reportable (mean, sd) pairs for a small design
#' grid <- umbrella_data(n = 12, l = 1, u = 3, digits = 1)
#'
#' # how many pairs survive every test
#' table(grid$consistent)
#' @export
umbrella_data <- function(
  n,
  l,
  u,
  digits = 2,
  granularity = "integer",
  scoring = "singleitem",
  n_items = 1,
  alpha = NULL,
  rounding = "up_or_down"
) {
  step <- 10^(-digits)
  half_step <- step / 2
  # round() the grids back onto the nearest double to each decimal. seq() by a
  # decimal step accumulates error (seq(0, 6, by = 0.1)[4] is 0.3 + 5.6e-17),
  # and the granularity tests read the value as a decimal string, so an
  # off-by-one-ulp grid point can flip a GRIMMER verdict.
  means <- round(seq(l, u, by = step), digits)
  use_grimmer <- granularity == "integer"

  # the column shape a design with no feasible mean still has to return
  empty <- tibble::tibble(
    mean = numeric(0),
    sd = numeric(0),
    min_sd = numeric(0),
    max_sd = numeric(0),
    in_bounds = logical(0),
    grimmer = logical(0),
    consistent = logical(0)
  )

  means |>
    purrr::map(function(mean_current) {
      # sd_bounds() embeds the rounding-aware GRIM prefilter: a mean whose
      # rounding interval admits no integer sum is infeasible under strict
      # granularity, so none of its SDs need testing.
      bounds <- sd_bounds(
        l = l,
        u = u,
        n = n,
        mean = mean_current,
        digits_mean = digits,
        rounding = rounding,
        granularity = granularity,
        scoring = scoring,
        n_items = n_items,
        alpha = alpha
      )
      if (!isTRUE(bounds$feasible) || is.na(bounds$max_sd)) {
        return(NULL)
      }

      sds <- round(
        seq(0, ceiling(bounds$max_sd / step) * step, by = step),
        digits
      )
      in_bounds <- (sds + half_step) >= bounds$min_sd - TOLERANCE &
        (sds - half_step) <= bounds$max_sd + TOLERANCE

      # GRIMMER (the expensive per-tuple test) runs only where the SD is inside
      # the sharp bounds; outside, the tuple is already inconsistent.
      grimmer <- rep(NA, length(sds))
      if (use_grimmer && any(in_bounds)) {
        grimmer[in_bounds] <- as.logical(scrutiny::grimmer(
          x = mean_current,
          sd = sds[in_bounds],
          n = n,
          digits_x = digits,
          digits_sd = digits,
          items = n_items,
          rounding = rounding
        ))
      }

      tibble::tibble(
        mean = mean_current,
        sd = sds,
        min_sd = bounds$min_sd,
        max_sd = bounds$max_sd,
        in_bounds = in_bounds,
        grimmer = grimmer,
        consistent = in_bounds &
          (if (use_grimmer) !is.na(grimmer) & grimmer else TRUE)
      )
    }) |>
    purrr::compact() |>
    purrr::list_rbind() |>
    dplyr::bind_rows(empty)
}
