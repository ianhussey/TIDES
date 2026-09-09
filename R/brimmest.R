# Exact certification of reported (mean, SD) tuples against the attainable
# lattice of strictly integer data. See R/plot_sd_region.R for the dynamic
# program that enumerates the lattice.

#' BRIMMEST: certify whether a reported (mean, SD) is attainable
#'
#' Bounds-Related Inconsistency of Means and Errors, Settled Test. An exact
#' possible / impossible certificate for reported summary statistics on a
#' bounded integer scale, obtained without reconstructing any dataset.
#'
#' [brim()] and [brimmer()] — like GRIM and GRIMMER before them — test
#' conditions that are *necessary but not sufficient*: failing one proves a
#' report impossible, but passing them all proves nothing. `brimmest()` is
#' necessary *and* sufficient. It enumerates the exact attainable `(mean, sd)`
#' lattice for the design, rounds it to the reporting precision, and asks
#' whether the reported tuple is a member. A hit means some integer sample on
#' `[l, u]` of size `n` rounds to exactly this report; a miss means none does.
#'
#' | test | asks | verdict |
#' | --- | --- | --- |
#' | [brim()] | is the reported mean attainable? | necessary only |
#' | [brimmer()] | and is the reported SD attainable with it? | necessary only |
#' | `brimmest()` | is the pair jointly attainable by integer data? | necessary and sufficient |
#'
#' This is the question `unsum::closure_generate()` answers by search, reached
#' analytically instead. No witness datasets are produced, but the cost depends
#' only on `l`, `u` and `n`.
#'
#' Two routes reach the same verdict, chosen automatically by design size. For
#' a small design the whole lattice is enumerated once, cached per `(l, u, n)`,
#' and every report matched against it. For a large design asked about a
#' handful of reports, the report itself pins both coordinates almost
#' completely, and a closed-form sandwich plus a pruned constructive search
#' settles it without the lattice — which is the only option on a wide scale,
#' where the state table exceeds `max_cells`.
#'
#' @section Rounding and the direction of proof:
#' A hit certifies possibility under whichever rounding rule produced it. A
#' miss certifies impossibility only *relative to the rules in `rounding`*: a
#' report unreachable under one convention may be reachable under another. The
#' default takes the union of `"half_up"` and `"half_down"`. Narrow it only
#' when the source's convention is actually known, since certifying against one
#' rule when the paper used another manufactures a false impossibility. Unlike
#' [brimmer()] there is no `"up_or_down"` option, because rounding the lattice
#' forward needs a definite direction.
#'
#' @param l,u Numeric scalars, the scale limits. Integer-valued (in mean-score
#'   units when `n_items > 1` and `scoring` is `"meanscored"`).
#' @param n Integer scalar, sample size.
#' @param mean,sd Numeric vectors of reported means and SDs, recycled against
#'   each other. All are certified against one lattice.
#' @param digits_mean,digits_sd Integer scalars, the reported decimal places of
#'   the mean and of the SD.
#' @param rounding Character vector of rounding rules to admit, from
#'   `"half_up"`, `"half_down"`, `"native"`, `"ceiling"`, `"floor"`,
#'   `"trunc"`, `"anti_trunc"`. A tuple is possible if it is reachable under
#'   any of them.
#' @param scoring "singleitem" (default), "sumscored", or "meanscored".
#' @param n_items Positive whole number of response items (default 1).
#' @param max_cells Guard on the lattice enumeration's state space; a design
#'   above it takes the targeted route instead.
#' @param search_budget Node budget for the constructive search. A report that
#'   exhausts it errors rather than returning a guess.
#' @return A tibble with one row per reported tuple: `mean`, `sd`, `possible`
#'   (logical), and `rules` — the rounding rules under which the tuple is
#'   reachable, comma-separated and `""` when none.
#' @seealso [brim()] and [brimmer()] for the closed-form screens to run first,
#'   and [sd_region_data()] for the lattice itself.
#' @examples
#' # a report a real 1-5 scale sample can produce
#' brimmest(
#'   l = 1, u = 5, n = 9, mean = 3.0, digits_mean = 1,
#'   sd = 1.0, digits_sd = 1
#' )
#'
#' # inside the SD bounds and passing GRIM and GRIMMER, yet no integer
#' # sample produces it: the residual blind spot of the closed-form screen
#' brimmest(
#'   l = 1, u = 5, n = 9, mean = 1.3, digits_mean = 1,
#'   sd = 0.9, digits_sd = 1
#' )
#'
#' # one lattice certifies many reports at once
#' brimmest(
#'   l = 1, u = 5, n = 9,
#'   mean = c(3.0, 1.3, 2.5), digits_mean = 1,
#'   sd = c(1.0, 0.9, 1.2), digits_sd = 1
#' )
#' @export
brimmest <- function(
  l,
  u,
  n,
  mean,
  digits_mean,
  sd,
  digits_sd,
  rounding = c("half_up", "half_down"),
  scoring = c("singleitem", "sumscored", "meanscored"),
  n_items = 1,
  max_cells = MAX_LATTICE_CELLS,
  search_budget = SEARCH_BUDGET
) {
  scoring <- rlang::arg_match(scoring)

  valid <- c(
    "half_up",
    "half_down",
    "native",
    "ceiling",
    "floor",
    "trunc",
    "anti_trunc"
  )
  bad <- setdiff(rounding, valid)
  if (length(bad)) {
    cli::cli_abort(c(
      "Unknown rounding rule{?s}: {.val {bad}}.",
      i = "Choose from {.val {valid}}."
    ))
  }
  if (!length(rounding)) {
    cli::cli_abort("At least one rounding rule is required.")
  }
  if (is.null(n) || n < 2) {
    cli::cli_abort("{.arg n} must be >= 2 for a sample SD.")
  }

  geometry <- scoring_geometry(scoring, as.integer(round(n_items)), l, u)

  n_reports <- max(length(mean), length(sd))
  mean <- rep(mean, length.out = n_reports)
  sd <- rep(sd, length.out = n_reports)

  # Two routes to the same verdict, chosen by how much work each implies.
  #
  # Targeted (see R/attainable-target.R for the states, R/certify-sandwich.R
  # for the decision): arithmetic and a constructive search over just the
  # states the report pins down. Cost is per tuple, so it wins for a handful of
  # reports, and it is the only route on wide scales.
  #
  # Lattice: enumerate every attainable pair once and match against it. Cost is
  # per design regardless of how many tuples are asked about, so it wins once
  # there are enough of them. A small lattice is so cheap to build, and cached,
  # that it beats the targeted route even for one tuple.
  span <- as.integer(round(geometry$multiplier * (u - l)))
  cells <- lattice_cells(l, u, n, geometry$multiplier)
  affordable <- !is.na(cells) && cells <= max_cells
  use_targeted <- !affordable ||
    (cells > 1e6 && n_reports * length(rounding) <= 64)

  hits <- if (use_targeted) {
    rounding |>
      purrr::map(function(rounding_current) {
        purrr::map_lgl(seq_len(n_reports), function(report_current) {
          targets <- target_states(
            l,
            u,
            n,
            geometry$multiplier,
            mean[report_current],
            sd[report_current],
            digits_mean,
            digits_sd,
            rounding_current
          )
          # Exhausting the search tree is a proof of impossibility, so a
          # verdict is withheld only when the node budget cuts it short.
          settled <- certify_fast(span, n, targets, budget = search_budget)
          if (is.na(settled)) {
            cli::cli_abort(c(
              "This report was not settled within \\
               {.code search_budget = {search_budget}} nodes.",
              i = "Raise {.arg search_budget}, or report to fewer decimal places."
            ))
          }
          settled
        })
      }) |>
      (\(columns) do.call(cbind, columns))()
  } else {
    # one lattice per design, reused across every reported tuple and cached so
    # that repeated calls on the same design pay for it once
    lattice <- lattice_cached(l, u, n, geometry$multiplier, max_cells)
    # Both sides are on the 10^-digits reporting grid, so compare them as whole
    # numbers of grid steps rather than as strings: exact, and it keeps a
    # multi-million-row lattice cheap to match against. Integers also sidestep
    # the negative zero some rounding rules return at 0, which compares equal
    # but formats differently.
    scale_mean <- 10^digits_mean
    scale_sd <- 10^digits_sd
    reported_mean_steps <- round(mean * scale_mean)
    reported_sd_steps <- round(sd * scale_sd)

    rounding |>
      purrr::map(function(rounding_current) {
        lattice_mean_steps <- round(
          round_reported(lattice$mean, digits_mean, rounding_current) *
            scale_mean
        )
        lattice_sd_steps <- round(
          round_reported(lattice$sd, digits_sd, rounding_current) * scale_sd
        )
        # same packing for both sides
        span_sd <- max(c(lattice_sd_steps, reported_sd_steps), 0) + 1
        (reported_mean_steps * span_sd + reported_sd_steps) %in%
          (lattice_mean_steps * span_sd + lattice_sd_steps)
      }) |>
      (\(columns) do.call(cbind, columns))()
  }

  dim(hits) <- c(n_reports, length(rounding))

  tibble::tibble(
    mean = mean,
    sd = sd,
    possible = rowSums(hits) > 0,
    rules = purrr::map_chr(
      seq_len(n_reports),
      \(report_current) {
        paste(rounding[hits[report_current, ]], collapse = ",")
      }
    )
  )
}

# ---- Batch certification -----------------------------------------------------

#' Certify many reported (mean, SD) rows exactly
#'
#' Applies [brimmest()] to each row of a data frame. Columns of `data` whose
#' names match [brimmest()] arguments are taken per row; any argument given in
#' `...` is a constant broadcast to every row. Supplying one name both ways is
#' an error. The counterpart of [brimmer_map()] for the exact test.
#'
#' Rows are grouped by design and each group certified in a single [brimmest()]
#' call, which is the point of the function: one attainable lattice serves
#' every report of a design and is built once. Repeated `(mean, sd)` pairs
#' within a design are computed once. Row order is preserved.
#'
#' Recognised per row (column or constant): `l`, `u`, `n`, `mean`, `sd`,
#' `digits_mean`, `digits_sd`, `scoring`, `n_items`. `rounding`, `max_cells`
#' and `search_budget` may be given only as constants: `rounding` is a *set* of
#' admitted rules rather than one value per row, and the other two are cost
#' controls on the call.
#'
#' @param data A data frame, one reported statistic set per row.
#' @param ... Constant arguments applied to all rows.
#' @param include_inputs If TRUE (default), returns `data` column-bound to the
#'   results; if FALSE, only the result columns, in the same row order.
#' @return A tibble with `possible` (logical) and `rules` (the rounding rules
#'   under which the row is reachable, comma-separated and `""` when none),
#'   optionally with `data` prepended. Unlike [brimmest()], the reported `mean`
#'   and `sd` are not echoed back: they are already columns of `data` when
#'   supplied that way, and repeating them would collide.
#' @seealso [brimmest()] for one design at a time, [brimmer_map()] for the
#'   closed-form screen to run first.
#' @examples
#' reports <- tibble::tibble(mean = c(3.0, 1.3, 2.5), sd = c(1.0, 0.9, 1.2))
#'
#' # one design, one lattice. Rows 2 and 3 are both impossible, for different
#' # reasons: row 3 already fails GRIM and GRIMMER, while row 2 clears every
#' # closed-form screen and is caught only here.
#' brimmest_map(reports, l = 1, u = 5, n = 9, digits_mean = 1, digits_sd = 1)
#'
#' # designs may vary per row; rows sharing one are certified together
#' mixed <- tibble::tibble(
#'   mean = c(3.0, 1.3, 4.0),
#'   sd = c(1.0, 0.9, 1.5),
#'   n = c(9, 9, 12)
#' )
#' brimmest_map(
#'   mixed,
#'   l = 1, u = 5, digits_mean = 1, digits_sd = 1,
#'   include_inputs = FALSE
#' )
#' @export
brimmest_map <- function(data, ..., include_inputs = TRUE) {
  row_args <- c(
    "l",
    "u",
    "n",
    "mean",
    "sd",
    "digits_mean",
    "digits_sd",
    "scoring",
    "n_items"
  )
  call_args <- c("rounding", "max_cells", "search_budget")
  consts <- list(...)
  cols <- resolve_row_args(data, consts, row_args, call_args)
  present <- names(cols)

  if (!nrow(data)) {
    cli::cli_abort("{.arg data} has no rows.")
  }
  required <- c("l", "u", "n", "mean", "sd", "digits_mean", "digits_sd")
  missing <- setdiff(required, present)
  if (length(missing)) {
    cli::cli_abort(
      "{cli::qty(length(missing))}{.arg {missing}} {?is/are} required, \
       as a column of {.arg data} or a constant."
    )
  }

  passthrough <- consts[intersect(call_args, names(consts))]

  # Group by design, not by row. brimmest() amortises one lattice across every
  # tuple it is handed, and it also uses the number of tuples to choose its
  # route, so handing it a whole group at once is both cheaper and better
  # routed than calling it per row would be.
  design <- setdiff(present, c("mean", "sd"))
  design_key <- row_key(cols[design])

  possible <- logical(nrow(data))
  rules <- character(nrow(data))

  for (design_current in unique(design_key)) {
    rows <- which(design_key == design_current)
    args <- purrr::map(cols[design], \(column) column[rows[1]])
    means_group <- cols$mean[rows]
    sds_group <- cols$sd[rows]
    tuple_key <- row_key(list(means_group, sds_group))
    unique_tuples <- !duplicated(tuple_key)

    certified <- rlang::exec(
      brimmest,
      !!!args,
      mean = means_group[unique_tuples],
      sd = sds_group[unique_tuples],
      !!!passthrough
    )

    back <- match(tuple_key, tuple_key[unique_tuples])
    possible[rows] <- certified$possible[back]
    rules[rows] <- certified$rules[back]
  }

  out <- tibble::tibble(possible = possible, rules = rules)

  if (include_inputs) dplyr::bind_cols(data, out) else out
}
