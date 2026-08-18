# Exact certification of reported (mean, SD) tuples against the attainable
# lattice of strictly integer data. See R/plot_sd_region.R for the dynamic
# program that enumerates the lattice.

#' BRIMMEST: certify whether a reported (mean, SD) is attainable
#'
#' Bounds-Related Inconsistency of Means and Errors, Settled Test. An exact
#' possible / impossible certificate for reported summary statistics on a
#' bounded integer scale, obtained without reconstructing any dataset.
#'
#' `brimmest()` completes the family, and the superlative is the point.
#' [brim()] and [brimmer()] — like GRIM and GRIMMER before them — test
#' conditions that are *necessary but not sufficient*: failing one proves a
#' report impossible, but passing them all proves nothing. A handful of tuples
#' clear the SD bounds, GRIM and GRIMMER together and still correspond to no
#' integer dataset.
#'
#' `brimmest()` is *necessary and sufficient*. It enumerates the exact
#' attainable `(mean, sd)` lattice for the design by dynamic programming over
#' the reachable `(sum, sum of squares)` states, rounds that lattice to the
#' reporting precision, and asks whether the reported tuple is a member. The
#' verdict is therefore exact in both directions: a hit means some integer
#' sample on `[l, u]` of size `n` rounds to exactly this report, so it is
#' **possible**; a miss means no such sample exists, so it is **impossible** —
#' a proof, not a heuristic flag. Nothing in the GRIM family can make the
#' second claim.
#'
#' The three tests are nested, cheapest first:
#'
#' | test | asks | verdict |
#' | --- | --- | --- |
#' | [brim()] | is the reported mean attainable? | necessary only |
#' | [brimmer()] | and is the reported SD attainable with it? | necessary only |
#' | `brimmest()` | is the pair jointly attainable by real integer data? | necessary and sufficient |
#'
#' This is the same question the CLOSURE algorithm answers (see
#' `unsum::closure_generate()`), reached analytically rather than by search —
#' an *analytic CLOSURE certification*. The trade is that no witness datasets
#' are produced: the lattice records which `(mean, sd)` pairs are reachable,
#' not how. In exchange the cost depends only on `l`, `u` and `n`, not on how
#' many datasets happen to satisfy the constraints. The certification
#' validation document in `validation/` reports the comparison against CLOSURE
#' in full.
#'
#' @section How the work is done:
#' Two routes reach the same verdict, chosen automatically by design size.
#'
#' For a small design the whole attainable lattice is enumerated once and the
#' reports matched against it, so a vector of tuples costs barely more than a
#' single one. The lattice is cached per `(l, u, n)`, so repeated calls on one
#' design pay for it once.
#'
#' For a large design asked about a handful of reports, enumerating everything
#' to answer a membership question is wasteful — and on a wide scale it is
#' impossible, since the state table outgrows `max_cells`. A reported tuple
#' pins both coordinates almost completely: the mean's rounding interval
#' admits only a few integer sums, and each of those pins the sum of squares
#' to a narrow window. That reduces certification to a question about one
#' sample sum at a time, which is answered in three escalating steps.
#'
#' A closed-form sandwich comes first. For a given sum the achievable sums of
#' squares lie between the clustered configuration and the Structure-S one,
#' and all share the sum's parity; a window meeting no admissible integer
#' between them is impossible on arithmetic alone. What survives goes to a
#' constructive search over non-increasing samples — the partitions of the sum
#' — with the same sandwich re-applied at every step to what is left to place.
#' Reaching a complete sample proves possibility and produces a witness;
#' exhausting the search proves impossibility. Only when neither happens
#' inside `search_budget` does the corridor sweep run: a dynamic program
#' confined to the states that can still reach a target, roughly a tenth of
#' the full table.
#'
#' The practical effect is that designs the lattice refuses outright become
#' answerable — a 0-63 inventory at `n = 50` needs a 3151 x 24801 table — and
#' that the overwhelming majority of single reports are settled without any
#' sweep at all.
#'
#' @section Rounding and the direction of proof:
#' A hit certifies possibility under whichever rounding rule produced it. A
#' miss certifies impossibility only *relative to the rules in `rounding`*,
#' because a report unreachable under one convention may be reachable under
#' another. The default takes the union of `"half_up"` and `"half_down"`,
#' which spans the usual ambiguity about how a source rounded its halves.
#' Narrow it only when the source's convention is actually known: certifying
#' against one rule when the paper used another manufactures a false
#' impossibility. Note that unlike [brimmer()] there is no `"up_or_down"`
#' option, since rounding the lattice forward needs a definite direction.
#'
#' @param l,u Numeric scalars, the scale limits. Integer-valued (in mean-score
#'   units when `n_items > 1` and `scoring` is `"meanscored"`).
#' @param n Integer scalar, sample size.
#' @param mean,sd Numeric vectors of reported means and SDs, recycled against
#'   each other. All are certified against one lattice.
#' @param digits Integer, the reported decimal places, used for both `mean`
#'   and `sd` unless the next two are given.
#' @param mean_digits,sd_digits Integer scalars or NULL; reported decimal
#'   places when the mean and SD are reported to different precision.
#' @param rounding Character vector of rounding rules to admit, from
#'   `"half_up"`, `"half_down"`, `"native"`, `"ceiling"`, `"floor"`,
#'   `"trunc"`, `"anti_trunc"`. A tuple is possible if it is reachable under
#'   any of them.
#' @param scoring "singleitem" (default), "sumscored", or "meanscored".
#' @param n_items Positive whole number of response items (default 1).
#' @param max_cells Guard on the dynamic program's state space (default
#'   `2e7`); the enumeration errors rather than exhausting memory above it.
#' @param search_budget Node budget for the constructive search (default
#'   `2e5`). Raising it lets the search settle more reports on its own;
#'   lowering it hands them to the corridor sweep sooner. The verdict is the
#'   same either way — only which route produces it changes.
#' @return A data.frame with one row per reported tuple: `mean`, `sd`,
#'   `possible` (logical), and `rules` — the rounding rules under which the
#'   tuple is reachable, comma-separated and `""` when none.
#' @seealso [brim()] and [brimmer()] for the closed-form screens to run
#'   first, and [sd_region_data()] for the lattice itself.
#' @examples
#' # a report a real 1-5 scale sample can produce
#' brimmest(l = 1, u = 5, n = 9, mean = 3.0, sd = 1.0, digits = 1)
#'
#' # inside the SD bounds and passing GRIM and GRIMMER, yet no integer
#' # sample produces it: the residual blind spot of the closed-form screen
#' brimmest(l = 1, u = 5, n = 9, mean = 1.3, sd = 0.9, digits = 1)
#'
#' # one lattice certifies many reports at once
#' brimmest(l = 1, u = 5, n = 9, digits = 1,
#'         mean = c(3.0, 1.3, 2.5), sd = c(1.0, 0.9, 1.2))
#' @export
brimmest <- function(
  l,
  u,
  n,
  mean,
  sd,
  digits = NULL,
  mean_digits = NULL,
  sd_digits = NULL,
  rounding = c("half_up", "half_down"),
  scoring = c("singleitem", "sumscored", "meanscored"),
  n_items = 1,
  max_cells = 2e7,
  search_budget = 2e5
) {
  scoring <- match.arg(scoring)
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
    stop(
      "unknown rounding rule(s): ",
      paste(bad, collapse = ", "),
      ". Choose from: ",
      paste(valid, collapse = ", ")
    )
  }
  if (!length(rounding)) {
    stop("at least one rounding rule is required")
  }

  md <- if (!is.null(mean_digits)) mean_digits else digits
  sdd <- if (!is.null(sd_digits)) sd_digits else digits
  if (is.null(md) || is.null(sdd)) {
    stop(
      "reported decimal places are required: give digits, or both ",
      "mean_digits and sd_digits"
    )
  }
  if (is.null(n) || n < 2) {
    stop("n must be >= 2 for a sample SD")
  }

  k <- as.integer(round(n_items))
  g <- .scoring_geometry(scoring, k, l, u)

  nn <- max(length(mean), length(sd))
  mean <- rep(mean, length.out = nn)
  sd <- rep(sd, length.out = nn)

  # Two routes to the same verdict, chosen by how much work each implies.
  #
  # Targeted (see R/attainable-target.R): reachability of just the states the
  # report pins down. Cost is per tuple, so it wins for a handful of reports,
  # and it is the only route on wide scales, where the full lattice exceeds
  # its own size guard.
  #
  # Lattice: enumerate every attainable pair once and match against it. Cost
  # is per design regardless of how many tuples are asked about, so it wins
  # once there are enough of them, as when brimmest() is handed a whole grid.
  # A small lattice is so cheap to build, and now cached, that it beats the
  # targeted route even for one tuple; the targeted route earns its overhead
  # only once the lattice is genuinely large, and is the sole option once the
  # lattice will not fit at all.
  W <- as.integer(round(g$mg * (u - l)))
  cells <- .lattice_cells(l, u, n, g$mg)
  affordable <- !is.na(cells) && cells <= max_cells
  use_target <- !affordable ||
    (cells > 1e6 && nn * length(rounding) <= 64)

  if (use_target) {
    hit_mat <- vapply(
      rounding,
      function(rr) {
        vapply(
          seq_len(nn),
          function(i) {
            tg <- .target_states(l, u, n, g$mg, mean[i], sd[i], md, sdd, rr)
            # Arithmetic and a pruned constructive search first (see
            # R/certify-sandwich.R); the corridor sweep only for what they leave
            # undecided. Both give the same verdict, so the order is purely a
            # matter of cost.
            got <- .certify_fast(W, n, tg, budget = search_budget)
            if (is.na(got)) {
              got <- .attainable_target(W, n, tg, max_cells = max_cells)
            }
            if (is.na(got)) {
              stop(
                "this design is too large to certify at the requested precision; ",
                "raise max_cells, or report to fewer decimal places"
              )
            }
            got
          },
          logical(1)
        )
      },
      logical(nn)
    )
  } else {
    # one lattice per design, reused across every reported tuple and cached so
    # that repeated calls on the same design pay for it once
    lat <- .lattice_cached(l, u, n, g$mg, max_cells)
    # Both sides are on the 10^-digits reporting grid, so compare them as
    # whole numbers of grid steps rather than as strings: exact, and it keeps
    # a multi-million-row lattice cheap to match against. Integers also
    # sidestep the negative zero some rounding rules return at 0, which
    # compares equal but formats differently.
    mult_m <- 10^md
    mult_s <- 10^sdd
    qm <- round(mean * mult_m)
    qs <- round(sd * mult_s)
    hit_mat <- vapply(
      rounding,
      function(rr) {
        lm <- round(.round_reported(lat$mean, md, rr) * mult_m)
        ls <- round(.round_reported(lat$sd, sdd, rr) * mult_s)
        span <- max(c(ls, qs), 0) + 1 # same packing for both sides
        (qm * span + qs) %in% (lm * span + ls)
      },
      logical(nn)
    )
  }
  dim(hit_mat) <- c(nn, length(rounding))

  rules <- apply(hit_mat, 1L, function(z) paste(rounding[z], collapse = ","))
  data.frame(
    mean = mean,
    sd = sd,
    possible = as.logical(rowSums(hit_mat) > 0),
    rules = rules,
    stringsAsFactors = FALSE
  )
}

# ---- Batch certification -----------------------------------------------------

#' Certify many reported (mean, SD) rows exactly
#'
#' Applies [brimmest()] to each row of a data frame. Columns of `data` whose
#' names match [brimmest()] arguments are taken per row; any argument given in
#' `...` is a constant broadcast to every row. Supplying one name both ways is
#' an error. The counterpart of [brimmer_multiple()] for the exact test.
#'
#' Rows are **grouped by design** and each group certified in a single
#' [brimmest()] call, which is the whole point of the function: one attainable
#' lattice serves every report of a design and is built once, so certifying a
#' 200-row table on one design costs barely more than certifying one row of it.
#' Repeated `(mean, sd)` pairs within a design are computed once and reused.
#' Row order is preserved.
#'
#' Recognised per row (column or constant): `l`, `u`, `n`, `mean`, `sd`,
#' `digits`, `mean_digits`, `sd_digits`, `scoring`, `n_items`.
#'
#' `rounding`, `max_cells` and `search_budget` may be given only as constants,
#' never as columns. `rounding` is a *set* of admitted rules rather than one
#' value per row, and the other two are cost controls on the call; passing any
#' of them as a column is an error rather than a silent per-row reading.
#'
#' @param data A data frame, one reported statistic set per row.
#' @param ... Constant arguments applied to all rows.
#' @param include_inputs If TRUE (default), returns `data` column-bound to the
#'   results; if FALSE, only the result columns (same row order), for drop-in
#'   use inside a `dplyr::mutate()` / `purrr` pipeline.
#' @return A data frame with `possible` (logical) and `rules` (the rounding
#'   rules under which the row is reachable, comma-separated and `""` when
#'   none), optionally with `data` prepended. Unlike [brimmest()], the reported
#'   `mean` and `sd` are not echoed back: they are already columns of `data`
#'   when supplied that way, and repeating them would collide.
#' @seealso [brimmest()] for one design at a time, [brimmer_multiple()] for the
#'   closed-form screen to run first.
#' @examples
#' reports <- data.frame(mean = c(3.0, 1.3, 2.5),
#'                       sd   = c(1.0, 0.9, 1.2))
#'
#' # one design, one lattice. Rows 2 and 3 are both impossible, for different
#' # reasons: row 3 already fails GRIM and GRIMMER, while row 2 clears every
#' # closed-form screen and is caught only here.
#' brimmest_multiple(reports, l = 1, u = 5, n = 9, digits = 1)
#'
#' # designs may vary per row; rows sharing one are certified together
#' mixed <- data.frame(mean = c(3.0, 1.3, 4.0),
#'                     sd   = c(1.0, 0.9, 1.5),
#'                     n    = c(9, 9, 12))
#' brimmest_multiple(mixed, l = 1, u = 5, digits = 1,
#'                   include_inputs = FALSE)
#' @export
brimmest_multiple <- function(data, ..., include_inputs = TRUE) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  row_args <- c(
    "l",
    "u",
    "n",
    "mean",
    "sd",
    "digits",
    "mean_digits",
    "sd_digits",
    "scoring",
    "n_items"
  )
  call_args <- c("rounding", "max_cells", "search_budget")
  consts <- list(...)
  unknown <- setdiff(names(consts), c(row_args, call_args))
  if (length(unknown)) {
    stop("unknown constant argument(s): ", paste(unknown, collapse = ", "))
  }
  clash <- intersect(call_args, names(data))
  if (length(clash)) {
    stop(
      paste(clash, collapse = ", "),
      " applies to the whole call, not to ",
      "one row; supply it as a constant rather than a column"
    )
  }

  N <- nrow(data)
  if (!N) {
    stop("data has no rows")
  }
  resolve <- function(nm) {
    incol <- nm %in% names(data)
    incon <- nm %in% names(consts) && !is.null(consts[[nm]])
    if (incol && incon) {
      stop(sprintf("'%s' supplied as both a column and a constant", nm))
    }
    if (incol) {
      data[[nm]]
    } else if (incon) {
      rep(consts[[nm]], length.out = N)
    } else {
      NULL
    }
  }
  cols <- lapply(row_args, resolve)
  names(cols) <- row_args
  present <- row_args[!vapply(cols, is.null, logical(1))]
  for (nm in c("l", "u", "n", "mean", "sd")) {
    if (!(nm %in% present)) {
      stop(sprintf("'%s' is required (as a column of data or a constant)", nm))
    }
  }
  if (!any(c("digits", "mean_digits", "sd_digits") %in% present)) {
    stop(
      "reported decimal places are required: give digits, or both ",
      "mean_digits and sd_digits"
    )
  }

  passthru <- consts[intersect(call_args, names(consts))]

  # Group by design, not by row. brimmest() amortises one lattice across every
  # tuple it is handed, and it also uses the number of tuples to choose its
  # route, so handing it a whole group at once is both cheaper and better
  # routed than calling it per row would be.
  design <- setdiff(present, c("mean", "sd"))
  fmt <- function(nm) format(cols[[nm]], nsmall = 6, trim = TRUE)
  gkey <- do.call(paste, c(lapply(design, fmt), sep = "\r"))

  possible <- logical(N)
  rules <- character(N)
  for (g in unique(gkey)) {
    ix <- which(gkey == g)
    args <- lapply(design, function(nm) cols[[nm]][ix[1]])
    names(args) <- design
    mu <- cols$mean[ix]
    sg <- cols$sd[ix]
    tk <- paste(
      format(mu, nsmall = 6, trim = TRUE),
      format(sg, nsmall = 6, trim = TRUE),
      sep = "\r"
    )
    uk <- !duplicated(tk)
    r <- do.call(brimmest, c(args, list(mean = mu[uk], sd = sg[uk]), passthru))
    back <- match(tk, tk[uk])
    possible[ix] <- r$possible[back]
    rules[ix] <- r$rules[back]
  }

  res <- data.frame(
    possible = possible,
    rules = rules,
    stringsAsFactors = FALSE
  )
  if (include_inputs) cbind(data, res) else res
}
