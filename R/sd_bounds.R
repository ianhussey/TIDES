# sd_bounds() dispatcher and its exact-mean core, plus the POMP transforms.
# Part of the strait package; see bounds-primitives.R for the bound formulas.

# ---- Layer 4: dispatcher -----------------------------------------------------

# Internal: bounds for one EXACT mean (or no mean), with sides already
# resolved. Returns a list; feasibility gates for band/GRIM are applied here.
.sd_bounds_core <- function(
  lower,
  upper,
  lower_att,
  upper_att,
  n,
  mean,
  Z,
  alpha,
  k_items
) {
  res <- list(
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
    band <- feasible_mean_band(lower, upper, lower_att, upper_att, n)
    if (mean < band[1] - 1e-9 || mean > band[2] + 1e-9) {
      return(fail(
        if (lower_att || upper_att) {
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
    identical(Z, "integer") && !is.null(mean) && !is_grim_consistent(mean, n)
  ) {
    return(fail(
      "mean is GRIM-inconsistent: no strictly integer sample has this mean"
    ))
  }

  # alpha branch (walls only; enforced by the dispatcher)
  if (!is.null(alpha)) {
    ab <- sd_bounds_alpha(lower, upper, n, mean, Z, alpha, k_items)
    if (!ab$feasible) {
      return(fail(
        "alpha floor exceeds alpha ceiling: no sample satisfies all constraints"
      ))
    }
    res$min_sd <- ab$min_sd
    res$max_sd <- ab$max_sd
    res$min_rule <- if (ab$min_sd > 0) ab$min_rule else "s >= 0"
    res$max_rule <- "alpha ceiling"
    if (!is.na(ab$note)) {
      res$note <- ab$note
    }
    return(res)
  }

  # ceiling: min over applicable ceilings
  if (!is.null(lower) && !is.null(upper)) {
    if (is.null(n)) {
      res$max_sd <- sd_max_span(lower, upper)
      res$max_rule <- "span/sqrt(2) (n = 2 maximum)"
    } else if (is.null(mean)) {
      res$max_sd <- sd_max_span_n(lower, upper, n)
      res$max_rule <- "parity ceiling"
    } else {
      res$max_sd <- sd_max_structure_s(mean, n, lower, upper)
      res$max_rule <- "Structure S (sharp mean-conditional ceiling)"
    }
  }

  # floor: max over applicable floors
  bump <- function(candidate, rule) {
    if (candidate > res$min_sd + 1e-12) {
      res$min_sd <<- candidate
      res$min_rule <<- rule
    }
  }
  if (Z %in% c("integer", "quasiinteger") && !is.null(mean)) {
    fl <- if (Z == "integer") {
      sd_min_integer(mean, n)
    } else {
      sd_min_quasi_integer(mean, n)
    }
    bump(fl, sprintf("%s floor (range-free)", Z))
  }
  if (!is.null(n)) {
    if (lower_att && upper_att) {
      bump(
        sd_min_two_pin(lower, upper, n, mean, Z),
        if (is.null(mean)) {
          "two-pin floor W/sqrt(2(n-1))"
        } else {
          "two-pin attained floor"
        }
      )
    } else if (upper_att && !is.null(mean)) {
      bump(
        sd_min_one_pin(upper, n, mean, Z, side = "max"),
        "one-pin attained floor (max)"
      )
    } else if (lower_att && !is.null(mean)) {
      bump(
        sd_min_one_pin(lower, n, mean, Z, side = "min"),
        "one-pin attained floor (min)"
      )
    }
  }

  if (res$min_sd > res$max_sd + 1e-9) {
    return(fail("floor exceeds ceiling: no sample satisfies all constraints"))
  }
  res
}

# Internal: candidate exact means inside [m_lo, m_hi] for the envelope —
# a dense grid plus every analytic breakpoint of the piecewise bounds
# (multiples of 1/n, where both floors and Structure S can kink).
.candidate_means <- function(m_lo, m_hi, n, grid_n = 401) {
  cand <- seq(m_lo, m_hi, length.out = grid_n)
  if (!is.null(n)) {
    t_range <- ceiling(n * m_lo - 1e-9):floor(n * m_hi + 1e-9)
    if (length(t_range) && length(t_range) <= 5000) {
      cand <- c(
        cand,
        t_range / n,
        pmax(m_lo, pmin(m_hi, t_range / n + 1e-12)),
        pmax(m_lo, pmin(m_hi, t_range / n - 1e-12))
      )
    }
  }
  sort(unique(pmax(m_lo, pmin(m_hi, cand))))
}

#' Bounds of the sample SD under a chosen set of constraints
#'
#' Computes the smallest and largest sample standard deviations consistent with
#' whichever constraints are supplied, following the nested framework and the
#' attained-extremes extensions of the STRAIT article. All constraints default
#' to NULL; supplying more of them can only narrow the bounds. With none
#' supplied the only bound is `0 <= s < Inf`.
#'
#' Constraint semantics:
#' * `l`, `u` — logical scale limits (walls): observations lie in `[l, u]` but
#'   need not touch either limit.
#' * `a`, `b` — observed extremes (attained): at least one observation equals
#'   each. `a` supersedes `l` and `b` supersedes `u` (a wall beyond an attained
#'   extreme is vacuous). Attainment leaves ceilings unchanged but creates
#'   nonzero floors.
#' * `n` — sample size (required for any mean-conditional bound).
#' * `mean`, `mean_digits` — the reported mean, numeric, with its number of
#'   reported decimal places. This mirrors scrutiny's post-string API (numeric
#'   `x` plus `digits_x`); a character-parsing wrapper can be layered above via
#'   [infer_digits()], but is not required here. With `rounding = NULL` the
#'   mean is treated as exact and `mean_digits` is ignored.
#' * `rounding` — NULL (mean exact) or one of `"up_or_down"`, `"up"`,
#'   `"down"`, `"even"`, `"ceiling"`, `"floor"`, `"trunc"`, `"anti_trunc"`
#'   (see [unround_interval()]; requires `mean_digits`). The returned bounds
#'   are then the ENVELOPE over all exact means consistent with the report,
#'   intersected with the feasible mean band. Under `Z = "integer"` the
#'   envelope enumerates the GRIM-consistent means in the interval exactly.
#' * `sd`, `sd_digits` — optionally, the reported SD and its decimal places.
#'   The SD never changes the bounds; it enables the report-consistency
#'   columns: `sd_in_bounds` (does the SD's own rounding interval overlap
#'   `[min_sd, max_sd]`?) and, under `Z = "integer"` with `rounding` set, the
#'   GRIMMER verdict.
#' * `Z` — granularity: `"continuous"` (none; the default), `"integer"` (all
#'   observations on the grid), or `"quasiinteger"` (all but one on the grid;
#'   defined at every mean, the GRIM-free relaxation). For mean-scored data the
#'   grid is the 1/`n_items` mean-score lattice, not the integers.
#' * `scoring`, `n_items` — how the `n_items` response items enter the scale.
#'   `"singleitem"` (default, requires `n_items = 1`): observations are single
#'   integer responses. `"sumscored"`: observations are sums of `n_items`
#'   integer items, so still integer-gridded; `n_items` is used only to size a
#'   reported `alpha`. `"meanscored"`: observations are means of `n_items`
#'   integer items, so they sit on a 1/`n_items` grid; `n_items` sets that
#'   granularity (and, with `alpha`, doubles as the composite item count).
#' * `alpha` — reported Cronbach's alpha of the composite; requires
#'   `scoring = "sumscored"` or `"meanscored"` (a single item has no internal
#'   consistency) and `l`, `u`, `mean`, `n`. For `"sumscored"`, `l`, `u`, `mean`
#'   are in SUM-score units; for `"meanscored"`, in mean-score units. Not
#'   combinable with `a`/`b` (open problem).
#'
#' GRIM and GRIMMER verdicts are DEFERRED to the scrutiny package
#' (`scrutiny::grim()`, `scrutiny::grimmer()`) rather than reimplemented:
#' under `Z = "integer"` with `rounding` set, the `grim` column carries
#' scrutiny's verdict on the mean and (when `sd` is supplied) the `grimmer`
#' column scrutiny's verdict on the (mean, sd) pair. The bounds themselves are
#' computed from this package's own enumeration of GRIM-consistent means; when
#' scrutiny's verdict disagrees with that enumeration the disagreement is
#' surfaced in `note`. Divergences are expected in known edge cases: scrutiny's
#' GRIM/GRIMMER have documented floating-point boundary bugs (false passes
#' under `rounding = "up"`; false flags, including under the default rounding,
#' when a candidate sum sits exactly on an interval endpoint), so verify such a
#' mismatch against a full enumeration before attributing it to the bounds.
#'
#' @param l,u Numeric scalars or NULL, logical limits.
#' @param a,b Numeric scalars or NULL, observed (attained) extremes.
#' @param n Integer scalar or NULL, sample size (`n >= 2` when given).
#' @param mean Numeric scalar or NULL, the reported sample mean.
#' @param mean_digits Integer scalar or NULL, decimal places of the reported
#'   mean (required with `rounding`).
#' @param sd Numeric scalar or NULL, the reported sample SD (verdict columns
#'   only; never affects the bounds).
#' @param sd_digits Integer scalar or NULL, decimal places of the reported SD
#'   (required for GRIMMER and for interval-based `sd_in_bounds` when
#'   `rounding` is set).
#' @param rounding NULL (mean is exact) or a rounding/truncation rule; see
#'   [unround_interval()].
#' @param Z "continuous" (default), "integer", or "quasiinteger".
#' @param scoring "singleitem" (default), "sumscored", or "meanscored".
#' @param n_items Positive whole number of response items (default 1). Must be
#'   1 for `scoring = "singleitem"`.
#' @param alpha Numeric scalar or NULL, reported Cronbach's alpha (requires
#'   `scoring = "sumscored"` or `"meanscored"`).
#'
#' @return A one-row data.frame: `min_sd`, `max_sd` (Inf if unbounded),
#'   `feasible` (FALSE when the constraint set admits no sample, with
#'   `min_sd`/`max_sd` NA), `min_rule`, `max_rule` (which bound binds),
#'   `grim`, `grimmer` (scrutiny's verdicts; NA when not applicable),
#'   `sd_in_bounds` (NA when no `sd` supplied), `note`.
#'
#' @examples
#' sd_bounds()                                     # 0 <= s < Inf
#' sd_bounds(l = 1, u = 5)                         # range-only ceiling
#' sd_bounds(l = 1, u = 5, n = 9)                  # parity ceiling
#' sd_bounds(l = 1, u = 5, n = 9, mean = 2)        # sharp mean-conditional
#' sd_bounds(l = 1, u = 5, n = 9, mean = 2.44, Z = "quasiinteger")
#' sd_bounds(a = 2, b = 4, n = 9, mean = 2.44, Z = "integer")   # attained
#' sd_bounds(l = 1, b = 4, n = 9, mean = 3.5, Z = "integer")    # mixed
#'
#' # reported mean rounded to 2 dp: envelope over [2.965, 2.975];
#' # GRIM verdict deferred to scrutiny
#' sd_bounds(l = 1, u = 7, n = 30, mean = 2.97, mean_digits = 2,
#'           rounding = "up_or_down", Z = "integer")
#' # with a reported SD: adds scrutiny's GRIMMER verdict and the bounds check
#' sd_bounds(l = 1, u = 7, n = 30, mean = 2.97, mean_digits = 2,
#'           sd = 2.83, sd_digits = 2, rounding = "up_or_down", Z = "integer")
#' # truncated reporting
#' sd_bounds(l = 1, u = 7, n = 30, mean = 2.97, mean_digits = 2,
#'           rounding = "trunc", Z = "quasiinteger")
#'
#' # sum-scored 3-item composite (l, u, mean in sum units) with alpha
#' sd_bounds(l = 3, u = 15, n = 20, mean = 9, Z = "integer",
#'           scoring = "sumscored", n_items = 3, alpha = 0.8)
#' # the same composite reported as an item mean (l, u, mean in mean units)
#' sd_bounds(l = 1, u = 5, n = 20, mean = 3, Z = "integer",
#'           scoring = "meanscored", n_items = 3, alpha = 0.8)
#' @export
sd_bounds <- function(
  l = NULL,
  u = NULL,
  a = NULL,
  b = NULL,
  n = NULL,
  mean = NULL,
  mean_digits = NULL,
  sd = NULL,
  sd_digits = NULL,
  rounding = NULL,
  Z = c("continuous", "integer", "quasiinteger"),
  scoring = c("singleitem", "sumscored", "meanscored"),
  n_items = 1,
  alpha = NULL
) {
  out <- function(r) {
    data.frame(
      min_sd = r$min_sd,
      max_sd = r$max_sd,
      feasible = r$feasible,
      min_rule = r$min_rule,
      max_rule = r$max_rule,
      grim = r$grim,
      grimmer = r$grimmer,
      sd_in_bounds = r$sd_in_bounds,
      note = r$note,
      stringsAsFactors = FALSE
    )
  }

  # -- validate ----------------------------------------------------------------
  Z <- match.arg(Z)
  scoring <- match.arg(scoring)
  if (is.character(mean) || is.character(sd)) {
    stop(
      "mean and sd must be numeric (mirroring scrutiny's post-string API); ",
      "parse reported strings upstream, e.g. with infer_digits() + as.numeric()"
    )
  }
  if (
    is.null(n_items) ||
      length(n_items) != 1 ||
      n_items < 1 ||
      abs(n_items - round(n_items)) > 1e-9
  ) {
    stop("n_items must be a positive whole number")
  }
  n_items <- as.integer(round(n_items))
  if (scoring == "singleitem") {
    if (n_items != 1L) {
      stop("scoring = 'singleitem' requires n_items = 1")
    }
    if (!is.null(alpha)) {
      stop(
        "scoring = 'singleitem' cannot take alpha (a single item has no ",
        "internal consistency); use scoring = 'sumscored' or 'meanscored'"
      )
    }
  }
  if (!is.null(n) && n < 2) {
    stop("n must be >= 2 for a sample SD")
  }
  if (!is.null(mean) && is.null(n)) {
    stop("mean-conditional bounds require n")
  }
  if (!is.null(rounding) && is.null(mean)) {
    stop("rounding requires a mean")
  }
  if (!is.null(rounding) && is.null(mean_digits)) {
    stop("rounding requires mean_digits (decimal places of the reported mean)")
  }
  if (!is.null(sd) && !is.null(rounding) && is.null(sd_digits)) {
    stop("a reported sd with rounding requires sd_digits")
  }
  if (!is.null(l) && !is.null(u) && u <= l) {
    stop("need u > l")
  }
  if (!is.null(a) && !is.null(b) && b < a) {
    stop("need b >= a")
  }
  if (!is.null(a) && !is.null(l) && a < l) {
    stop("observed minimum a cannot lie below the wall l")
  }
  if (!is.null(b) && !is.null(u) && b > u) {
    stop("observed maximum b cannot lie above the wall u")
  }
  # granularity limits live on the reported scale; the 1/n_items mean-score grid
  # maps to integers under w = n_items * x, so reported limits must be integers.
  if (Z %in% c("integer", "quasiinteger")) {
    lims <- c(l, u, a, b)
    if (length(lims) && any(abs(lims - round(lims)) > 1e-9)) {
      stop("integer/quasiinteger constraints require integer-valued limits")
    }
  }

  # -- granularity multiplier (m) and alpha item count (k) ---------------------
  # m rescales the reported (mean-score) scale to the integer sum scale on which
  # the granularity and GRIM machinery operate: w = m * x. Single-item and
  # sum-scored observations already sit on an integer grid, so m = 1; mean-scored
  # data sit on a 1/n_items grid, so m = n_items. k is the item count the alpha
  # layer uses: sum-scored takes it as reported, and after w = m * x a mean score
  # becomes a k = n_items item sum, so the same value serves both the granularity
  # divisor and the composite item count (a mean-scored composite with alpha).
  m <- if (scoring == "meanscored") n_items else 1L
  k <- n_items
  if (!is.null(alpha)) {
    if (alpha >= 1) {
      stop("alpha must be < 1")
    }
    if (k == 1L) {
      # alpha is inert at k = 1
      alpha <- NULL
    } else if (is.null(l) || is.null(u) || is.null(n) || is.null(mean)) {
      stop("alpha bounds require l, u, n, and mean")
    } else if (!is.null(a) || !is.null(b)) {
      stop("alpha with attained extremes is not supported (open problem)")
    }
  }

  # -- resolve sides: attained supersedes wall ---------------------------------
  lower <- if (!is.null(a)) a else l
  upper <- if (!is.null(b)) b else u
  lower_att <- !is.null(a)
  upper_att <- !is.null(b)

  # scaled (w = m * x) copies for the integer-grid core; SD outputs divide by m.
  scl <- function(x) if (is.null(x)) NULL else x * m
  lower_w <- scl(lower)
  upper_w <- scl(upper)
  core_w <- function(mm) {
    r <- .sd_bounds_core(
      lower_w,
      upper_w,
      lower_att,
      upper_att,
      n,
      if (is.null(mm)) NULL else mm * m,
      Z,
      alpha,
      k
    )
    if (!is.na(r$min_sd)) {
      r$min_sd <- r$min_sd / m
    }
    if (!is.na(r$max_sd) && is.finite(r$max_sd)) {
      r$max_sd <- r$max_sd / m
    }
    r
  }

  # -- scrutiny verdicts (deferred, never reimplemented) -----------------------
  # GRIM and GRIMMER are report-level consistency tests: they apply when the
  # mean is a rounded report (rounding set) and the data are strictly integer.
  # scrutiny has documented floating-point boundary bugs; verdicts are reported
  # verbatim and any divergence from our own enumeration is surfaced in `note`.
  grim_v <- NA
  grimmer_v <- NA
  if (Z == "integer" && !is.null(rounding) && !is.null(mean)) {
    if (!requireNamespace("scrutiny", quietly = TRUE)) {
      stop(
        "the scrutiny package is required for GRIM/GRIMMER verdicts under Z = 'integer' with rounding"
      )
    }
    grim_v <- isTRUE(as.logical(unname(
      .grim_compat(
        x = mean,
        n = n,
        digits = mean_digits,
        items = n_items,
        rounding = rounding
      )
    ))[1])
    if (!is.null(sd)) {
      grimmer_v <- isTRUE(as.logical(unname(
        .grimmer_compat(
          x = mean,
          sd = sd,
          n = n,
          digits_x = mean_digits,
          digits_sd = sd_digits,
          items = n_items,
          rounding = rounding
        )
      ))[1])
    }
  }

  # helper: does the reported sd (as an interval if rounded) overlap the bounds
  sd_check <- function(min_sd, max_sd) {
    if (is.null(sd) || is.na(min_sd)) {
      return(NA)
    }
    if (!is.null(rounding) && !is.null(sd_digits)) {
      iv <- unround_interval(sd, sd_digits, rounding)
      iv$hi >= min_sd - 1e-9 && iv$lo <= max_sd + 1e-9
    } else {
      sd >= min_sd - 1e-9 && sd <= max_sd + 1e-9
    }
  }
  finish <- function(r, extra_note = NULL) {
    r$grim <- grim_v
    r$grimmer <- grimmer_v
    r$sd_in_bounds <- sd_check(r$min_sd, r$max_sd)
    if (!is.null(extra_note)) {
      r$note <- if (is.na(r$note)) {
        extra_note
      } else {
        paste(r$note, extra_note, sep = "; ")
      }
    }
    out(r)
  }

  # -- exact-mean path ---------------------------------------------------------
  if (is.null(rounding) || is.null(mean)) {
    r <- core_w(mean)
    return(finish(r))
  }

  # -- rounded/truncated mean: envelope over the rounding interval -------------
  iv <- unround_interval(mean, mean_digits, rounding)
  m_lo <- iv$lo
  m_hi <- iv$hi
  band <- feasible_mean_band(lower, upper, lower_att, upper_att, n)
  m_lo2 <- max(m_lo, band[1])
  m_hi2 <- min(m_hi, band[2])
  if (m_lo2 > m_hi2 + 1e-12) {
    return(finish(list(
      min_sd = NA_real_,
      max_sd = NA_real_,
      feasible = FALSE,
      min_rule = NA_character_,
      max_rule = NA_character_,
      note = sprintf(
        "no mean in the rounding interval [%.6g, %.6g] lies in the feasible band [%.6g, %.6g]",
        m_lo,
        m_hi,
        band[1],
        band[2]
      )
    )))
  }

  divergence <- NULL
  if (Z == "integer") {
    # exact: enumerate the GRIM-consistent means in the interval. Under the
    # 1/n_items mean grid a mean is GRIM-consistent when n * m * mean is an
    # integer, so enumerate over nm = n * m.
    nm <- n * m
    t_range <- ceiling(nm * m_lo2 - 1e-9):floor(nm * m_hi2 + 1e-9)
    cand <- t_range[
      t_range >= nm * m_lo2 - 1e-9 & t_range <= nm * m_hi2 + 1e-9
    ] /
      nm
    if (!iv$lo_incl) {
      cand <- cand[abs(cand - m_lo) > 1e-9]
    }
    if (!iv$hi_incl) {
      cand <- cand[abs(cand - m_hi) > 1e-9]
    }
    # surface any disagreement with scrutiny's deferred GRIM verdict
    if (!is.na(grim_v) && grim_v != (length(cand) > 0)) {
      divergence <- "scrutiny::grim verdict disagrees with the package's GRIM-mean enumeration (a known scrutiny floating-point boundary case)"
    }
    if (!length(cand)) {
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
    cand <- .candidate_means(m_lo2, m_hi2, n * m)
  }

  results <- lapply(cand, core_w)
  ok <- vapply(results, function(r) r$feasible, logical(1))
  if (!any(ok)) {
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
  results <- results[ok]

  mins <- vapply(results, function(r) r$min_sd, numeric(1))
  maxs <- vapply(results, function(r) r$max_sd, numeric(1))
  i_min <- which.min(mins)
  i_max <- which.max(maxs)
  finish(
    list(
      min_sd = mins[i_min],
      max_sd = maxs[i_max],
      feasible = TRUE,
      min_rule = sprintf(
        "%s (envelope over rounding interval)",
        results[[i_min]]$min_rule
      ),
      max_rule = sprintf(
        "%s (envelope over rounding interval)",
        results[[i_max]]$max_rule
      ),
      note = sprintf(
        "mean treated as %s-%s to %d dp: envelope over [%.6g, %.6g]%s",
        if (rounding %in% c("trunc", "anti_trunc", "ceiling", "floor")) {
          "truncated"
        } else {
          "rounded"
        },
        rounding,
        iv$digits,
        m_lo,
        m_hi,
        if (m_lo2 > m_lo || m_hi2 < m_hi) {
          sprintf(", clipped to feasible [%.6g, %.6g]", m_lo2, m_hi2)
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
#   pomp_mean       (mean - lower) / (upper - lower)  in [0, 1]
#   pomp_sd_parity  s / s_max_parity, the mean-agnostic parity (Popoviciu) max.
#                   A LINEAR rescaling: floor 0, mean-independent, so points that
#                   differ only in s stay ordered and the umbrella geometry is
#                   undistorted. Comparable across scales that share only their
#                   own (range, n).
#   pomp_sd_sharp   (s - min_sd) / (max_sd - min_sd), the position of s within the
#                   SHARP mean-conditional band actually returned for this
#                   constraint set (Structure S ceiling; quasi-integer/Bernoulli
#                   or pinned floor as applicable, not merely zero). Non-linear
#                   and mean-dependent, but "1" means attained-at-this-mean and it
#                   pools scales onto one relative-dispersion axis. NA without a
#                   mean (the band is then not mean-conditional) or a degenerate
#                   band.
# lower/upper are the EFFECTIVE limits (a supersedes l, b supersedes u).
.pomp_cols <- function(mean, sd, min_sd, max_sd, lower, upper, n, has_mean) {
  pm <- if (
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
  ps_parity <- if (!is.null(sd) && !is.na(parity_max) && parity_max > 0) {
    sd / parity_max
  } else {
    NA_real_
  }
  ps_sharp <- if (
    !is.null(sd) &&
      has_mean &&
      !is.na(min_sd) &&
      !is.na(max_sd) &&
      is.finite(max_sd) &&
      (max_sd - min_sd) > 1e-12
  ) {
    (sd - min_sd) / (max_sd - min_sd)
  } else {
    NA_real_
  }
  list(pomp_mean = pm, pomp_sd_parity = ps_parity, pomp_sd_sharp = ps_sharp)
}
