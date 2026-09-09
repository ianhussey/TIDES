#' strait: Truncation-Induced Dependency among Summary Statistics
#'
#' A forensic meta-science toolkit for checking whether reported means,
#' standard deviations and sample sizes measured on a bounded (truncated) scale
#' are mutually consistent. When a measure has a known minimum and maximum, the
#' mean it can take constrains the standard deviation that is arithmetically
#' possible: the two summary statistics are dependent.
#'
#' `sd_bounds()` computes the smallest and largest sample SD consistent with a
#' chosen set of constraints, in closed form. `brimmer()` turns those bounds
#' into a report-level verdict and `brimmer_map()` applies it across a data
#' frame; `brimmest()` settles a reported tuple exactly. Companion functions
#' trace the feasible envelope (`sd_bounds_curve()`), build the jointly GRIM-,
#' GRIMMER- and bounds-consistent grid (`umbrella_data()`), and visualise both
#' (`plot_sd_bounds()`, `plot_sd_bounds_pomp()`, `plot_umbrella()`).
#'
#' @seealso The package README and `vignette("strait")` for a worked example
#'   and the method background.
#' @keywords internal
"_PACKAGE"
