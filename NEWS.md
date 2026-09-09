# strait 0.5.0

## Breaking changes

* The `Z` argument is renamed `granularity` in `sd_bounds()`, `brimmer()`,
  `brim()`, `sd_bounds_curve()`, `umbrella_data()`, `sd_bounds_sample()`,
  `sd_bounds_alpha()`, `sd_min_two_pin()` and `sd_min_one_pin()`.

* The `mean_digits` and `sd_digits` arguments are renamed `digits_mean` and
  `digits_sd` throughout, matching the `digits_*` convention of the rest of
  the scrutinyverse.

* `brimmer_multiple()` and `brimmest_multiple()` are renamed `brimmer_map()`
  and `brimmest_map()`, matching `scrutiny::grim_map()` and
  `scrutiny::grimmer_map()`.

* `brimmest()` and `brimmest_map()` take `digits_mean` and `digits_sd`
  instead of a shared `digits`, and both are now required.

* `infer_digits()` is removed. It was an alias for
  `scrutiny::decimal_places()`; call that directly.

* `plot_sd_region(reference =)` is renamed `show_reference`.

* Every function that returned a data frame now returns a tibble.

* `sd_region_data()` returns a fixed set of columns for every rule: `mean`,
  `lo`, `hi`, `sd`. Band rules leave `sd` `NA`; the lattice rules leave `lo`
  and `hi` `NA`. The `type` attribute still says which pair to read.

* `band_polygon()` and `umbrella_contour()` return a zero-row tibble rather
  than `NULL` when nothing is feasible.

## Dependencies

* `scrutiny (>= 1.0.0)` is now required, and `R/scrutiny-compat.R` is gone
  with it. It existed only to dispatch between the pre- and post-1.0
  GRIM/GRIMMER interfaces.

* Rounding and unrounding now go through `roundwork` rather than `scrutiny`.
  `unround_interval()` corrects two rules where `roundwork::unround()` does
  not yet agree with roundwork's own rounding functions: `"even"` returns `NA`
  inclusion flags (read here as included, the conservative direction), and
  `"anti_trunc"` returns them inverted. Verdicts are unchanged.

* `cli`, `dplyr (>= 1.1.0)`, `purrr (>= 1.0.0)` and `tibble` are now imports;
  `ggplot2 (>= 3.4.0)` is version-pinned for `linewidth`.

## Internal

* Errors are signalled with `cli::cli_abort()` and enum arguments matched with
  `rlang::arg_match()`.

* Tolerances, budgets and plot colours are named constants in `R/aaa.R` rather
  than literals repeated across the sources.

* No verdict, bound or plotted region changes in this release; the refactor is
  checked output-for-output against 0.4.10 across the bounds, certification,
  curve, lattice, umbrella and ring builders.

# strait 0.4.10

## Bug fixes

* `plot_sd_bounds(shade = "outside")` drew the feasible region as thousands of
  one-column slivers instead of a single band, for any design with a dense kink
  lattice. `sd_bounds_curve()` now records its grid spacing in a `"step"`
  attribute and `plot_sd_bounds()` reads it, rather than inferring the spacing
  from the median gap between means, which fails once the kinks outnumber the
  uniform grid. Verdicts were never affected.

# strait 0.4.9

## New features

* `brimmest_multiple()` applies `brimmest()` across a data frame, grouping rows
  by design so that one attainable lattice serves every report of a design.
  `rounding`, `max_cells` and `search_budget` may be given only as constants.

* `brimmest()` gains a fast certification path that decides most single reports
  by a closed-form sandwich screen and a pruned constructive search, instead of
  sweeping a state space sized by the design. Verdicts are unchanged; the cost
  is now set by the report rather than by the scale width and sample size.

* `brimmest()` gains `search_budget`, bounding that search.

## Performance

* A single report on a 0-63 inventory at `n = 50` — a design the full lattice
  refuses outright — falls from 78 s to 0.63 ms. Grid-sized workloads are
  unaffected: a whole reporting grid still routes to the cached lattice.

## Validation

* The new path agreed with the lattice route on all 2,719,364 reporting-grid
  cells spanning six designs, two precisions and two rounding rules, and was
  never cut short by its node budget.

# strait 0.4.8

## Breaking changes

* `sd_bounds_check()` is renamed `brimmer()` and `sd_bounds_check_multiple()`
  to `brimmer_multiple()`, naming the SD-side bounds test by analogy with
  GRIM / GRIMMER.

* `brimmer()` no longer requires a reported `sd`. Supplying only a mean runs
  the mean-side tests and returns `NA` in the SD columns.

* `brimmer()` gains an `in_scale_range` column and no longer reports a
  granularity failure as a bounds failure. `feasibility` is now a residual,
  reserved for infeasibility no other test accounts for.

## New features

* `brimmest()` completes the `brim()` / `brimmer()` / `brimmest()` family with
  an exact possible / impossible certificate for reported `(mean, sd)` tuples
  on a bounded integer scale. Where `brim()` and `brimmer()` test conditions
  that are necessary but not sufficient, `brimmest()` is both, so it closes the
  blind spot the closed-form screen admits. Verified cell-for-cell against
  `unsum::closure_generate()` across six designs.

* `brim()` is the mean-side bounds test: is the reported mean attainable at
  all, given the scale limits, `n` and any attained extremes? It reports the
  feasible mean band alongside the verdict.

* `band_polygon()` turns a `(mean, lo, hi)` band into closed rings, and the
  outside-shading branches of `plot_sd_bounds()` and `plot_sd_region()` knock
  those rings out of a shaded panel rather than assembling the shading around
  them. Constraint sets that leave stretches with no feasible sample previously
  left those means unshaded, asserting that any SD at all was possible there.

* `plot_sd_bounds()` and `plot_sd_region()` gain `shade`, defaulting to
  `"outside"`, so both read the same way as
  `plot_sd_bounds_pomp(reference = "sharp")`: shaded means unreachable.

* `plot_sd_bounds()` gains `expand`, padding the plotted region as a proportion
  of the scale width rather than a fixed number of SD units.

* `plot_umbrella()` gains `style`, defaulting to `"points"`. The previous
  `"tiles"` view is retained for methods exposition.

* `sd_delta()` exposes the count-parity correction relating Muilwijk's
  mean-conditional ceiling to the sharp one.

* `rule = "muilwijk"` is an alias for `rule = "mean"`.

## Documentation

* Every exported function carries a running example.

* New validation document `validation/certification.qmd`, comparing CLOSURE
  reconstruction with analytic certification.

* `Language` is now `en-GB`, with `inst/WORDLIST` regenerated.

# strait 0.4.4

* The package was renamed from `tides` to `strait`, to avoid a collision with
  the existing CRAN package `Tides`. The exported API was unchanged. Entries
  below describe the package under its former name.

# strait 0.4.0

A ground-up rewrite of the engine around the closed-form standard-deviation
bounds derived in the STRAIT article. The API is new throughout.

## Breaking changes

* Removed `tides()`, `tides_df()`, `umbrella()`, `approximate_sd_bounds()`,
  `plot_tides()` and `plot_tides_relative()`.

* `sd_bounds()` has a new signature and semantics: logical scale limits `l`,
  `u` (formerly `min`, `max`), optional attained extremes `a`, `b`, a
  granularity argument, `scoring` with `n_items`, and an optional Cronbach's
  `alpha`. It returns bounds in closed form together with the binding rule for
  each.

* `plot_umbrella()` now takes the output of `umbrella_data()`.

## New features

* Bounds are sharp under many more constraint sets: attained observed extremes,
  the quasi-integer floor (defined at every mean, so bound curves are
  hole-free), the sharp Structure-S ceiling, and a reported Cronbach's alpha
  for sum- or mean-scored composites.

* `n_items` greater than 1 is fully supported via the granularity grid.

* Rounded and truncated reported values are handled explicitly
  (`unround_interval()` plus a `rounding` argument); GRIM and GRIMMER verdicts
  are deferred to `scrutiny`.

* `sd_bounds_check()` gives a single verdict, the failing tests, and
  percent-of-maximum-possible transforms. `sd_bounds_check_multiple()` applies
  it across a data frame.

* `sd_bounds_curve()` and `umbrella_data()` build the plotting data;
  `plot_sd_bounds()`, `plot_sd_bounds_pomp()` and `plot_umbrella()` draw it.

* The alpha-conditional floor is sharpened for strictly integer composites via
  the composite's Gini mean difference (`sd_min_alpha_gini()`), so a reported
  positive alpha yields a strictly positive minimum SD even at whole-number
  sum-score means.

* `plot_sd_region()` and `sd_region_data()` draw the feasible SD region for any
  one constraint set in the nested framework, reproducing the panels of the
  article's nested-constraints figure.

* The single-purpose bound primitives, `sd_bounds_sample()`,
  `sd_max_muilwijk()` and `v_max_alpha()` are exported and documented.

* `umbrella_data()` gains an `alpha` argument.

# strait 0.3.2

## Bug fixes

* Core functions no longer error when the package is installed: functions that
  were used but never imported are now fully qualified.

* `umbrella()` and the vignette work across `scrutiny` versions, detecting the
  installed GRIM/GRIMMER interface at run time.

* `plot_umbrella()` no longer warns about a missing `digits` column; it gains a
  `digits` argument.

## New features

* `tides_df()` falls back to the `tides()` defaults when the optional columns
  are absent, so an input data frame need not carry them.

## Documentation

* Rewrote the README, added a testthat suite, a GitHub Actions `R-CMD-check`
  workflow, `cran-comments.md`, `inst/WORDLIST` and `LICENSE.md`.
