# Build-time constants. Named here rather than inlined so that the tolerances
# in particular are one decision rather than sixty.

# Slack on floating-point comparisons of quantities that are whole numbers,
# grid steps or bounds by construction. Any legitimate difference in this
# package is at least 1/n, orders of magnitude larger.
TOLERANCE <- 1e-9

# The same idea one step tighter, for comparisons between two quantities that
# were computed the same way and so differ only by accumulated dust.
TOLERANCE_TIGHT <- 1e-12

# Guard on the attainable lattice's state table, in cells. A design above it
# takes brimmest()'s targeted route instead.
MAX_LATTICE_CELLS <- 2e7

# Node budget for the constructive witness search. Measured across 2,719,364
# reporting-grid cells and six designs it was never reached, so it is a safety
# valve rather than a tuning knob.
SEARCH_BUDGET <- 2e5

# Enumeration budget for the exact alpha-conditional Gini envelope, in integer
# sum-score profiles.
MAX_GINI_PROFILES <- 5e5

# Default density of a mean grid, as divisions of the scale width.
MEAN_GRID_DIVISIONS <- 1000

# Points in the dense part of the rounding-envelope candidate grid, before the
# analytic breakpoints are added to it.
ENVELOPE_GRID_POINTS <- 401

# Plot palette. Green and red are the consistency verdict; amber marks a tuple
# the bounds admit but GRIMMER does not; blue is the plain lattice point.
COLOUR_CONSISTENT <- "#43BF71"
COLOUR_INCONSISTENT <- "#D7191C"
COLOUR_GRIMMER <- "#FDAE61"
COLOUR_LATTICE <- "#1d4ed8"

# The rlang data pronoun, bound here so that tidy-eval code throughout the
# package can name columns as `.data$x` without an @importFrom directive. A
# data mask shadows this binding with its own pronoun; the binding exists so
# that R CMD check sees a definition and so that `strait::` calls need no
# attached rlang.
.data <- rlang::.data
