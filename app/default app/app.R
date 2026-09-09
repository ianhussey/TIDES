# strait: a browser front end -------------------------------------------------
#
# One reported (mean, SD, n) on a bounded scale, checked against the smallest
# and largest sample standard deviations that scale can produce. All of the
# computation is done by the strait package - the bounds, the verdicts and the
# plots alike - so this app only handles input, layout, and presentation.
#
# Results can be read in the measure's own units or in the scale-free (POMP)
# units of the STRAIT article, against either of two reference denominators;
# see the "About" panel and `?plot_sd_bounds_pomp`.
#
# Hussey, I. (2026). strait: Bounds checks for reported summary statistics.
# https://doi.org/10.5281/zenodo.21439905
#
# Run with: shiny::runApp("app/default app")

library(shiny)
library(bslib)
library(ggplot2)
library(strait)

# The shading of the feasible region relies on sd_bounds_curve() recording the
# mean-grid spacing that plot_sd_bounds() needs, which arrived in 0.4.10. An
# older install draws the band as thousands of slivers instead of one region,
# and does so silently, so refuse to start rather than deploy a wrong figure.
if (utils::packageVersion("strait") < "0.4.10") {
  stop(
    "This app needs strait >= 0.4.10 (installed: ",
    utils::packageVersion("strait"),
    "). Update with remotes::install_github('ianhussey/strait').",
    call. = FALSE
  )
}

strait_blue <- "#1d5a8e"
strait_cream <- "#f4ede0"

strait_citation <- paste0(
  "Hussey, I. (2026). strait: Bounds checks for reported summary statistics. ",
  "[Computer software] https://github.com/ianhussey/strait ",
  "doi:10.5281/zenodo.21439905"
)

# ---- presentation -----------------------------------------------------------

fmt <- function(x, digits = 3) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) {
    return("—")
  }
  formatC(x, format = "f", digits = digits)
}

# The bound rules read "<rule> (envelope over rounding interval)" whenever the
# mean is unrounded, which is nearly always. The suffix is already implied by
# the rounding rule in the sidebar and doubles the width of a glance-level card,
# so it is dropped there and kept in full in the "Full output" table.
short_rule <- function(x) {
  if (length(x) != 1L || is.na(x)) {
    return("—")
  }
  sub(" \\(envelope over rounding interval\\)$", "", x)
}

pct <- function(x, digits = 1) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) {
    return("—")
  }
  paste0(formatC(x, format = "f", digits = digits), "%")
}

# A verdict is a triple: banner class, headline, gloss. "Consistent" is the
# settled term in the consistency-check literature, so the headline uses it;
# the gloss carries the epistemic asymmetry the word alone does not - a failure
# is a proof of impossibility, a pass only a failure to rule the report out
# (see brimmest() for the test that can make the positive claim).
verdict_of <- function(r) {
  failed <- if (nzchar(r$failed_tests)) {
    strsplit(r$failed_tests, ",", fixed = TRUE)[[1]]
  } else {
    character(0)
  }
  gloss <- c(
    in_scale_range = paste(
      "the reported mean lies outside the range attainable on this scale at",
      "this sample size, so no data set has it, at any SD"
    ),
    bounds = paste(
      "the reported SD lies outside the feasible band the reported mean",
      "implies, so no data set on this scale has both"
    ),
    grim = paste(
      "GRIM: no strictly integer sample of this size has a mean that rounds to",
      "the reported value"
    ),
    grimmer = paste(
      "GRIMMER: the reported SD is not attainable by integer data at this mean",
      "(a deferred scrutiny verdict, with documented false-flag cases)"
    ),
    feasibility = paste(
      "the constraint set admits no sample, for a reason none of the named",
      "tests accounts for"
    )
  )
  if (!length(failed)) {
    return(list(
      "strait-good",
      "CONSISTENT",
      paste(
        "The reported values pass every test applied. These tests are",
        "necessary but not sufficient, so consistent means the report could not",
        "be ruled out, not that a data set producing it is known to exist. Use",
        "the exact certificate to settle that."
      )
    ))
  }
  only_grimmer <- identical(failed, "grimmer")
  list(
    if (only_grimmer) "strait-warn" else "strait-bad",
    if (only_grimmer) "FLAGGED (GRIMMER ONLY)" else "INCONSISTENT",
    paste0(
      sub("^(.)", "\\U\\1", paste(unname(gloss[failed]), collapse = "; "),
          perl = TRUE), ". ",
      if (only_grimmer) {
        paste(
          "Feasibility and the bounds both pass, and scrutiny's GRIMMER has",
          "documented false-flag cases, so verify this before acting on it."
        )
      } else {
        "The reported values cannot all be correct as printed."
      }
    )
  )
}

verdict_banner <- function(r) {
  v <- verdict_of(r)
  div(
    class = paste("strait-banner", v[[1]]),
    div(class = "strait-verdict", v[[2]]),
    div(class = "strait-gloss", v[[3]])
  )
}

stat_cards <- function(cells) {
  div(
    class = "strait-stats",
    lapply(cells, function(x) {
      div(class = "strait-stat", div(class = "k", x[[1]]), div(class = "v", x[[2]]))
    })
  )
}

# ---- bound curves -----------------------------------------------------------

# sd_bounds_curve() is the package's own envelope, and it is what draws the
# curve wherever it can. Two cases it does not cover:
#
#  - attained extremes, for which it has no argument;
#  - a wide scale with a large n, where its grid has a point at every kink of
#    the floor's sawtooth (the 1/(n * items) lattice). At 0-100 with n = 200
#    that is 60,001 evaluations and 23 seconds, far more resolution than a
#    460-pixel panel can show and slow enough to make the app look hung.
#
# In both, the same grid is walked through sd_bounds() directly - thinning the
# kinks to CURVE_MAX_KINKS in the second case - and the POMP columns rebuilt
# exactly as sd_bounds_curve() builds them, so either branch returns an object
# plot_sd_bounds() and plot_sd_bounds_pomp() can draw. Only the drawing is
# approximated: the verdict is computed from the reported values themselves and
# never from this curve.
CURVE_MAX_KINKS <- 600L

bounds_curve <- function(p) {
  m <- if (identical(p$scoring, "meanscored")) p$k else 1L
  nm <- p$n * m
  # counted rather than enumerated: at a large n this lattice is the thing
  # being avoided, so it must not be materialised to find out how big it is
  thinned <- (floor(p$u * nm) - ceiling(p$l * nm) + 1) > CURVE_MAX_KINKS

  if (!p$use_ab && !thinned) {
    # `by` is left at the package default: at half a second for a typical
    # design it is cheap enough, and it resolves the floor's sawtooth exactly.
    # (sd_bounds_curve() records the spacing it used, so plot_sd_bounds() no
    # longer has to guess it, and a coarser `by` would be safe here too.)
    out <- sd_bounds_curve(
      l = p$l, u = p$u, n = p$n, granularity = p$Z,
      scoring = p$scoring, n_items = p$k, alpha = p$alpha
    )
    attr(out, "thinned") <- FALSE
    return(out)
  }

  # A strictly uniform grid: the kinks are sampled rather than landed on, which
  # is what "thinned" warns about below the plot, and the spacing is recorded in
  # the same "step" attribute sd_bounds_curve() sets, so plot_sd_bounds() can
  # tell a genuine gap in the band from an ordinary interval.
  means <- seq(p$l, p$u, length.out = 1200)
  out <- do.call(rbind, lapply(means, function(mu) {
    d <- sd_bounds(
      l = p$l, u = p$u, a = p$a, b = p$b, n = p$n, mean = mu,
      granularity = p$Z, scoring = p$scoring, n_items = p$k, alpha = p$alpha
    )
    data.frame(
      mean = mu, min_sd = d$min_sd, max_sd = d$max_sd, feasible = d$feasible
    )
  }))
  # the effective limits are the attained ones (a supersedes l, b supersedes u),
  # which is also what brimmer() normalises the reported point by
  parity_max <- sd_max_span_n(p$eff_l, p$eff_u, p$n)
  out$pomp_mean <- (out$mean - p$eff_l) / (p$eff_u - p$eff_l)
  out$parity_max <- parity_max
  out$ceil_parity <- out$max_sd / parity_max
  out$floor_parity <- out$min_sd / parity_max
  attr(out, "step") <- (p$u - p$l) / 1199
  attr(out, "thinned") <- thinned
  out
}

# brimmer()'s unrounding vocabulary against brimmest()'s forward-rounding one.
# "up_or_down" spans the usual ambiguity about how a source rounded its halves,
# so it maps to the union of the two half rules.
brimmest_rules <- function(rounding) {
  switch(rounding,
    up_or_down = c("half_up", "half_down"),
    up = "half_up",
    down = "half_down",
    even = "native",
    rounding
  )
}

# ---- UI ---------------------------------------------------------------------

strait_theme <- bs_theme(
  version = 5,
  primary = strait_blue,
  bg = "#ffffff",
  fg = "#1c1c1c",
  base_font = font_link(
    "Space Grotesk",
    href = paste0(
      "https://fonts.googleapis.com/css2?",
      "family=Space+Grotesk:wght@400;500;600;700&display=swap"
    )
  )
)

strait_css <- sprintf(
  "
  .strait-banner { border-radius: 10px; padding: 1rem 1.25rem; margin: 0 0 1rem 0;
                   border-left: 10px solid; }
  .strait-banner .strait-verdict { font-size: 1.5rem; font-weight: 700;
                                   letter-spacing: 0.02em; }
  .strait-banner .strait-gloss { margin-top: 0.35rem; font-size: 0.95rem; }
  .strait-good { background: #eef6ee; border-color: #2e7d32; }
  .strait-good .strait-verdict { color: #2e7d32; }
  .strait-bad  { background: #fdeeed; border-color: #b3261e; }
  .strait-bad .strait-verdict { color: #b3261e; }
  .strait-warn { background: #fdf6e6; border-color: #b26a00; }
  .strait-warn .strait-verdict { color: #b26a00; }
  .strait-stats { display: flex; flex-wrap: wrap; gap: 0.75rem; margin-bottom: 1rem; }
  .strait-stat { background: %s; border-radius: 8px; padding: 0.5rem 0.9rem;
                 min-width: 8rem; }
  .strait-stat .k { font-size: 0.75rem; text-transform: uppercase;
                    letter-spacing: 0.06em; color: #5b5b5b; }
  .strait-stat .v { font-size: 1.15rem; font-weight: 600; color: %s; }
  .strait-note { font-size: 0.9rem; color: #5b5b5b; }
  .navbar-brand img { margin-right: 0.5rem; }
  .shiny-table-output table td, .shiny-table-output table th {
    white-space: nowrap;
  }
  pre { background: #f7f7f7; border-radius: 8px; padding: 0.9rem;
        overflow-x: auto; }
  ",
  strait_cream,
  strait_blue
)

report_panel <- nav_panel(
  title = "Check a report",
  layout_sidebar(
    fillable = FALSE,
    sidebar = sidebar(
      width = 340,
      accordion(
        open = c("The design", "The reported values", "Units"),
        accordion_panel(
          "The design",
          numericInput("l", "Scale minimum", value = 1, step = 1),
          numericInput("u", "Scale maximum", value = 7, step = 1),
          numericInput("n", "n (sample size)", value = 30, min = 2, step = 1),
          selectInput(
            "Z",
            "Granularity of the responses",
            choices = c(
              "quasi-integer (all but one response whole)" = "quasiinteger",
              "integer (every response whole)" = "integer",
              "continuous" = "continuous"
            ),
            selected = "quasiinteger"
          ),
          selectInput(
            "scoring",
            "How the items form the score",
            choices = c(
              "a single item" = "singleitem",
              "a sum of items" = "sumscored",
              "a mean of items" = "meanscored"
            )
          ),
          conditionalPanel(
            "input.scoring != 'singleitem'",
            numericInput("n_items", "Number of items", value = 3, min = 2, step = 1),
            tags$p(
              class = "strait-note",
              "For a sum score the limits above are the composite's; for an",
              "item mean they are the item's, and the responses then sit on a",
              "1/items grid rather than on the integers."
            ),
            checkboxInput(
              "use_alpha",
              HTML("Condition on a reported Cronbach&rsquo;s &alpha;"),
              value = FALSE
            ),
            conditionalPanel(
              "input.use_alpha == true",
              numericInput(
                "alpha",
                HTML("Cronbach&rsquo;s &alpha;"),
                value = 0.80, min = -0.99, max = 0.99, step = 0.05
              ),
              tags$p(
                class = "strait-note",
                "Internal consistency caps how much the items can disagree,",
                "which lowers the ceiling and, for integer items, raises the",
                "floor. If the observed value is not reported, a defensible",
                "upper value for the scale (say .90 to .95) gives a",
                "conservative check."
              )
            )
          )
        ),
        accordion_panel(
          "The reported values",
          numericInput("mean", "Reported mean", value = 2.97, step = 0.01),
          numericInput("digits_mean", "Its decimal places", value = 2, min = 0, max = 8, step = 1),
          numericInput("sd", "Reported SD", value = 2.83, min = 0, step = 0.01),
          numericInput("digits_sd", "Its decimal places", value = 2, min = 0, max = 8, step = 1),
          selectInput(
            "rounding",
            "Rounding rule used by the authors",
            choices = c(
              "nearest, halves either way" = "up_or_down",
              "nearest, halves up" = "up",
              "nearest, halves down" = "down",
              "nearest, halves to even" = "even",
              "ceiling" = "ceiling",
              "floor" = "floor",
              "truncate (toward zero)" = "trunc",
              "anti-truncate (away from zero)" = "anti_trunc"
            )
          ),
          tags$p(
            class = "strait-note",
            "A reported value designates an interval, not a point. Admitting",
            "fewer rounding rules than the source may have used manufactures",
            "false impossibilities."
          )
        ),
        accordion_panel(
          "Attained extremes",
          checkboxInput(
            "use_ab",
            "Condition on an observed minimum and maximum",
            value = FALSE
          ),
          conditionalPanel(
            "input.use_ab == true",
            numericInput("a", "Observed minimum", value = 1, step = 1),
            numericInput("b", "Observed maximum", value = 7, step = 1),
            tags$p(
              class = "strait-note",
              "Taken to be attained: at least one observation equals each. That",
              "forces the sample to span the range, which lifts the floor and",
              "narrows the means reachable at all."
            )
          )
        ),
        accordion_panel(
          "Units",
          radioButtons(
            "units",
            "Report the results in",
            choices = c(
              "the measure's own units" = "native",
              "scale-free (standardised) units" = "standardised"
            )
          ),
          conditionalPanel(
            "input.units == 'standardised'",
            radioButtons(
              "reference",
              "Standardise the SD against",
              choices = c(
                "the Popoviciu maximum: the largest SD any sample of this n can have" = "parity",
                "the sharp quasi-integer band: 0 is its floor, 1 its ceiling" = "sharp"
              ),
              selected = "parity"
            )
          )
        )
      ),
      tags$hr(),
      uiOutput("download_ui")
    ),
    uiOutput("verdict_ui"),
    uiOutput("stats_ui"),
    navset_underline(
      nav_panel(
        "Feasible region",
        plotOutput("region_plot", height = "460px"),
        uiOutput("plot_note")
      ),
      nav_panel(
        "Full output",
        tags$p(
          class = "strait-note",
          "Every column brimmer() returns, plus the rules that bind each bound",
          "from sd_bounds()."
        ),
        tableOutput("detail_table")
      ),
      nav_panel(
        "Exact certificate",
        tags$p(
          class = "strait-note",
          "The bounds, GRIM and GRIMMER are all necessary but not sufficient:",
          "failing one proves a report impossible, passing them all proves",
          "nothing. brimmest() settles the question outright for strictly",
          "integer data, by enumerating the attainable (mean, SD) pairs of the",
          "design rather than by reconstructing any data set."
        ),
        actionButton("certify", "Certify this report", class = "btn-primary"),
        tags$p(
          class = "strait-note",
          style = "margin-top: 0.75rem;",
          "This is the slow step, and on a wide scale with a large n it can",
          "take several seconds."
        ),
        uiOutput("certificate_ui")
      )
    )
  )
)

about_panel <- nav_panel(
  title = "About",
  div(
    class = "container",
    style = "max-width: 52rem;",
    withMathJax(),
    tags$h3("Bounds checks for reported summary statistics"),
    tags$p(
      "strait is a forensic meta-science tool for auditing reported summary",
      "statistics measured on a bounded scale. When a measure has a known",
      "minimum and maximum, the reported mean constrains the standard",
      "deviation that is arithmetically possible: the two statistics are",
      tags$em("not"),
      "independent. Given a reported mean, SD and sample size, the package",
      "computes the smallest and largest SDs that could have produced that",
      "mean under the scale's limits, and flags a report whose SD falls",
      "outside them."
    ),
    tags$p(
      "Source code, documentation and the R package:",
      a(
        "https://github.com/ianhussey/strait",
        href = "https://github.com/ianhussey/strait",
        target = "_blank",
        rel = "noopener"
      )
    ),
    tags$h4("The bounds"),
    tags$p(
      "The constraints are nested, and each one can only tighten the bounds.",
      "Write \\(R = u - l\\) for the scale range and \\(\\bar{x}\\) for the",
      "reported mean."
    ),
    tags$p(tags$b("Scale limits alone."), "The widest a bounded sample can be",
           "is one observation at each end, which happens at \\(n = 2\\):"),
    tags$p("$$s \\le \\frac{R}{\\sqrt{2}}.$$"),
    tags$p(
      tags$b("Adding the sample size."), "Dispersion is maximised by half the",
      "observations at each limit, which only parity permits when \\(n\\) is",
      "even. Popoviciu (1935) stated both cases; the odd-\\(n\\) case is often",
      "credited to Petocz (2005). They combine into one closed form:"
    ),
    tags$p("$$s_{\\max}(l,u,n) = \\frac{R}{2}\\sqrt{\\frac{n}{n-1}\\cdot\\frac{2n^{2}-1+(-1)^{n}}{2n^{2}}}.$$"),
    tags$p(
      tags$b("Adding the mean."), "Muilwijk (1966) sharpened Popoviciu's bound",
      "by conditioning on the mean, in the inequality better known through its",
      "rediscovery by Bhatia and Davis (2000):"
    ),
    tags$p("$$s \\le \\sqrt{\\frac{n}{n-1}(u-\\bar{x})(\\bar{x}-l)}.$$"),
    tags$p(
      "That is a valid ceiling at every mean, but it is attained only where",
      "the counts it requires at each limit,",
      "\\(n_l = n(u-\\bar{x})/R\\) and \\(n_u = n(\\bar{x}-l)/R\\), come out",
      "whole. When they do not, the best feasible configuration places",
      "\\(\\lfloor n_l \\rfloor\\) observations at \\(l\\),",
      "\\(\\lfloor n_u \\rfloor\\) at \\(u\\), and exactly one at the interior",
      "point \\(x_r = \\{n_l\\}\\,l + \\{n_u\\}\\,u\\)  -  the configuration",
      "Mestdagh et al. (2018) call Structure \\(S\\)  -  and the sharp ceiling",
      "is Muilwijk's bound times a single count-parity factor \\(\\delta\\):"
    ),
    tags$p("$$s_{\\max}(l,u,n,\\bar{x}) = \\sqrt{\\frac{n}{n-1}(u-\\bar{x})(\\bar{x}-l)\\,\\delta},\\qquad \\delta = 1 - \\frac{n\\,\\{n_l\\}\\{n_u\\}}{n_l\\,n_u}.$$"),
    tags$p(
      tags$b("Adding granularity."), "Up to here the floor is zero. It rises",
      "above zero once the responses are whole numbers, because a mean with a",
      "fractional part cannot be matched by identical values.",
      tags$em("Quasi-integer"), "data (\\(\\mathbb{Z}_{n-1}\\): all but one",
      "response whole) is the relaxation the maximising configuration already",
      "makes, and it is the sharpest floor defined at",
      tags$em("every"), "mean. With \\(d = \\bar{x} - \\lfloor\\bar{x}\\rfloor\\)",
      "and \\(g = \\{n\\bar{x}\\}\\)  -  precisely the quantity GRIM tests, so",
      "\\(g = 0\\) exactly when the mean is GRIM-consistent:"
    ),
    tags$p("$$s_{\\min}(l,u,n,\\bar{x},\\mathbb{Z}_{n-1}) = \\sqrt{\\frac{n\\,d(1-d) - g(1-g)}{n-1}}.$$"),
    tags$p(
      "The scale limits drop out: this floor is range-free. The ceiling is",
      "unchanged, since Structure \\(S\\)'s only interior observation is the",
      "one the constraint leaves free. At a GRIM-consistent mean \\(g = 0\\)",
      "and the expression reduces to the strictly integer floor",
      "\\(\\sqrt{n\\,d(1-d)/(n-1)}\\); evaluating that strict form at a mean",
      "with \\(g \\ne 0\\) overstates the true minimum by exactly \\(g(1-g)\\)",
      "units of sum of squares, which errs toward false flags  -  the one",
      "failure mode a screening tool must not have."
    ),
    tags$p(
      tags$b("Adding a reported Cronbach's alpha."), "For a composite of",
      "\\(k\\) items, a reported \\(\\alpha\\) fixes how much of the composite",
      "variance the items can contribute independently, and so caps that",
      "variance from above. Writing \\(D = 1 - \\frac{k-1}{k}\\alpha\\), the",
      "ceiling falls to the internal-consistency one and, for integer items,",
      "the quasi-integer floor is amplified to",
      "\\(s_{\\min}/\\sqrt{D}\\). In the scale-free units below the",
      "\\(\\alpha\\)-ceiling is again a single universal curve:"
    ),
    tags$p("$$\\tilde{s}_{\\max}(p,\\alpha) = \\frac{2\\sqrt{p\\,(100-p)}}{\\sqrt{\\,k-(k-1)\\alpha\\,}}.$$"),
    tags$p(
      "That splits the panel into zones. Above the mean-conditional ceiling a",
      "report is impossible for any bounded data of this size, whatever the",
      "internal consistency. Between the two ceilings it is impossible",
      tags$em("given the reported"), "\\(\\alpha\\) - achievable only with",
      "lower internal consistency than was reported, which is itself an",
      "inconsistency among the reported statistics. Below the",
      "\\(\\alpha\\)-amplified floor it is likewise inconsistent with the",
      "reported \\(\\alpha\\), and because that floor is GRIM-free the check",
      "applies whatever the granularity status of the mean. Near each scale",
      "limit no composite with the reported \\(\\alpha\\) exists at all, so the",
      "region is drawn in more than one piece. A reported \\(\\alpha\\) cannot",
      "currently be combined with attained extremes: the two sharpen the",
      "bounds by different mechanisms and no joint form is known."
    ),
    tags$p(
      tags$b("Attained extremes."), "Reporting an observed minimum \\(a\\) and",
      "maximum \\(b\\) says more than the scale limits do: at least one",
      "observation sits at each. The sample is then forced to span \\(b - a\\),",
      "which lifts the floor sharply and restricts the reachable means to",
      "\\([a + (b-a)/n,\\; b - (b-a)/n]\\)."
    ),
    tags$h4("Scale-free (standardised) units"),
    tags$p(
      "Bounds stated in the measure's units are not comparable across",
      "instruments with different ranges. Re-expressing each statistic as a",
      "percent of the way from its minimum to its maximum possible value gives",
      "the POMP mean (Cohen et al., 1999) and the percent of maximum possible",
      "SD:"
    ),
    tags$p("$$p = 100\\cdot\\frac{\\bar{x}-l}{R},\\qquad \\tilde{s} = \\frac{200\\,s}{R}\\sqrt{\\frac{n-1}{n}}.$$"),
    tags$p(
      "Under this transformation \\(l\\), \\(u\\) and \\(n\\) absorb into the",
      "units and the mean-conditional ceiling becomes a single universal arch,",
      "the same curve for a 1-to-7 Likert item and a 0-to-100 slider:"
    ),
    tags$p("$$\\tilde{s}_{\\max}(p) = 2\\sqrt{p\\,(100-p)}.$$"),
    tags$p(
      "This is the", tags$b("first"), "of the two standardisations the app",
      "offers: the SD is divided by the parity-corrected Popoviciu maximum,",
      "the largest SD any sample of this size on this scale can have whatever",
      "its mean. It is a linear rescaling, so the geometry is undistorted and",
      "the Muilwijk arch stays visible as the ceiling beneath 1. The",
      tags$b("second"), "goes further and divides by the sharp quasi-integer",
      "band itself,"
    ),
    tags$p("$$\\frac{s - s_{\\min}(l,u,n,\\bar{x},\\mathbb{Z}_{n-1})}{s_{\\max}(l,u,n,\\bar{x}) - s_{\\min}(l,u,n,\\bar{x},\\mathbb{Z}_{n-1})} \\in [0, 1],$$"),
    tags$p(
      "so the feasible region becomes exactly the unit square and a point's",
      "height is its position within its own band. Each denominator absorbs",
      "more of the reported information than the last; whichever is used, a",
      "value above 1 is impossible given the information in the denominator.",
      "Integer granularity does",
      tags$em("not"), "carry over  -  it is a property of the original units,",
      "and rescaling does not erase it  -  and near the scale ends the",
      "denominators approach zero, so rounding error in the mean produces large",
      "swings in the ratio. The scale-free forms are therefore for",
      "communication and cross-scale comparison; the consistency check itself",
      "is always computed in raw units over the rounding intervals, as it is",
      "here."
    ),
    tags$h4("What a verdict does and does not prove"),
    tags$p(
      "The tests are", tags$em("necessary"), "but not", tags$em("sufficient."),
      "Failing one proves the report impossible: no data set on this scale has",
      "those statistics. Passing them all proves only that this test could not",
      "rule the report out. A small residue of reports clears the bounds, GRIM",
      "and GRIMMER together and still corresponds to no integer data set."
    ),
    tags$ul(
      tags$li(tags$b("brim()"), " - is the reported mean attainable? Necessary only."),
      tags$li(tags$b("brimmer()"), " - and is the reported SD attainable with it? Necessary only."),
      tags$li(tags$b("brimmest()"), " - is the pair jointly attainable by real integer data? Necessary and sufficient.")
    ),
    tags$p(
      "The 'Exact certificate' tab runs the last of these. It enumerates the",
      "attainable (mean, SD) lattice for the design analytically, so a miss is",
      "a proof of impossibility rather than a flag  -  the same certificate the",
      "CLOSURE algorithm provides, reached without reconstructing any data set."
    ),
    tags$h4("Citation"),
    tags$p(strait_citation),
    tags$h4("Links"),
    tags$ul(
      tags$li(a(
        "Source code and documentation on GitHub",
        href = "https://github.com/ianhussey/strait",
        target = "_blank",
        rel = "noopener"
      )),
      tags$li(a(
        "Archived release (doi:10.5281/zenodo.21439905)",
        href = "https://doi.org/10.5281/zenodo.21439905",
        target = "_blank",
        rel = "noopener"
      )),
      tags$li(a(
        "Report a bug",
        href = "https://github.com/ianhussey/strait/issues",
        target = "_blank",
        rel = "noopener"
      ))
    ),
    tags$h4("References"),
    tags$ul(
      tags$li("Bhatia, R., & Davis, C. (2000). A better bound on the variance. American Mathematical Monthly, 107(4), 353-357."),
      tags$li("Brown, N. J. L., & Heathers, J. A. J. (2017). The GRIM test. Social Psychological and Personality Science, 8(4), 363-369."),
      tags$li("Cohen, P., Cohen, J., Aiken, L. S., & West, S. G. (1999). The problem of units and the circumstance for POMP. Multivariate Behavioral Research, 34(3), 315-346."),
      tags$li("Mestdagh, M., Pe, M., Pestman, W., Verdonck, S., Kuppens, P., & Tuerlinckx, F. (2018). Sidelining the mean. Behavior Research Methods, 50(5), 1953-1972."),
      tags$li("Muilwijk, J. (1966). Note on a theorem of M. N. Murthy and V. K. Sethi. Sankhya B, 28, 183."),
      tags$li("Petocz, P. (2005). Upper bounds on the standard deviation. Teaching Statistics, 27(2), 62."),
      tags$li("Popoviciu, T. (1935). Sur les equations algebriques ayant toutes leurs racines reelles. Mathematica, 9, 129-145.")
    ),
    tags$p(
      class = "strait-note",
      sprintf("This app runs strait version %s.", utils::packageVersion("strait"))
    ),
    tags$p(
      class = "strait-note",
      "Code is MIT licensed © Ian Hussey. Images are CC BY 4.0."
    )
  )
)

ui <- page_navbar(
  title = tags$span(
    tags$img(src = "logo.png", height = 34, alt = "strait logo"),
    "strait::BRIMMER"
  ),
  window_title = "strait::BRIMMER",
  theme = strait_theme,
  # a fillable page pins the body to the viewport height, which clips the About
  # panel - longer than any viewport - with nothing able to scroll it
  fillable = FALSE,
  header = tags$head(tags$style(HTML(strait_css))),
  report_panel,
  about_panel
)

# ---- server -----------------------------------------------------------------

server <- function(input, output, session) {
  params <- reactive({
    num <- list(
      l = input$l, u = input$u, n = input$n, mean = input$mean, sd = input$sd,
      digits_mean = input$digits_mean, digits_sd = input$digits_sd
    )
    validate(need(
      all(vapply(num, function(x) length(x) == 1L && !is.na(x) && is.finite(x), logical(1))),
      "Every field needs a value."
    ))
    validate(need(num$u > num$l, "The scale maximum must exceed the scale minimum."))
    validate(need(
      num$n >= 2 && num$n <= 1e6,
      "A sample standard deviation needs n of at least 2, and at most a million."
    ))
    validate(need(num$sd >= 0, "A standard deviation cannot be negative."))
    validate(need(
      num$digits_mean >= 0 && num$digits_mean <= 8 &&
        num$digits_sd >= 0 && num$digits_sd <= 8,
      "Decimal places must be between 0 and 8."
    ))

    scoring <- input$scoring
    k <- if (identical(scoring, "singleitem")) 1L else as.integer(input$n_items)
    validate(need(
      length(k) == 1L && !is.na(k) && k >= 1,
      "The number of items must be a whole number of at least 1."
    ))
    validate(need(
      identical(scoring, "singleitem") || k >= 2,
      "A composite needs at least two items; otherwise choose 'a single item'."
    ))

    Z <- input$Z
    use_ab <- isTRUE(input$use_ab)
    a <- b <- NULL
    if (use_ab) {
      validate(need(
        all(vapply(list(input$a, input$b), function(x) length(x) == 1L && !is.na(x) && is.finite(x), logical(1))),
        "Enter both the observed minimum and the observed maximum."
      ))
      a <- input$a
      b <- input$b
      validate(need(b > a, "The observed maximum must exceed the observed minimum."))
      validate(need(
        a >= num$l && b <= num$u,
        "The observed extremes must lie within the scale limits."
      ))
    }

    # a reported alpha needs a composite to be about, and the package has no
    # combined form for alpha and attained extremes (an open problem)
    alpha <- NULL
    if (!identical(scoring, "singleitem") && isTRUE(input$use_alpha)) {
      validate(need(
        length(input$alpha) == 1L && !is.na(input$alpha) && is.finite(input$alpha),
        "Enter a value for Cronbach's alpha."
      ))
      alpha <- input$alpha
      validate(need(
        alpha > -1 && alpha < 1,
        "Cronbach's alpha must lie strictly between -1 and 1."
      ))
      validate(need(
        !use_ab,
        paste(
          "A reported alpha cannot be combined with attained extremes: the two",
          "sharpen the bounds by different mechanisms and no joint form is",
          "known. Clear one of them."
        )
      ))
    }

    if (Z %in% c("integer", "quasiinteger")) {
      lims <- c(num$l, num$u, a, b)
      validate(need(
        all(abs(lims - round(lims)) < 1e-9),
        paste(
          "Integer and quasi-integer granularity need whole-number limits.",
          "For a mean-scored composite these are the item's limits, not the",
          "composite's."
        )
      ))
    }

    c(num, list(
      granularity = Z, scoring = scoring, k = k, rounding = input$rounding,
      use_ab = use_ab, a = a, b = b, alpha = alpha,
      eff_l = if (use_ab) a else num$l,
      eff_u = if (use_ab) b else num$u,
      units = input$units,
      reference = if (identical(input$units, "standardised")) input$reference else "parity"
    ))
  })

  # the constraint arguments shared by every strait entry point below
  constraints <- function(p) {
    x <- list(
      l = p$l, u = p$u, n = p$n,
      digits_mean = p$digits_mean, digits_sd = p$digits_sd,
      rounding = p$rounding, granularity = p$Z, scoring = p$scoring, n_items = p$k
    )
    if (p$use_ab) {
      x$a <- p$a
      x$b <- p$b
    }
    if (!is.null(p$alpha)) {
      x$alpha <- p$alpha
    }
    x
  }

  # one row, in exactly the shape plot_sd_bounds() and plot_sd_bounds_pomp()
  # take as `points`
  report <- reactive({
    p <- params()
    do.call(
      brimmer_map,
      c(list(data = data.frame(mean = p$mean, sd = p$sd)), constraints(p))
    )
  })

  # the same bounds again, for the rules that name which constraint binds
  bounds <- reactive({
    p <- params()
    do.call(sd_bounds, c(list(mean = p$mean, sd = p$sd), constraints(p)))
  })

  curve <- reactive({
    bounds_curve(params())
  })

  output$verdict_ui <- renderUI({
    verdict_banner(report())
  })

  output$stats_ui <- renderUI({
    p <- params()
    r <- report()
    b <- bounds()

    if (identical(p$units, "native")) {
      loc <- c("reported mean", fmt(p$mean, p$digits_mean))
      val <- c("reported SD", fmt(p$sd, p$digits_sd))
      lo <- c("smallest possible SD", fmt(r$min_sd))
      hi <- c("largest possible SD", fmt(r$max_sd))
    } else if (identical(p$reference, "parity")) {
      pm <- sd_max_span_n(p$eff_l, p$eff_u, p$n)
      loc <- c("POMP mean", pct(100 * r$pomp_mean))
      val <- c("reported SD (% of max)", pct(100 * r$pomp_sd_parity))
      lo <- c("smallest possible (% of max)", pct(100 * r$min_sd / pm))
      hi <- c("largest possible (% of max)", pct(100 * r$max_sd / pm))
    } else {
      loc <- c("POMP mean", pct(100 * r$pomp_mean))
      val <- c("position in the band", fmt(r$pomp_sd_sharp))
      lo <- c("band floor", if (is.na(r$min_sd)) "—" else "0.000")
      hi <- c("band ceiling", if (is.na(r$max_sd)) "—" else "1.000")
    }

    cells <- list(loc, val, lo, hi,
                  c("binding floor", short_rule(b$min_rule)),
                  c("binding ceiling", short_rule(b$max_rule)))
    if (!is.null(p$alpha)) {
      cells <- c(cells, list(c("reported alpha", fmt(p$alpha, 2))))
    }
    if (!is.na(r$grim)) {
      cells <- c(cells, list(c("GRIM", if (isTRUE(r$grim)) "pass" else "fail")))
    }
    if (!is.na(r$grimmer)) {
      cells <- c(cells, list(c("GRIMMER", if (isTRUE(r$grimmer)) "pass" else "fail")))
    }
    stat_cards(cells)
  })

  output$region_plot <- renderPlot({
    p <- params()
    cur <- curve()
    pts <- report()
    validate(need(
      any(cur$feasible & is.finite(cur$max_sd)),
      "No mean on this scale admits any SD under these constraints, so there is no region to draw."
    ))
    title <- sprintf(
      "%g to %g, n = %d, %s%s",
      p$l, p$u, as.integer(p$n),
      switch(p$Z,
        quasiinteger = "quasi-integer responses",
        integer = "integer responses",
        continuous = "continuous responses"
      ),
      paste0(
        if (p$use_ab) sprintf(", extremes %g and %g attained", p$a, p$b) else "",
        if (!is.null(p$alpha)) {
          sprintf(", %d items with alpha = %s", p$k, fmt(p$alpha, 2))
        } else {
          ""
        }
      )
    )
    if (identical(p$units, "native")) {
      plot_sd_bounds(cur, points = pts, title = title)
    } else {
      # both scale-free panels are drawn on one square frame, so the two
      # denominators can be compared by eye. coord_fixed() replaces the
      # coordinate system plot_sd_bounds_pomp() already set, which ggplot2
      # announces; the replacement is the point, so the message is muted.
      # No limits are imposed with it: a reported SD above its denominator
      # belongs on the plot, not clipped off the top of it.
      suppressMessages(
        plot_sd_bounds_pomp(cur, points = pts, reference = p$reference,
                            title = title) +
          ggplot2::coord_fixed(ratio = 1)
      )
    }
  })

  output$plot_note <- renderUI({
    p <- params()
    r <- report()
    cur <- curve()
    base <- paste(
      "Shading marks what is ruled out, so the clear region is what the",
      "constraints allow. The point is the reported (mean, SD), green if the",
      "bounds admit it and red if not."
    )
    extra <- if (identical(p$units, "native")) {
      NULL
    } else if (identical(p$reference, "parity")) {
      paste(
        "Both axes are scale-free: the SD is divided by the parity-corrected",
        "Popoviciu maximum, the largest SD any sample of this size on this",
        "scale can have whatever its mean. That is a linear rescaling, so the",
        "geometry is undistorted and the ceiling appears as the",
        "mean-conditional (Muilwijk) arch beneath 1."
      )
    } else {
      paste(
        "The SD is divided by its own sharp quasi-integer band, so the feasible",
        "region is exactly the unit square and the point's height is its",
        "position within the band."
      )
    }
    dropped <- !identical(p$units, "native") &&
      is.na(if (identical(p$reference, "parity")) r$pomp_sd_parity else r$pomp_sd_sharp)
    tagList(
      tags$p(class = "strait-note", base, extra),
      if (dropped) {
        tags$p(
          class = "strait-note",
          tags$b("The reported point is not drawn here:"),
          "it has no position on this standardised scale, because the bounds",
          "it would be divided by are undefined or degenerate at this mean.",
          "Switch to the measure's own units to see it."
        )
      },
      if (!identical(p$units, "native")) {
        tags$p(
          class = "strait-note",
          "Scale-free units are for communication and cross-scale comparison.",
          "Near the scale ends the denominators approach zero, so rounding",
          "error in the mean swings the ratio wildly. The verdict above is",
          "computed in raw units over the rounding intervals either way."
        )
      },
      if (identical(p$Z, "integer")) {
        tags$p(
          class = "strait-note",
          "Under strictly integer responses the band exists only at the means",
          "an integer sum can produce, so the region is drawn with gaps.",
          "Quasi-integer granularity is the same constraint relaxed by one",
          "observation, which is defined at every mean and never rules out",
          "anything integer data could produce."
        )
      },
      if (isTRUE(attr(cur, "thinned"))) {
        tags$p(
          class = "strait-note",
          "This design puts a kink of the floor at every 1/(n x items) step,",
          "more than the panel can resolve, so the curve is drawn on a thinned",
          "grid. Only the drawing is approximate: the verdict above comes from",
          "the reported values themselves, not from this curve."
        )
      }
    )
  })

  output$detail_table <- renderTable(
    {
      p <- params()
      r <- report()
      b <- bounds()
      logical_chr <- function(x) if (is.na(x)) "—" else as.character(isTRUE(x))
      # renderTable() renders a non-ASCII character as a literal <u+2014>
      # pseudo-tag, which the browser then swallows as an unknown element, so
      # the em dash used elsewhere becomes an empty cell here
      plain <- function(x) ifelse(x == "—", "n/a", x)
      data.frame(
        field = c(
          "consistent", "failed tests", "smallest possible SD",
          "largest possible SD", "binding floor", "binding ceiling",
          "constraint set feasible", "mean within the attainable range",
          "reported alpha", "GRIM", "GRIMMER", "reported SD within the bounds",
          "POMP mean", "SD as a percent of the maximum",
          "position in the sharp band", "note"
        ),
        value = plain(c(
          logical_chr(r$consistent),
          if (nzchar(r$failed_tests)) r$failed_tests else "none",
          fmt(r$min_sd, 4),
          fmt(r$max_sd, 4),
          if (is.na(b$min_rule)) "—" else b$min_rule,
          if (is.na(b$max_rule)) "—" else b$max_rule,
          logical_chr(r$feasible),
          logical_chr(r$in_scale_range),
          if (is.null(p$alpha)) "—" else fmt(p$alpha, 2),
          logical_chr(r$grim),
          logical_chr(r$grimmer),
          logical_chr(r$sd_in_bounds),
          fmt(100 * r$pomp_mean, 1),
          fmt(100 * r$pomp_sd_parity, 1),
          fmt(r$pomp_sd_sharp, 3),
          if (is.na(r$note)) "—" else r$note
        )),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    },
    striped = TRUE,
    rownames = FALSE,
    colnames = FALSE,
    width = "100%"
  )

  # -- the exact certificate --------------------------------------------------

  # brimmest() certifies strictly integer data against one lattice, so it needs
  # whole-number limits and has no attained-extremes form.
  certifiable <- reactive({
    p <- params()
    lims <- c(p$l, p$u)
    if (p$use_ab) {
      return("Exact certification does not take attained extremes; clear that option to use it.")
    }
    if (any(abs(lims - round(lims)) > 1e-9)) {
      return("Exact certification needs whole-number scale limits, since it enumerates integer data.")
    }
    NULL
  })

  certificate <- eventReactive(input$certify, {
    p <- params()
    withProgress(
      message = "Certifying against the attainable lattice",
      detail = "closed-form sandwich, then a constructive search",
      value = 0.5,
      tryCatch(
        brimmest(
          l = p$l, u = p$u, n = p$n, mean = p$mean, sd = p$sd,
          digits_mean = p$digits_mean, digits_sd = p$digits_sd,
          rounding = brimmest_rules(p$rounding),
          scoring = p$scoring, n_items = p$k
        ),
        error = function(e) list(error = conditionMessage(e))
      )
    )
  })

  output$certificate_ui <- renderUI({
    blocked <- certifiable()
    if (!is.null(blocked)) {
      return(div(class = "strait-banner strait-warn", div(
        class = "strait-verdict",
        "NOT AVAILABLE"
      ), div(class = "strait-gloss", blocked)))
    }
    alpha_note <- if (!is.null(params()$alpha)) {
      tags$p(
        class = "strait-note",
        tags$b("Note:"),
        "brimmest() certifies against the attainable lattice of integer data",
        "alone; it does not use the reported alpha. A verdict of possible here",
        "therefore means possible for some integer sample, which may not be one",
        "with the reported internal consistency. The bounds test above is the",
        "one that uses alpha."
      )
    }
    if (input$certify == 0L) {
      return(alpha_note)
    }
    cert <- certificate()
    if (!is.null(cert$error)) {
      return(div(class = "strait-banner strait-warn", div(
        class = "strait-verdict",
        "CERTIFICATION FAILED"
      ), div(class = "strait-gloss", cert$error)))
    }
    possible <- isTRUE(cert$possible[1])
    tagList(
      alpha_note,
      div(
        class = paste("strait-banner", if (possible) "strait-good" else "strait-bad"),
        div(class = "strait-verdict", if (possible) "POSSIBLE" else "IMPOSSIBLE"),
        div(
          class = "strait-gloss",
          if (possible) {
            sprintf(
              paste(
                "Some integer sample of size %d on this scale rounds to exactly",
                "this report, under %s. This is a proof of possibility, not a",
                "failure to rule it out."
              ),
              as.integer(params()$n), cert$rules[1]
            )
          } else {
            paste(
              "No integer sample on this scale produces this report under any",
              "of the rounding rules admitted. This is a proof, not a flag.",
              "A miss is relative to those rules: if the source rounded some",
              "other way, widen the rounding rule in the sidebar and re-run."
            )
          }
        )
      ),
      stat_cards(list(
        c("reported mean", fmt(cert$mean[1], params()$digits_mean)),
        c("reported SD", fmt(cert$sd[1], params()$digits_sd)),
        c("verdict", if (possible) "possible" else "impossible"),
        c("rounding rules admitted", paste(brimmest_rules(params()$rounding), collapse = ", "))
      ))
    )
  })

  # -- download ---------------------------------------------------------------

  output$download_ui <- renderUI({
    # validate()'s shiny.silent.error inherits from error, so this hides the
    # button exactly while the inputs are unusable
    ok <- tryCatch({
      report()
      TRUE
    }, error = function(e) FALSE)
    if (!isTRUE(ok)) {
      return(div(class = "strait-note", "The results download appears once the inputs are valid."))
    }
    downloadButton("download_results", "Download results (.csv)")
  })

  output$download_results <- downloadHandler(
    filename = function() paste0("strait_results_", Sys.Date(), ".csv"),
    content = function(file) {
      p <- params()
      r <- report()
      b <- bounds()
      out <- data.frame(
        l = p$l, u = p$u, n = p$n,
        a = if (p$use_ab) p$a else NA_real_,
        b = if (p$use_ab) p$b else NA_real_,
        mean = p$mean, digits_mean = p$digits_mean,
        sd = p$sd, digits_sd = p$digits_sd,
        rounding = p$rounding, granularity = p$Z, scoring = p$scoring, n_items = p$k,
        alpha = if (is.null(p$alpha)) NA_real_ else p$alpha,
        consistent = r$consistent, failed_tests = r$failed_tests,
        min_sd = r$min_sd, max_sd = r$max_sd,
        min_rule = b$min_rule, max_rule = b$max_rule,
        feasible = r$feasible, in_scale_range = r$in_scale_range,
        grim = r$grim, grimmer = r$grimmer, sd_in_bounds = r$sd_in_bounds,
        pomp_mean = r$pomp_mean, pomp_sd_parity = r$pomp_sd_parity,
        pomp_sd_sharp = r$pomp_sd_sharp,
        note = r$note,
        strait_version = as.character(utils::packageVersion("strait")),
        citation = strait_citation,
        stringsAsFactors = FALSE
      )
      utils::write.csv(out, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
