# ggplot2 visualisations of the SD bounds, POMP-normalised bounds, and umbrella.

.bounds_point_layer <- function(pts, xvar, yvar) {
  if (is.null(pts[["consistent"]])) {
    pts$consistent <- NA
  }
  pts$.consistent <- as.character(pts$consistent)
  list(
    ggplot2::geom_point(
      data = pts,
      ggplot2::aes(
        x = .data[[xvar]],
        y = .data[[yvar]],
        fill = .data$.consistent
      ),
      shape = 21,
      colour = "black",
      size = 2.4,
      na.rm = TRUE
    ),
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#43BF71", "FALSE" = "#D7191C"),
      na.value = "grey50",
      name = "Consistent",
      labels = c("TRUE" = "consistent", "FALSE" = "inconsistent")
    )
  )
}

#' Plot SD bounds on the native scale
#'
#' The feasible SD band (floor to ceiling) against the mean, optionally with
#' reported points coloured by consistency (green/red outlined dots).
#'
#' `shade = "outside"` (default) shades the infeasible region and leaves the
#' feasible one clear, matching [plot_sd_bounds_pomp()] with
#' `reference = "sharp"`, so the two scales read the same way: shaded means
#' unreachable. `shade = "inside"` fills the feasible band instead, which was
#' the behaviour before this argument existed.
#'
#' @param curve Output of [sd_bounds_curve()]. Its `"step"` attribute, when
#'   present, is the mean-grid spacing used to tell a sampling gap from a
#'   genuine one under `shade = "outside"`; a curve built by hand, without that
#'   attribute, should be on a uniform grid, from which the spacing is inferred.
#' @param points Optional data.frame with `mean`, `sd`, and (optionally)
#'   `consistent`; e.g. the output of [brimmer_multiple()].
#' @param title Optional plot title.
#' @param fill,line_colour Band fill and outline colours. `fill` is used only
#'   by `shade = "inside"`; the infeasible shading has its own fixed grey.
#' @param shade `"outside"` (default) shades the infeasible region;
#'   `"inside"` fills the feasible band.
#' @param expand Padding around the plotted region, as a proportion of the
#'   scale width `u - l` rather than a fixed number of SD units, so that the
#'   margin looks the same on a 1-5 scale and a 0-100 one. Applied to both
#'   axes. The limits always stretch to include `points`, so an out-of-bounds
#'   report is never clipped out of view.
#' @return A ggplot object.
#' @examples
#' curve <- sd_bounds_curve(l = 1, u = 7, n = 30, by = 0.1)
#' plot_sd_bounds(curve, title = "Feasible SDs, 1-7 scale, n = 30")
#'
#' # overlay reported values, coloured by consistency (the second is
#' # above the ceiling, so it plots as inconsistent)
#' reports <- data.frame(mean = c(2.97, 3.51), sd = c(2.83, 3.50))
#' checked <- brimmer_multiple(reports, l = 1, u = 7, n = 30,
#'                             mean_digits = 2, sd_digits = 2)
#' plot_sd_bounds(curve, points = checked)
#'
#' # the previous look, with the feasible band filled
#' plot_sd_bounds(curve, points = checked, shade = "inside")
#' @export
plot_sd_bounds <- function(
  curve,
  points = NULL,
  title = NULL,
  fill = "grey85",
  line_colour = "grey30",
  shade = c("outside", "inside"),
  expand = 0.03
) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  shade <- match.arg(shade)
  cur <- curve[curve$feasible & is.finite(curve$max_sd), ]

  # Padding is a proportion of the scale width, not a fixed SD amount, so the
  # margin is visually constant across scales of very different widths.
  lo_m <- min(cur$mean)
  hi_m <- max(cur$mean)
  pad <- expand * (hi_m - lo_m)
  # Shading the outside needs finite limits, and finite limits will silently
  # clip an out-of-bounds point - exactly the case the plot exists to show -
  # so the ceiling of the view must account for the reported points as well.
  y_hi <- max(c(cur$max_sd, points$sd), na.rm = TRUE)

  p <- ggplot2::ggplot(cur, ggplot2::aes(x = .data$mean))
  if (shade == "inside") {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$min_sd, ymax = .data$max_sd),
        fill = fill
      )
  } else {
    # Shade the whole panel, then knock the feasible rings out of it. Doing it
    # this way rather than assembling the grey from side rectangles and ribbons
    # means a mean at which the band is undefined stays shaded by construction:
    # there is nothing there to knock out. Assembling the grey instead leaves
    # such a gap unshaded, which would assert that any SD at all is possible
    # there. See band_polygon().
    # band_polygon() needs the grid spacing to tell a sampling gap from a
    # genuine one, and sd_bounds_curve() records the spacing it used because
    # its grid is not uniform: it adds each kink of the 1/(n * n_items) lattice
    # plus a pair of neighbours 1e-9 away. Guessing the spacing from the means
    # is what this used to do, and it fails outright once the kinks outnumber
    # the uniform grid (about n * (u - l) > 333 at the default `by`): the
    # median then falls BELOW the plain grid spacing, so every ordinary
    # interval reads as a gap and the band is drawn as thousands of slivers.
    # The median remains the fallback for a hand-built curve, where a uniform
    # grid makes it right.
    step <- attr(curve, "step")
    rings <- band_polygon(
      data.frame(mean = cur$mean, lo = cur$min_sd, hi = cur$max_sd),
      by = if (!is.null(step)) {
        step
      } else if (nrow(cur) > 1) {
        stats::median(diff(cur$mean))
      } else {
        1
      }
    )
    p <- p +
      ggplot2::annotate(
        "rect",
        xmin = -Inf,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf,
        fill = "grey10",
        alpha = 0.12
      )
    if (!is.null(rings)) {
      p <- p +
        ggplot2::geom_polygon(
          data = rings,
          ggplot2::aes(x = .data$mean, y = .data$y, group = .data$ring),
          inherit.aes = FALSE,
          fill = "white"
        )
    }
  }
  p <- p +
    ggplot2::geom_line(ggplot2::aes(y = .data$max_sd), colour = line_colour) +
    ggplot2::geom_line(ggplot2::aes(y = .data$min_sd), colour = line_colour) +
    ggplot2::coord_cartesian(
      xlim = c(lo_m - pad, hi_m + pad),
      ylim = c(-pad, y_hi + pad),
      expand = FALSE
    ) +
    ggplot2::labs(x = "Mean", y = "SD", title = title) +
    ggplot2::theme_minimal()
  if (!is.null(points)) {
    p <- p + .bounds_point_layer(points, "mean", "sd")
  }
  p
}

#' Plot SD bounds on a percent-of-maximum-possible (POMP) scale
#'
#' `reference = "parity"` normalises every SD by the mean-agnostic parity
#' (Popoviciu) ceiling: a linear rescaling, so the Structure S ceiling appears
#' as a dome under 1 and the umbrella geometry is undistorted. `reference =
#' "sharp"` normalises each SD by its own sharp mean-conditional band, so the
#' feasible region is exactly the unit square and a point's height is its
#' position within the band (`pomp_sd_sharp`); regions outside `[0, 1]^2` are
#' shaded infeasible.
#'
#' @param curve Output of [sd_bounds_curve()] (used for the parity band).
#' @param points Optional data.frame with `pomp_mean` and `pomp_sd_parity` /
#'   `pomp_sd_sharp` and `consistent` (e.g. from [brimmer_multiple()]).
#' @param reference `"sharp"` (default) or `"parity"`.
#' @param title Optional plot title.
#' @return A ggplot object.
#' @examples
#' curve <- sd_bounds_curve(l = 1, u = 7, n = 30, by = 0.1)
#'
#' # "sharp" makes the feasible region exactly the unit square
#' plot_sd_bounds_pomp(curve)
#'
#' # "parity" is a linear rescaling, so the umbrella keeps its shape
#' plot_sd_bounds_pomp(curve, reference = "parity")
#' @export
plot_sd_bounds_pomp <- function(
  curve,
  points = NULL,
  reference = c("sharp", "parity"),
  title = NULL
) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  reference <- match.arg(reference)
  if (reference == "parity") {
    cur <- curve[curve$feasible & is.finite(curve$max_sd), ]
    p <- ggplot2::ggplot(cur, ggplot2::aes(x = .data$pomp_mean)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$floor_parity, ymax = .data$ceil_parity),
        fill = "grey85"
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = .data$ceil_parity),
        colour = "grey30"
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = .data$floor_parity),
        colour = "grey30"
      ) +
      ggplot2::labs(
        x = "Relative location (POMP mean)",
        y = "Relative dispersion (parity-normalised SD)",
        title = title
      )
    if (!is.null(points)) {
      p <- p + .bounds_point_layer(points, "pomp_mean", "pomp_sd_parity")
    }
  } else {
    shade <- data.frame(
      xmin = c(-Inf, 1, 0, 0),
      xmax = c(0, Inf, 1, 1),
      ymin = c(-Inf, -Inf, 1, -Inf),
      ymax = c(Inf, Inf, Inf, 0)
    )
    p <- ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = shade,
        ggplot2::aes(
          xmin = .data$xmin,
          xmax = .data$xmax,
          ymin = .data$ymin,
          ymax = .data$ymax
        ),
        fill = "grey10",
        alpha = 0.12
      ) +
      ggplot2::geom_rect(
        data = data.frame(x = 0),
        ggplot2::aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1),
        fill = NA,
        colour = "black",
        linewidth = 0.3
      ) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(-0.1, 1.1)) +
      ggplot2::labs(
        x = "Relative location (POMP mean)",
        y = "Position in sharp SD band",
        title = title
      )
    if (!is.null(points)) {
      p <- p + .bounds_point_layer(points, "pomp_mean", "pomp_sd_sharp")
    }
  }
  p + ggplot2::theme_minimal()
}

#' Plot the umbrella grid
#'
#' Renders the reportable `(mean, sd)` tuples of a design, in either of two
#' styles. Optionally overlays the bound curves from [sd_bounds_curve()].
#'
#' `style = "points"` (default) greys the whole panel and draws only the
#' consistent tuples. Nothing else is drawn, because nothing else exists: every
#' other cell of the grid is a value that cannot be reported. This matches the
#' convention of [plot_sd_bounds()] and [plot_sd_region()], where shading marks
#' what is ruled out, and it makes the real structure legible — the vertical
#' striping at means an integer sum can round to.
#'
#' `style = "tiles"` draws every cell of the reporting grid, coloured as
#' consistent, GRIMMER-inconsistent, or out of bounds. That is useful for
#' methods exposition, since it separates what the bounds rule out from what
#' GRIMMER additionally rules out, but it inverts the emphasis: at
#' `n = 14` on a 1-7 scale at two decimal places it spends two thirds of its ink, in the most
#' saturated colour, on impossible tuples, and GRIMMER's verdict alternating
#' between adjacent SDs produces interference banding that obscures the striping.
#'
#' `style = "contour"` draws only the outline of the umbrella: the panel is
#' shaded, the region between the lowest and highest consistent SD at each mean
#' is knocked out of the shading, and that region is stroked. It is the look of
#' [plot_sd_bounds()] — and of the app in `app/default app/` — but taken from
#' the lattice itself rather than from the continuous bounds, so what it outlines
#' is what the tests actually admit. Use it when the point cloud is too dense to
#' read, when the figure is a backdrop for reported points, or to see the
#' discrete envelope against the continuous one by passing `curve` as well.
#'
#' The outline itself comes from [umbrella_contour()], which is exported, so
#' the rings can also be had as data. Being an envelope, it is drawn across the
#' vertical striping rather than around it, and its interior therefore claims
#' less than the points do: see that function for what `contour_by` controls
#' and what the region does and does not assert.
#'
#' The grey under `"points"` and `"contour"` is a layer, not a theme element:
#' it is a rectangle over the whole panel with the feasible region knocked out
#' of it, so it says "no tuple here" as data rather than as decoration. That is
#' why no theme call removes it — `theme_void()` strips the axes and panel and
#' the shading stays, as it must, since the statement it makes survives the
#' loss of the axes. Use `shade = "none"` to drop it, which also leaves the
#' contour unfilled, so the plot can sit over something else.
#'
#' Note what the points represent. The consistent set is the GRIM- and
#' GRIMMER-consistent one, which is what [brimmer()] applies. It is strictly
#' larger than the set of attainable tuples, so an umbrella plot shows what the
#' test admits rather than what exists; use [brimmest()] to certify a tuple.
#' The gap is not always small — 9\% at `l = 1, u = 5, n = 10`, and nearly a
#' factor of three at `l = 0, u = 6, n = 23` with two items.
#'
#' @param umbrella Output of [umbrella_data()], or an already-filtered lattice
#'   from `sd_region_data(rule = "integer")` — anything with `mean` and `sd`,
#'   filtered by `consistent` when that column is present.
#' @param curve Optional [sd_bounds_curve()] output to overlay as bound lines.
#' @param title Optional plot title.
#' @param style `"points"` (default), `"tiles"` or `"contour"`; see Details.
#' @param point_colour,point_size Point appearance, used by `style = "points"`.
#'   `point_size` defaults to a value chosen from the number of points, since a
#'   size that reads well for a few hundred is a solid mass at twenty thousand.
#' @param reference_colour Colour of the overlaid bound curves.
#' @param line_colour Colour of the contour outline, used by
#'   `style = "contour"`.
#' @param contour_by Passed to [umbrella_contour()] as `by`: the mean spacing
#'   that `style = "contour"` uses to tell a genuine gap in the umbrella from
#'   its ordinary striping (see Details). Defaults to the median spacing between
#'   the means that have a consistent tuple; give it explicitly to keep a narrow
#'   gap — a stretch near a scale limit that a reported alpha rules out, say —
#'   from being bridged.
#' @param shade `"outside"` (default) greys everything the design rules out;
#'   `"none"` draws no shading at all, and leaves the contour unfilled. Applies
#'   to `"points"` and `"contour"`; `"tiles"` colours every cell already, so
#'   there is nothing left to shade.
#' @param expand Padding as a proportion of the scale width, as in
#'   [plot_sd_bounds()].
#' @return A ggplot object.
#' @examples
#' grid <- umbrella_data(n = 12, l = 1, u = 7, digits = 2)
#' plot_umbrella(grid, title = "n = 12, 1-3 scale")
#'
#' # The previous look, separating the two ways a tuple can be ruled out
#' plot_umbrella(grid, style = "tiles")
#'
#' # Just the outline of the region, as the Shiny app draws it, with the
#' # continuous bounds dashed over it for comparison
#' curve <- sd_bounds_curve(l = 1, u = 3, n = 12, by = 0.05)
#' plot_umbrella(grid, curve = curve, style = "contour")
#'
#' # The outline alone, on whatever background the theme gives it
#' plot_umbrella(grid, style = "contour", shade = "none") + ggplot2::theme_void()
#' @export
plot_umbrella <- function(
  umbrella,
  curve = NULL,
  title = NULL,
  style = c("points", "tiles", "contour"),
  point_colour = "black",
  point_size = NULL,
  reference_colour = "grey35",
  line_colour = "grey30",
  contour_by = NULL,
  shade = c("outside", "none"),
  expand = 0.03
) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  style <- rlang::arg_match(style)
  shade <- rlang::arg_match(shade)

  if (style == "tiles") {
    umbrella$category <- ifelse(
      umbrella$consistent,
      "consistent",
      ifelse(
        umbrella$in_bounds & !is.na(umbrella$grimmer) & !umbrella$grimmer,
        "GRIMMER-inconsistent",
        "out of bounds"
      )
    )
    p <- ggplot2::ggplot(
      umbrella,
      ggplot2::aes(x = .data$mean, y = .data$sd, fill = .data$category)
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_manual(
        values = c(
          "consistent" = "#43BF71",
          "GRIMMER-inconsistent" = "#FDAE61",
          "out of bounds" = "grey80"
        ),
        name = NULL
      ) +
      ggplot2::labs(x = "Mean", y = "SD", title = title) +
      ggplot2::theme_minimal()
    if (!is.null(curve)) {
      cur <- curve[curve$feasible & is.finite(curve$max_sd), ]
      p <- p +
        ggplot2::geom_line(
          data = cur,
          ggplot2::aes(x = .data$mean, y = .data$max_sd),
          inherit.aes = FALSE,
          colour = "grey20"
        ) +
        ggplot2::geom_line(
          data = cur,
          ggplot2::aes(x = .data$mean, y = .data$min_sd),
          inherit.aes = FALSE,
          colour = "grey20"
        )
    }
    return(p)
  }

  # The tuples both remaining styles draw: the consistent ones, or every row
  # when the lattice has already been filtered (e.g. `sd_region_data(rule =
  # "integer")`) and carries no verdict column
  pts <- if ("consistent" %in% names(umbrella)) {
    umbrella[
      !is.na(umbrella$consistent) & umbrella$consistent,
      c("mean", "sd"),
      drop = FALSE
    ]
  } else {
    umbrella[, c("mean", "sd"), drop = FALSE]
  }

  lo_m <- min(umbrella$mean)
  hi_m <- max(umbrella$mean)
  pad <- expand * (hi_m - lo_m)
  y_hi <- max(c(pts$sd, curve$max_sd), na.rm = TRUE)

  p <- ggplot2::ggplot(pts, ggplot2::aes(.data$mean, .data$sd))
  if (shade == "outside") {
    # the panel is infeasible everywhere except at the points themselves
    p <- p +
      ggplot2::annotate(
        "rect",
        xmin = -Inf,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf,
        fill = "grey10",
        alpha = 0.12
      )
  }

  if (style == "contour") {
    # Rings rather than a ribbon, for the reason band_polygon() exists: where
    # the umbrella genuinely stops, the shading must stay, and a ring leaves it
    # there by construction. The knockout is white, as in `plot_sd_bounds()`;
    # with nothing to knock out of, the ring is left unfilled instead, so
    # whatever is underneath shows through.
    rings <- umbrella_contour(pts, by = contour_by)
    if (!is.null(rings)) {
      p <- p +
        ggplot2::geom_polygon(
          data = rings,
          ggplot2::aes(x = .data$mean, y = .data$y, group = .data$ring),
          inherit.aes = FALSE,
          fill = if (shade == "outside") "white" else NA,
          colour = line_colour,
          linewidth = 0.35
        )
    }
  } else {
    if (is.null(point_size)) {
      point_size <- if (nrow(pts) > 12000) {
        0.045
      } else if (nrow(pts) > 2000) {
        0.12
      } else {
        0.35
      }
    }
    p <- p +
      ggplot2::geom_point(
        colour = point_colour,
        size = point_size,
        shape = 16,
        na.rm = TRUE
      )
  }

  if (!is.null(curve)) {
    cur <- curve[curve$feasible & is.finite(curve$max_sd), ]
    p <- p +
      ggplot2::geom_line(
        data = cur,
        ggplot2::aes(.data$mean, .data$max_sd),
        inherit.aes = FALSE,
        colour = reference_colour,
        linetype = "dashed",
        linewidth = 0.35
      ) +
      ggplot2::geom_line(
        data = cur,
        ggplot2::aes(.data$mean, .data$min_sd),
        inherit.aes = FALSE,
        colour = reference_colour,
        linetype = "dashed",
        linewidth = 0.35
      )
  }
  p +
    ggplot2::coord_cartesian(
      xlim = c(lo_m - pad, hi_m + pad),
      ylim = c(-pad, y_hi + pad),
      expand = FALSE
    ) +
    ggplot2::labs(x = "Mean", y = "Sample standard deviation", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
