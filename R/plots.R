# ggplot2 visualisations of the SD bounds, POMP-normalised bounds, and umbrella.

# Internal: shade the whole panel, then knock the feasible region out of it.
# The shared construction behind plot_sd_bounds(), plot_sd_region() and
# plot_umbrella(shade = "outside").
#
# Doing it this way rather than assembling the grey from side rectangles and
# ribbons means a mean at which the region is undefined stays shaded by
# construction: there is nothing there to knock out. Assembling the grey
# instead leaves such a gap unshaded, which would assert that any SD at all is
# possible there. See band_polygon().
#
# `rings` is band_polygon()/umbrella_contour() output, empty or NULL when
# nothing is feasible. With `shade = FALSE` there is no grey to knock a hole
# in, so the rings are left unfilled and whatever sits underneath shows
# through; `colour` strokes their outline and NA leaves them unstroked.
knockout_layers <- function(
  rings,
  shade = TRUE,
  colour = NA,
  linewidth = 0.35
) {
  c(
    if (shade) {
      list(ggplot2::annotate(
        "rect",
        xmin = -Inf,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf,
        fill = "grey10",
        alpha = 0.12
      ))
    },
    if (!is.null(rings) && nrow(rings)) {
      args <- list(
        data = rings,
        mapping = ggplot2::aes(
          x = .data$mean,
          y = .data$y,
          group = .data$ring
        ),
        inherit.aes = FALSE,
        fill = if (shade) "white" else NA
      )
      # A plain knockout draws no outline, so it leaves linewidth at the geom
      # default rather than setting one nothing renders. Only the stroked
      # contour needs it.
      if (!is.na(colour)) {
        args <- c(args, list(colour = colour, linewidth = linewidth))
      }
      list(rlang::exec(ggplot2::geom_polygon, !!!args))
    }
  )
}

# Internal: the plotted window, padded by a proportion of the scale width
# rather than a fixed number of SD units, so the margin looks the same on a
# 1-5 scale and a 0-100 one. Limits are set explicitly rather than left to
# ggplot2 because the outside shading needs finite ones.
padded_coord <- function(lo, hi, y_hi, expand) {
  pad <- expand * (hi - lo)

  ggplot2::coord_cartesian(
    xlim = c(lo - pad, hi + pad),
    ylim = c(-pad, y_hi + pad),
    expand = FALSE
  )
}

# Internal: the shared reported-point layer, green/red outlined dots.
bounds_point_layer <- function(points, x_var, y_var) {
  if (is.null(points[["consistent"]])) {
    points$consistent <- NA
  }
  points$consistent_label <- as.character(points$consistent)

  list(
    ggplot2::geom_point(
      data = points,
      ggplot2::aes(
        x = .data[[x_var]],
        y = .data[[y_var]],
        fill = .data$consistent_label
      ),
      shape = 21,
      colour = "black",
      size = 2.4,
      na.rm = TRUE
    ),
    ggplot2::scale_fill_manual(
      values = c("TRUE" = COLOUR_CONSISTENT, "FALSE" = COLOUR_INCONSISTENT),
      na.value = "grey50",
      name = "Consistent",
      labels = c("TRUE" = "consistent", "FALSE" = "inconsistent")
    )
  )
}

#' Plot SD bounds on the native scale
#'
#' The feasible SD band (floor to ceiling) against the mean, optionally with
#' reported points coloured by consistency.
#'
#' `shade = "outside"` (default) shades the infeasible region and leaves the
#' feasible one clear, matching [plot_sd_bounds_pomp()] with
#' `reference = "sharp"`, so the two scales read the same way: shaded means
#' unreachable. `shade = "inside"` fills the feasible band instead.
#'
#' @param curve Output of [sd_bounds_curve()]. Its `"step"` attribute, when
#'   present, is the mean-grid spacing used to tell a sampling gap from a
#'   genuine one under `shade = "outside"`; a curve built by hand, without that
#'   attribute, should be on a uniform grid, from which the spacing is inferred.
#' @param points Optional data frame with `mean`, `sd`, and optionally
#'   `consistent`; e.g. the output of [brimmer_map()].
#' @param title Optional plot title.
#' @param fill,line_colour Band fill and outline colours. `fill` is used only
#'   by `shade = "inside"`; the infeasible shading has its own fixed grey.
#' @param shade `"outside"` (default) shades the infeasible region;
#'   `"inside"` fills the feasible band.
#' @param expand Padding around the plotted region, as a proportion of the
#'   scale width. The limits always stretch to include `points`, so an
#'   out-of-bounds report is never clipped out of view.
#' @return A ggplot object.
#' @examples
#' curve <- sd_bounds_curve(l = 1, u = 7, n = 30, by = 0.1)
#' plot_sd_bounds(curve, title = "Feasible SDs, 1-7 scale, n = 30")
#'
#' # overlay reported values, coloured by consistency (the second is
#' # above the ceiling, so it plots as inconsistent)
#' reports <- tibble::tibble(mean = c(2.97, 3.51), sd = c(2.83, 3.50))
#' checked <- brimmer_map(
#'   reports,
#'   l = 1, u = 7, n = 30, digits_mean = 2, digits_sd = 2
#' )
#' plot_sd_bounds(curve, points = checked)
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
  shade <- rlang::arg_match(shade)
  feasible <- curve[curve$feasible & is.finite(curve$max_sd), ]

  # Finite limits will silently clip an out-of-bounds point - exactly the case
  # the plot exists to show - so the ceiling of the view must account for the
  # reported points as well.
  y_hi <- max(c(feasible$max_sd, points$sd), na.rm = TRUE)

  plot <- ggplot2::ggplot(feasible, ggplot2::aes(x = .data$mean))

  plot <- plot +
    if (shade == "inside") {
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$min_sd, ymax = .data$max_sd),
        fill = fill
      )
    } else {
      # band_polygon() needs the grid spacing to tell a sampling gap from a
      # genuine one, and sd_bounds_curve() records the spacing it used because
      # its grid is not uniform: it adds each kink of the 1/(n * n_items)
      # lattice plus a pair of neighbours a tolerance away. The median of the
      # realised spacing is the fallback for a hand-built curve, where a
      # uniform grid makes it right.
      step <- attr(curve, "step")
      knockout_layers(band_polygon(
        tibble::tibble(
          mean = feasible$mean,
          lo = feasible$min_sd,
          hi = feasible$max_sd
        ),
        by = if (!is.null(step)) {
          step
        } else if (nrow(feasible) > 1) {
          stats::median(diff(feasible$mean))
        } else {
          1
        }
      ))
    }

  plot <- plot +
    ggplot2::geom_line(ggplot2::aes(y = .data$max_sd), colour = line_colour) +
    ggplot2::geom_line(ggplot2::aes(y = .data$min_sd), colour = line_colour) +
    padded_coord(min(feasible$mean), max(feasible$mean), y_hi, expand) +
    ggplot2::labs(x = "Mean", y = "SD", title = title) +
    ggplot2::theme_minimal()

  if (!is.null(points)) {
    plot <- plot + bounds_point_layer(points, "mean", "sd")
  }

  plot
}

#' Plot SD bounds on a percent-of-maximum-possible (POMP) scale
#'
#' `reference = "parity"` normalises every SD by the mean-agnostic parity
#' (Popoviciu) ceiling: a linear rescaling, so the Structure S ceiling appears
#' as a dome under 1 and the umbrella geometry is undistorted.
#' `reference = "sharp"` normalises each SD by its own sharp mean-conditional
#' band, so the feasible region is exactly the unit square and a point's height
#' is its position within the band (`pomp_sd_sharp`); regions outside
#' `[0, 1]^2` are shaded infeasible.
#'
#' @param curve Output of [sd_bounds_curve()] (used for the parity band).
#' @param points Optional data frame with `pomp_mean` and `pomp_sd_parity` /
#'   `pomp_sd_sharp` and `consistent`, e.g. from [brimmer_map()].
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
  reference <- rlang::arg_match(reference)

  if (reference == "parity") {
    feasible <- curve[curve$feasible & is.finite(curve$max_sd), ]
    plot <- ggplot2::ggplot(feasible, ggplot2::aes(x = .data$pomp_mean)) +
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
      plot <- plot + bounds_point_layer(points, "pomp_mean", "pomp_sd_parity")
    }

    return(plot + ggplot2::theme_minimal())
  }

  outside <- tibble::tibble(
    xmin = c(-Inf, 1, 0, 0),
    xmax = c(0, Inf, 1, 1),
    ymin = c(-Inf, -Inf, 1, -Inf),
    ymax = c(Inf, Inf, Inf, 0)
  )
  plot <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = outside,
      ggplot2::aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        ymin = .data$ymin,
        ymax = .data$ymax
      ),
      fill = "grey10",
      alpha = 0.12
    ) +
    ggplot2::annotate(
      "rect",
      xmin = 0,
      xmax = 1,
      ymin = 0,
      ymax = 1,
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
    plot <- plot + bounds_point_layer(points, "pomp_mean", "pomp_sd_sharp")
  }

  plot + ggplot2::theme_minimal()
}

#' Plot the umbrella grid
#'
#' Renders the reportable `(mean, sd)` tuples of a design, in any of three
#' styles. Optionally overlays the bound curves from [sd_bounds_curve()].
#'
#' `style = "points"` (default) greys the whole panel and draws only the
#' consistent tuples. Nothing else is drawn, because nothing else exists: every
#' other cell of the grid is a value that cannot be reported. This matches the
#' convention of [plot_sd_bounds()] and [plot_sd_region()], where shading marks
#' what is ruled out, and it makes the vertical striping legible.
#'
#' `style = "tiles"` draws every cell of the reporting grid, coloured as
#' consistent, GRIMMER-inconsistent, or out of bounds. Useful for methods
#' exposition, since it separates what the bounds rule out from what GRIMMER
#' additionally rules out, but it spends most of its ink on impossible tuples.
#'
#' `style = "contour"` draws only the outline of the umbrella, from
#' [umbrella_contour()]. Use it when the point cloud is too dense to read, when
#' the figure is a backdrop for reported points, or to see the discrete
#' envelope against the continuous one by passing `curve` as well. Being an
#' envelope, it is drawn across the vertical striping rather than around it, so
#' its interior claims less than the points do.
#'
#' The grey under `"points"` and `"contour"` is a layer, not a theme element:
#' it is a rectangle over the whole panel with the feasible region knocked out
#' of it, so it says "no tuple here" as data rather than as decoration. No
#' theme call removes it — `theme_void()` strips the axes and the shading
#' stays, as it must. Use `shade = "none"` to drop it.
#'
#' Note what the points represent: the GRIM- and GRIMMER-consistent set, which
#' is what [brimmer()] applies. It is strictly larger than the set of
#' attainable tuples, so an umbrella plot shows what the test admits rather
#' than what exists; use [brimmest()] to certify a tuple.
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
#'   that tells a genuine gap in the umbrella from its ordinary striping. Give
#'   it explicitly to keep a narrow gap from being bridged.
#' @param shade `"outside"` (default) greys everything the design rules out;
#'   `"none"` draws no shading at all, and leaves the contour unfilled. Applies
#'   to `"points"` and `"contour"`; `"tiles"` colours every cell already.
#' @param expand Padding as a proportion of the scale width, as in
#'   [plot_sd_bounds()].
#' @return A ggplot object.
#' @examples
#' grid <- umbrella_data(n = 12, l = 1, u = 3, digits = 2)
#' plot_umbrella(grid, title = "n = 12, 1-3 scale")
#'
#' # separating the two ways a tuple can be ruled out
#' plot_umbrella(grid, style = "tiles")
#'
#' # just the outline, with the continuous bounds dashed over it
#' curve <- sd_bounds_curve(l = 1, u = 3, n = 12, by = 0.05)
#' plot_umbrella(grid, curve = curve, style = "contour")
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
  style <- rlang::arg_match(style)
  shade <- rlang::arg_match(shade)

  bound_curves <- function(colour, linetype, linewidth) {
    if (is.null(curve)) {
      return(NULL)
    }
    feasible <- curve[curve$feasible & is.finite(curve$max_sd), ]
    list(
      ggplot2::geom_line(
        data = feasible,
        ggplot2::aes(.data$mean, .data$max_sd),
        inherit.aes = FALSE,
        colour = colour,
        linetype = linetype,
        linewidth = linewidth
      ),
      ggplot2::geom_line(
        data = feasible,
        ggplot2::aes(.data$mean, .data$min_sd),
        inherit.aes = FALSE,
        colour = colour,
        linetype = linetype,
        linewidth = linewidth
      )
    )
  }

  if (style == "tiles") {
    umbrella <- umbrella |>
      dplyr::mutate(
        category = dplyr::case_when(
          .data$consistent ~ "consistent",
          .data$in_bounds & !is.na(.data$grimmer) & !.data$grimmer ~
            "GRIMMER-inconsistent",
          .default = "out of bounds"
        )
      )

    return(
      ggplot2::ggplot(
        umbrella,
        ggplot2::aes(x = .data$mean, y = .data$sd, fill = .data$category)
      ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_manual(
          values = c(
            "consistent" = COLOUR_CONSISTENT,
            "GRIMMER-inconsistent" = COLOUR_GRIMMER,
            "out of bounds" = "grey80"
          ),
          name = NULL
        ) +
        bound_curves("grey20", "solid", 0.5) +
        ggplot2::labs(x = "Mean", y = "SD", title = title) +
        ggplot2::theme_minimal()
    )
  }

  # The tuples both remaining styles draw: the consistent ones, or every row
  # when the lattice has already been filtered and carries no verdict column
  points <- if ("consistent" %in% names(umbrella)) {
    umbrella[
      !is.na(umbrella$consistent) & umbrella$consistent,
      c("mean", "sd"),
      drop = FALSE
    ]
  } else {
    umbrella[, c("mean", "sd"), drop = FALSE]
  }

  plot <- ggplot2::ggplot(points, ggplot2::aes(.data$mean, .data$sd))

  if (style == "contour") {
    # Rings rather than a ribbon, for the reason band_polygon() exists: where
    # the umbrella genuinely stops, the shading must stay, and a ring leaves it
    # there by construction.
    plot <- plot +
      knockout_layers(
        umbrella_contour(points, by = contour_by),
        shade = shade == "outside",
        colour = line_colour
      )
  } else {
    # the panel is infeasible everywhere except at the points themselves, so
    # there is nothing to knock out of the shading
    if (shade == "outside") {
      plot <- plot + knockout_layers(NULL)
    }
    if (is.null(point_size)) {
      point_size <- if (nrow(points) > 12000) {
        0.045
      } else if (nrow(points) > 2000) {
        0.12
      } else {
        0.35
      }
    }
    plot <- plot +
      ggplot2::geom_point(
        colour = point_colour,
        size = point_size,
        shape = 16,
        na.rm = TRUE
      )
  }

  plot +
    bound_curves(reference_colour, "dashed", 0.35) +
    padded_coord(
      min(umbrella$mean),
      max(umbrella$mean),
      max(c(points$sd, curve$max_sd), na.rm = TRUE),
      expand
    ) +
    ggplot2::labs(x = "Mean", y = "Sample standard deviation", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
