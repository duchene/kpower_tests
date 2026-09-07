# kpower simulation tests -- shared publication figure theme
#
# Sourced by 03_-06_. Defines the palettes, theme, scales and save helpers so
# every figure in the manuscript is typographically and chromatically
# consistent. Reads and writes nothing.
#
# Palette provenance: the family and class palettes are Okabe-Ito steps chosen
# so that EVERY pair separates under normal vision and under simulated
# deuteranopia, protanopia and tritanopia (OKLab dE >= 8, normal-vision floor
# 15), and every swatch clears 3:1 contrast on white. The previous palette
# (#d6604d/#4393c3/#2ca02c) failed both checks: +R vs +T had dE 2.4 under
# deuteranopia -- indistinguishable to ~6% of male readers -- and +H vs +T had
# dE 3.8 under tritanopia.
#
# Environment variables:
#   KPOWER_FIG_TITLES=1  draw in-figure titles (default: off; captions belong
#                        in the manuscript). The multi-page diagnostic atlases
#                        (05, 06) always title their pages regardless.
#   KPOWER_FIG_DPI=<n>   raster resolution for the PNG copies (default 600).

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

# --- Palettes ---------------------------------------------------------------
# Categorical: assigned in fixed order, never cycled.
PALETTE_FAMILY <- c("+R" = "#D55E00",   # vermillion
                    "+H" = "#0072B2",   # blue
                    "+T" = "#009E73")   # bluish green
PALETTE_CLASS  <- c("sep" = "#0072B2",  # blue
                    "close" = "#B35900")# dark orange

# Sequential: ONE hue, light -> dark. Single-hue keeps the heatmap readable in
# greyscale reproduction, which a diverging or rainbow ramp does not. ColorBrewer
# Blues: the conventional scheme for a proportion, and neutral -- a red/orange
# ramp reads as heat/alarm on a quantity that is *good* when high, and would sit
# in the same hue family as vermillion "+R". The low end is Blues step 2 rather
# than the near-white step 1, so 0% cells still read as tiles.
SEQ_LOW  <- "#DEEBF7"
SEQ_HIGH <- "#08519C"

# Text and structural ink. Text never wears the series colour.
INK <- c(primary = "#1A1A1A", secondary = "#4D4D4D",
         muted   = "#767676", grid      = "#E4E4E4")

# --- Figure geometry --------------------------------------------------------
# Journal column widths in inches (88 / 140 / 183 mm).
FIG_WIDTH <- c(single = 3.46, onehalf = 5.51, double = 7.20)

fig_width <- function(w) if (is.character(w)) unname(FIG_WIDTH[[w]]) else w

# --- Theme ------------------------------------------------------------------
# Sized for final print size: no post-hoc scaling, so 8 pt here is 8 pt on the
# page. Recessive grid, no panel border, no grey strip boxes.
theme_kpower <- function(base_size = 8, grid = c("y", "both", "none")) {
  grid <- match.arg(grid)
  th <- theme_minimal(base_size = base_size) +
    theme(
      text             = element_text(colour = INK[["primary"]]),
      plot.title       = element_text(size = base_size + 1, face = "bold",
                                      hjust = 0, margin = margin(b = 3)),
      plot.subtitle    = element_text(size = base_size - 0.5,
                                      colour = INK[["secondary"]],
                                      margin = margin(b = 5)),
      plot.caption     = element_text(size = base_size - 1,
                                      colour = INK[["muted"]], hjust = 0),
      axis.title       = element_text(size = base_size,
                                      colour = INK[["secondary"]]),
      axis.text        = element_text(size = base_size - 1,
                                      colour = INK[["secondary"]]),
      axis.ticks       = element_line(colour = INK[["grid"]], linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = INK[["grid"]], linewidth = 0.3),
      panel.spacing    = unit(7, "pt"),
      strip.text       = element_text(size = base_size, face = "bold",
                                      colour = INK[["primary"]],
                                      margin = margin(2, 2, 2, 2)),
      strip.background = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_text(size = base_size - 1,
                                      colour = INK[["secondary"]]),
      legend.text      = element_text(size = base_size - 1),
      legend.key.height = unit(8, "pt"),
      legend.key.width  = unit(14, "pt"),
      legend.margin    = margin(t = 0),
      plot.margin      = margin(4, 6, 4, 4)
    )
  if (grid == "y")    th <- th + theme(panel.grid.major.x = element_blank())
  if (grid == "none") th <- th + theme(panel.grid.major = element_blank())
  th
}

# --- Shared scales ----------------------------------------------------------
scale_x_length <- function(...) {
  scale_x_log10(breaks = c(100, 300, 1000, 3000, 6000),
                labels = c("100", "300", "1000", "3000", "6000"), ...)
}

# Proportions are labelled as percentages: readers compare percentages faster
# than they compare 0.25 / 0.50 / 0.75.
scale_y_prop <- function(limits = c(0, 1), breaks = seq(0, 1, 0.25), ...) {
  scale_y_continuous(limits = limits, breaks = breaks,
                     labels = label_percent(accuracy = 1), expand = c(0.01, 0),
                     ...)
}

# --- Title control ----------------------------------------------------------
# Manuscript figures carry their explanation in the caption, not burned into
# the artwork. Set KPOWER_FIG_TITLES=1 to draw them while working.
FIG_TITLES <- Sys.getenv("KPOWER_FIG_TITLES", "0") == "1"

fig_labs <- function(..., title = NULL, subtitle = NULL) {
  if (!FIG_TITLES) { title <- NULL; subtitle <- NULL }
  labs(..., title = title, subtitle = subtitle)
}

# --- Saving -----------------------------------------------------------------
# Every single-page figure is written three ways from one plot object, at the
# same physical size, so the three copies are the same figure and not three
# re-layouts of it:
#   .pdf  cairo_pdf, which embeds fonts (the default pdf() device does not, and
#         unembedded fonts are the most common reason a journal bounces a
#         vector figure)
#   .svg  vector, for editing in Illustrator / Inkscape / a browser
#   .png  raster at FIG_PNG_DPI, for slides and for previewing in tools that
#         will not render a PDF
FIG_PNG_DPI <- as.numeric(Sys.getenv("KPOWER_FIG_DPI", "600"))

save_fig <- function(p, path, width = "single", height = 3.0) {
  w <- fig_width(width)
  ggsave(path, p, width = w, height = height, units = "in",
         device = grDevices::cairo_pdf)
  # svglite is not installed here; grDevices::svg() is the cairo SVG device.
  ggsave(sub("\\.pdf$", ".svg", path), p, width = w, height = height,
         units = "in", device = grDevices::svg, bg = "white")
  ggsave(sub("\\.pdf$", ".png", path), p, width = w, height = height,
         units = "in", dpi = FIG_PNG_DPI, bg = "white")
  message("  ", sub("\\.pdf$", ".{pdf,svg,png}", basename(path)),
          sprintf("  (%.2f x %.2f in)", w, height))
  invisible(path)
}

# Multi-page diagnostic atlases (05, 06). These stay PDF-only: SVG and PNG are
# one image per file, so there is no multi-page equivalent to export.
open_atlas <- function(path, width = 6.5, height = 4.5) {
  grDevices::cairo_pdf(path, width = width, height = height, onefile = TRUE)
}

# Readable label colour for a filled tile: flip to white once the fill is dark
# enough that near-black text would fail contrast.
label_ink_for_fill <- function(values, low = SEQ_LOW, high = SEQ_HIGH,
                               limits = c(0, 1)) {
  ramp <- grDevices::colorRamp(c(low, high))
  frac <- (values - limits[1]) / diff(limits)
  frac[is.na(frac)] <- 0
  rgb  <- ramp(pmin(pmax(frac, 0), 1)) / 255
  lin  <- ifelse(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
  lum  <- lin %*% c(0.2126, 0.7152, 0.0722)
  ifelse(lum < 0.42, "#FFFFFF", INK[["primary"]])
}
