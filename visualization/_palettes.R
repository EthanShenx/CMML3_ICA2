# Shared colour palettes and constants used by every panel script.
# Each panel script sources this file with `source("_palettes.R")` from the
# visualization/ directory. The script-level logic remains identical to the
# original Python 07_make_figures.py.

PROTEINS <- c("A8", "A9", "A10", "A11", "A12", "A13", "A14")

# RColorBrewer "Paired" palette for the seven isoforms
PROT_COLORS <- c(
  A8  = "#a6cee3", A9  = "#1f78b4", A10 = "#b2df8a",
  A11 = "#33a02c", A12 = "#fb9a99", A13 = "#e31a1c",
  A14 = "#fdbf6f"
)

CAT_ORDER <- c("Short-chain", "Medium-chain", "α,β-Unsaturated",
               "Aromatic", "Retinoid", "Semialdehyde", "Other")

CAT_COLORS <- c(
  "Short-chain"     = "#4CAF50",
  "Medium-chain"    = "#2196F3",
  "α,β-Unsaturated" = "#FF9800",
  "Aromatic"        = "#9C27B0",
  "Retinoid"        = "#F44336",
  "Semialdehyde"    = "#00BCD4",
  "Other"           = "#9E9E9E"
)

# Resolve paths relative to this file regardless of cwd
.this_dir   <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) getwd()
)
DATA2PLOT  <- normalizePath(file.path(.this_dir, "..", "data2plot"),
                            mustWork = FALSE)
FIG_OUT    <- normalizePath(file.path(.this_dir, "..", "figures",
                                       "panels_R"),
                            mustWork = FALSE)
dir.create(FIG_OUT, showWarnings = FALSE, recursive = TRUE)

format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) sprintf("p = %.2e", p) else sprintf("p = %.3f", p)
}

# Universal overlay applied last in save_panel(): all panels get black
# axes/ticks/tick labels and Arial typography.
.panel_overlay <- function() {
  ggplot2::theme(
    text          = ggplot2::element_text(family = "Arial", colour = "black",
                                           face = "plain"),
    axis.title    = ggplot2::element_text(family = "Arial", colour = "black",
                                           size = 12, face = "plain"),
    axis.title.x  = ggplot2::element_text(family = "Arial", colour = "black",
                                           size = 12, face = "plain",
                                           margin = ggplot2::margin(t = 6)),
    axis.title.y  = ggplot2::element_text(family = "Arial", colour = "black",
                                           size = 12, face = "plain",
                                           margin = ggplot2::margin(r = 6)),
    axis.text     = ggplot2::element_text(family = "Arial", colour = "black",
                                           size = 10, face = "plain"),
    axis.line     = ggplot2::element_line(colour = "black", linewidth = 0.4),
    axis.ticks    = ggplot2::element_line(colour = "black", linewidth = 0.4),
    plot.title    = ggplot2::element_text(family = "Arial", colour = "black",
                                           face = "plain"),
    plot.subtitle = ggplot2::element_text(family = "Arial", colour = "black",
                                           face = "plain"),
    legend.title  = ggplot2::element_text(family = "Arial", colour = "black",
                                           face = "plain"),
    legend.text   = ggplot2::element_text(family = "Arial", colour = "black",
                                           face = "plain"),
    strip.text    = ggplot2::element_text(family = "Arial", colour = "black",
                                           face = "plain")
  )
}

save_panel <- function(plot_obj, name, width = 5, height = 4, dpi = 300) {
  out <- file.path(FIG_OUT, paste0(name, ".png"))
  # Skip the overlay only when the caller explicitly passes a non-ggplot
  # object (e.g. tableGrob in 36_tableS1).
  if (inherits(plot_obj, "ggplot")) {
    plot_obj <- plot_obj + .panel_overlay()
  }
  # ragg::agg_png resolves system fonts (including Arial) reliably on macOS
  dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggplot2::ggsave(out, plot_obj, width = width, height = height,
                  dpi = dpi, bg = "white", device = dev)
  message("✓ ", basename(out))
}

# Convenience preview for RStudio: `panel_view(p)` renders the plot with the
# same overlay that save_panel() applies, so what you see line-by-line in
# RStudio matches what gets written to PNG.
panel_view <- function(p) {
  if (inherits(p, "ggplot")) print(p + .panel_overlay()) else print(p)
  invisible(p)
}

# Auto-apply the overlay to plot previews in interactive sessions
# (RStudio, R REPL). This is a no-op when run under Rscript/run_all.R.
if (interactive() && requireNamespace("ggplot2", quietly = TRUE)) {
  .orig_print_ggplot <- tryCatch(
    getS3method("print", "ggplot", envir = asNamespace("ggplot2")),
    error = function(e) NULL
  )
  if (!is.null(.orig_print_ggplot)) {
    print.ggplot <- function(x, ...) {
      if (inherits(x, "ggplot") &&
          !isTRUE(attr(x, ".panel_overlay_applied"))) {
        x <- x + .panel_overlay()
        attr(x, ".panel_overlay_applied") <- TRUE
      }
      .orig_print_ggplot(x, ...)
    }
    registerS3method("print", "ggplot", print.ggplot,
                     envir = globalenv())
  }
}
