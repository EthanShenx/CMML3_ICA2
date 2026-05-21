# Compose Supplementary Figure S2 (8 panels in 4x2 layout) into a single PDF.
# Panels (a)-(g): per-isoform ranked ΔG bar charts (A8..A14).
# Panel  (h):    inter-isoform Pearson ΔG correlation heatmap.
# Output: ../../ShenYuchen_overleaf_repo/figures/figS2_isoforms.pdf
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})
source("_palettes.R")

# ── Build a per-isoform panel as a ggplot (no save) ──
build_panel <- function(prot, panel_label) {
  d <- read.csv(file.path(DATA2PLOT, "figS2_data.csv"),
                fileEncoding = "UTF-8", check.names = FALSE)
  d$is_productive <- as.logical(d$is_productive)
  d$is_productive[is.na(d$is_productive)] <- FALSE
  sub <- d[d$protein == prot, ]
  sub <- sub[order(sub$DG), ]
  sub$y <- seq_len(nrow(sub))
  sub$category <- factor(sub$category, levels = CAT_ORDER)
  sub$bar_col  <- CAT_COLORS[as.character(sub$category)]

  med_dg  <- median(sub$DG, na.rm = TRUE)
  title   <- sprintf("%s  (median ΔG = %.2f)", prot, med_dg)

  top3 <- head(sub, 3)
  stars <- sub[which(sub$is_productive), ]

  ggplot(sub, aes(x = DG, y = y)) +
    geom_col(fill = sub$bar_col, width = 0.78) +
    geom_text(data = stars, aes(x = DG - 0.05, y = y),
              label = "★", colour = "#FFD700", size = 2.4,
              fontface = "plain", hjust = 1) +
    geom_text(data = top3, aes(x = DG + 0.05, y = y, label = short_name),
              hjust = 0, size = 2, colour = "#222") +
    geom_vline(xintercept = -6, linetype = "dashed", colour = "grey50",
               linewidth = 0.4) +
    labs(x = expression(Delta * G~"(kcal mol"^{-1}*")"),
         y = NULL, title = title, subtitle = panel_label) +
    theme_classic(base_size = 8) +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          plot.subtitle = element_text(face = "plain", hjust = -0.05),
          plot.title    = element_text(size = 8,
                                        colour = PROT_COLORS[[prot]],
                                        face = "plain", hjust = 0.5))
}

# ── Panel (h): inter-isoform Pearson correlation heatmap ──
build_panel_h <- function() {
  m <- read.csv(file.path(DATA2PLOT, "fig1h_corr.csv"), check.names = FALSE)
  names(m)[1] <- "row"
  long <- tidyr::pivot_longer(m, -row, names_to = "col", values_to = "r")
  long$row <- factor(long$row, levels = PROTEINS)
  long$col <- factor(long$col, levels = PROTEINS)
  ggplot(long, aes(x = col, y = row, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 2.6, colour = "#222") +
    scale_fill_distiller(palette = "RdYlBu", direction = -1,
                         limits = c(0, 1), name = "Pearson r") +
    scale_y_discrete(limits = rev(levels(long$row))) +
    coord_fixed() +
    labs(x = NULL, y = NULL,
         title = "Inter-isoform correlation",
         subtitle = "(h)") +
    theme_classic(base_size = 8) +
    theme(plot.subtitle = element_text(face = "plain", hjust = -0.05),
          plot.title    = element_text(size = 8, hjust = 0.5,
                                        colour = "#333"),
          axis.text     = element_text(size = 7),
          legend.key.height = grid::unit(0.6, "cm"),
          legend.key.width  = grid::unit(0.25, "cm"))
}

# ── Build all 8 panels ──
panels <- list(
  build_panel("A8",  "(a)"),
  build_panel("A9",  "(b)"),
  build_panel("A10", "(c)"),
  build_panel("A11", "(d)"),
  build_panel("A12", "(e)"),
  build_panel("A13", "(f)"),
  build_panel("A14", "(g)"),
  build_panel_h()
)

# ── Compose 4 rows x 2 cols ──
combined <- (panels[[1]] | panels[[2]]) /
            (panels[[3]] | panels[[4]]) /
            (panels[[5]] | panels[[6]]) /
            (panels[[7]] | panels[[8]])

# Save to the overleaf repo's figures dir
out_pdf <- normalizePath(
  file.path(.this_dir, "..", "..",
            "ShenYuchen_overleaf_repo", "figures",
            "figS2_isoforms.pdf"),
  mustWork = FALSE)

ggsave(out_pdf, combined, width = 13, height = 16, units = "in",
       device = cairo_pdf, bg = "white")
message("✓ Wrote ", out_pdf)
