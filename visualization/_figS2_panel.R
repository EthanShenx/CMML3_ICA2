# Shared per-isoform panel renderer for Supplementary Figure S2.
# Each of 25..31_figS2_<PROT>.R sources this and supplies (PROT, PANEL_LABEL).
suppressPackageStartupMessages({
  library(ggplot2)
})
source("_palettes.R")

render_figS2_panel <- function(prot, panel_label, save_name) {
  d <- read.csv(file.path(DATA2PLOT, "figS2_data.csv"),
                fileEncoding = "UTF-8", check.names = FALSE)
  d$is_productive <- as.logical(d$is_productive)
  d$is_productive[is.na(d$is_productive)] <- FALSE
  sub <- d[d$protein == prot, ]
  sub <- sub[order(sub$DG), ]
  sub$y <- seq_len(nrow(sub))
  sub$category <- factor(sub$category, levels = CAT_ORDER)
  sub$bar_col  <- CAT_COLORS[as.character(sub$category)]

  n_p     <- sum(sub$is_productive, na.rm = TRUE)
  med_dg  <- median(sub$DG, na.rm = TRUE)
  title   <- sprintf("%s — %d/%d productive  (median ΔG = %.2f)",
                     prot, n_p, nrow(sub), med_dg)

  top3 <- head(sub, 3)
  stars <- sub[which(sub$is_productive), ]

  p <- ggplot(sub, aes(x = DG, y = y)) +
    geom_col(fill = sub$bar_col, width = 0.78) +
    geom_text(data = stars, aes(x = DG - 0.05, y = y),
              label = "★", colour = "#FFD700", size = 2.4,
              fontface = "plain", hjust = 1) +
    geom_text(data = top3, aes(x = DG + 0.05, y = y, label = short_name),
              hjust = 0, size = 2, colour = "#222") +
    geom_vline(xintercept = -7, linetype = "dashed", colour = "grey50",
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
  save_panel(p, save_name, width = 5.5, height = 6)
}
