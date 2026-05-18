# Fig 2 panel (h) — top-15 dual-winner pairs (ΔG ≤ −7 AND Cys-SG ≤ 4 Å) lollipop.
suppressPackageStartupMessages({
  library(ggplot2)
})
source("_palettes.R")

d <- read.csv(file.path(DATA2PLOT, "fig2h_dual_lollipop.csv"),
              fileEncoding = "UTF-8", check.names = FALSE)
d$protein <- factor(d$protein, levels = PROTEINS)
d$label   <- factor(d$label, levels = rev(d$label))
d$dot_size <- 120 / (d$distance + 0.5)

p <- ggplot(d, aes(x = DG, y = label, colour = protein)) +
  geom_segment(aes(x = DG, xend = -6.8,
                   y = label, yend = label),
               linewidth = 0.8, alpha = 0.8) +
  geom_point(aes(size = dot_size), shape = 21, fill = NA, stroke = 0,
             show.legend = FALSE) +
  geom_point(aes(size = dot_size), shape = 21, fill = NULL,
             colour = "black", stroke = 0.4, show.legend = FALSE) +
  geom_point(aes(size = dot_size, fill = protein), shape = 21,
             colour = "black", stroke = 0.4) +
  geom_vline(xintercept = -7, linetype = "dashed", colour = "grey50",
             linewidth = 0.3) +
  scale_x_reverse() +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  scale_colour_manual(values = PROT_COLORS, guide = "none") +
  scale_fill_manual(values = PROT_COLORS, name = NULL) +
  labs(x = expression(Delta * G~"(kcal mol"^{-1}*")"), y = NULL,
       title = paste("Top-15 dual-winner pairs (ΔG ≤ −7 kcal mol⁻¹ AND Cys-SG ≤ 4 Å)\n",
                     "Dot size ∝ proximity (1/distance)"),
       subtitle = "(h)") +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 4))) +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size = 10),
        plot.subtitle = element_text(face = "plain", hjust = -0.02),
        plot.title    = element_text(size = 8),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 7))
save_panel(p, "18_fig2h_dual_lollipop", width = 9, height = 5)
